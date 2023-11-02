#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Add alignment based tags to the input vcf file.

Created at Monday 23 October 2023  12:54 by Kimmo Palin <kpalin@helsinki.fi>
"""

__version__ = "0.1"

from argparse import Namespace
import logging
from os import EX_OK
import sys
from typing import List, Optional, Tuple, Union
import pysam


def parse_args(args: List[str]) -> Namespace:
    import argparse

    description = "\n".join(__doc__.splitlines()[:-1]).strip()

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-v", "--input", help="Input vcf file [default:%(default)s]", required=True
    )

    parser.add_argument(
        "-a",
        "--alignment",
        help="Input file (aligned cram) [default:%(default)s]",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output vcf file  [default:%(default)s]",
        default="output.vcf",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use  [default:%(default)s]",
        default=7,
        type=int,
    )
    version = f"%(prog)s {__version__}"
    parser.add_argument("--version", action="version", version=version)

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        action="store_true",
        help="Be more verbose with output",
    )

    args = parser.parse_args()

    import logging

    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s:%(funcName)s:%(levelname)s:%(message)s",
        )

    return args


from collections import namedtuple

ContigSplit = namedtuple("ContigSplit", ("CHROM", "START", "END"))


def split_contig_regions(chromosome_lengths, max_region_length) -> List[ContigSplit]:
    """Split chromosomes to parts for parallel processing

    Args:
        chromosome_lengths (_type_): _description_
        max_region_length (_type_): _description_

    Returns:
        List[ContigSplit]: _description_
    """
    max_region_length = int(max_region_length)
    contig_regions = []

    for chromosome, length in chromosome_lengths.items():
        start = 0
        for start in range(0, length, max_region_length):
            end = min(start + max_region_length, length)
            contig_regions.append(ContigSplit(chromosome, start, end))

    return contig_regions


def _add_MQ_header_info(vcf_file: pysam.VariantFile) -> pysam.VariantHeader:
    # Create a new header with INFO records for MQ and MQ0
    new_header = vcf_file.header
    new_header.info.add("MQ", ".", "Float", "RMS Mapping Quality")
    new_header.info.add("MQ0", 1, "Integer", "Number of reads with mapping quality 0")
    new_header.info.add("DP", 1, "Integer", "Total read depth at variant.")
    new_header.info.add("SOMATIC", 0, "Flag", "Somatic variant.")
    new_header.filters.add(
        "MapQ0", None, None, r"MappingQ zero reads more than 10% of variant reads"
    )
    return new_header


def main(argv: List[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)

    # Input and output VCF file paths
    input_vcf = args.input
    output_vcf = args.output
    input_bam = args.alignment
    n_threads = args.threads

    # Open the input VCF file for reading
    vcf_reader = pysam.VariantFile(input_vcf, "r")

    # Create a dictionary to store chromosome names and their corresponding lengths.
    chromosome_lengths = {}

    # Iterate through the records in the VCF file to collect chromosome information.
    for chromosome_name, chromosome_length in vcf_reader.header.contigs.items():
        chromosome_lengths[chromosome_name] = chromosome_length.length

    chromosome_splits = split_contig_regions(chromosome_lengths, 20e6)

    new_header = _add_MQ_header_info(vcf_reader)

    # new_header.info.add("RMSMappingQuality", ".", "Float", "Root Mean Squared mapping quality")
    vcf_reader.close()

    # Open the output VCF file for writing with the updated header
    vcf_writer = pysam.VariantFile(output_vcf, "w", header=new_header)

    from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
    from functools import partial
    import itertools as it
    import logging

    assert (
        len(vcf_writer.header.samples) == 1
    ), "Must have exactly one sample in the input VCF"
    tumor_name = vcf_writer.header.samples[0]
    read_group = get_read_group(input_bam)

    with ThreadPoolExecutor(n_threads) as pool:
        add_mq_chunk_part = partial(add_mq_chunk, input_bam, input_vcf, read_group)
        out_records = pool.map(add_mq_chunk_part, chromosome_splits, chunksize=1)
        logging.info("Sent map")
        # for record in tqdm(it.chain.from_iterable(out_records)):
        # vcf_writer.write(record)
        # for record_block in tqdm(
        #     out_records, total=len(chromosome_splits), unit="chunks"
        # ):

        for record_block in out_records:
            if len(record_block) > 0:
                logging.info(
                    "Got a block of records start around %s:%d",
                    record_block[0].chrom,
                    record_block[0].pos,
                )

            for record in record_block:
                call_precise(record)

                call_somatic_qual(record, tumor_name)
                call_filters(record, tumor_name)
                vcf_writer.write(record)

    # Close the output VCF file

    vcf_writer.close()

    return EX_OK


def call_precise(record: pysam.VariantRecord) -> pysam.VariantRecord:
    """Set PRECISE if POS and LEN confidence intervals are empty. Otherwise set IMPRECISE

    Args:
        record (pysam.VariantRecord): _description_
        tumor_name (str): _description_

    Returns:
        pysam.VariantRecord: _description_
    """
    cipos = record.info.get("CIPOS", [-1e-9, 1e9])
    cilen = record.info.get("CILEN", [-1e-9, 1e9])
    if cipos[0] == 0 and cipos[1] == 0 and cilen[0] == cilen[1]:
        record.info["IMPRECISE"] = 0
        record.info["PRECISE"] = 1
    else:
        record.info["IMPRECISE"] = 1
        record.info["PRECISE"] = 0
    return record


def call_somatic(
    record: pysam.VariantRecord, tumor_rname: str
) -> Union[None, pysam.VariantRecord]:
    """Set SOMATIC flag if RNAMES contain only reads with tumor_rname.

    Args:
        record (pysam.VariantRecord): _description_
        tumor_rname (str): regex used to identify reads belongging to tumor.

    Returns:
        pysam.VariantRecord: _description_
    """
    import re

    rnames = record.info.get("RNAMES", list())
    read_sources = [rname.split(":")[-1] for rname in rnames]
    tumor_rname_re = re.compile(tumor_rname)

    n_reads_total = len(rnames)

    n_reads_tumor = sum(1 for rsrc in read_sources if tumor_rname_re.match(rsrc))
    n_reads_normal = n_reads_total - n_reads_tumor

    if n_reads_tumor > 0 and n_reads_normal == 0:
        record.info["SOMATIC"] = True
    elif n_reads_tumor > 0 and n_reads_normal > 0:
        record.info["SOMATIC"] = False
    elif n_reads_tumor == 0 and n_reads_normal > 0:
        return None
    else:
        raise ValueError(record.to_string())
    return record


def call_somatic_qual(
    record: pysam.VariantRecord, tumor_name: str
) -> pysam.VariantRecord:
    import numpy as np
    from scipy.stats import poisson

    # logging.info("Somatic calling %s", record.id)
    prior_variant_per_read = 0.05  # 5%

    fmt = record.samples[tumor_name]
    DepthVariant = fmt["DV"]
    if DepthVariant is None:
        logging.info("No DV format on %s", record.id)
        return record
    else:
        DepthVariant = int(DepthVariant)
    DepthReference = fmt["DR"]
    if DepthReference is None:
        logging.info("No DR format on %s", record.id)
        return record
    else:
        DepthReference = int(DepthReference)
    depth_total = DepthReference + DepthVariant
    expected_variants = max(2.0, depth_total * prior_variant_per_read)

    poisson_dist = poisson(expected_variants)
    # Here's magic number to convert logsf() to phred scale (i.e. -10/ln(10.))
    # Add extra 1 read for insertions since they are easier than the others to call.
    phred_qual = -4.3429448190325175 * poisson_dist.logsf(
        DepthVariant - 1 + (1 if record.info.get("SVTYPE", "None") == "INS" else 0)
    )
    phred_qual = min(phred_qual, 60.0)
    record.qual = int(np.around(phred_qual))

    return record


def call_filters(record: pysam.VariantRecord, tumor_name: str) -> pysam.VariantRecord:
    """Update variant filters q5, MapQ0. Set new heterozygous genotype if SOMATIC and PASS.

    Args:
        record (pysam.VariantRecord): _description_
        tumor_name (str): _description_

    Returns:
        pysam.VariantRecord: _description_
    """
    record.filter.clear()
    if record.qual is None or record.qual <= 5.0:
        record.filter.add("q5")

    DepthVariant = record.samples[tumor_name].get("DV", 1)
    if record.info["MQ0"] >= 0.1 * DepthVariant:
        record.filter.add("MapQ0")
    if len(record.filter) == 0:
        record.filter.add("PASS")
        if record.info.get("SOMATIC", False):
            record.samples[tumor_name]["GT"] = (0, 1)

    return record


def get_read_group(input_bam: str) -> str:
    bam = pysam.AlignmentFile(input_bam, "rb")
    rg = bam.header["RG"]
    rg_name = rg[0]["ID"]
    return rg_name


import pysam


class BamMapqFetcher:
    "Enable querying of mapping qualitites of reads on given region of a bam/cram file"

    def __init__(self, bam_file: pysam.AlignmentFile, chrom: str, start: int, end: int):
        from ncls import NCLS
        import numpy as np

        reads_chunk = bam_file.fetch(chrom, start, end)
        mapqs = np.fromiter(
            (
                (
                    r.reference_start,
                    r.reference_end,
                    BamMapqFetcher.read_hasher(r) | r.mapping_quality,
                )
                for r in reads_chunk
            ),
            dtype=(np.int64, 3),
        )
        self.read_quals = NCLS(mapqs[:, 0], mapqs[:, 1], mapqs[:, 2])

    def read_hasher(r: pysam.AlignedSegment) -> int:
        """Set to zero the 8 bottom bits of the hash of query name of the aligned segment

        Args:
            r (pysam.AlignedSegment): The segment (read) to hash

        Returns:
            int: Hash with bottom 8 bits set to zero.
        """
        return hash(r.query_name) & (~0xFF)

    def get_mapqs(self, start, end):
        import numpy as np

        it = self.read_quals.find_overlap(start, end)
        mapqs = np.fromiter((x[2] for x in it), dtype=int)
        mapqs = {b & (~0xFF): b & 0xFF for b in mapqs}
        return mapqs


def add_mq_chunk(input_bam, input_vcf, read_group: str, csplit: Tuple[str, int, int]):
    import pysam
    import logging
    import time
    import itertools as it

    start_time = time.time()

    logging.debug("Running %s:%d-%d, %s", *csplit, read_group)
    import re

    out_records = []
    # Open the input VCF file for reading
    vcf_reader = pysam.VariantFile(input_vcf, "r")
    _add_MQ_header_info(vcf_reader)
    # Open the BAM file for reading
    bam = pysam.AlignmentFile(input_bam, "rb")
    mapq_fetch = None
    SLOP = 10
    contig_len = bam.header.get_reference_length(csplit.CHROM)
    # Iterate over each variant in the VCF
    for record in vcf_reader.fetch(csplit.CHROM, csplit.START, csplit.END):
        # Fetch reads at the location of the record from the BAM file
        if record.pos < csplit.START or record.pos >= csplit.END:
            continue
        record = call_somatic(record, read_group)
        if not record:
            continue

        if mapq_fetch is None:
            mapq_fetch = BamMapqFetcher(bam, csplit.CHROM, csplit.START, csplit.END)
        start_contig, start_pos = record.contig, record.pos
        ci_pos = record.info.get("CIPOS", [-SLOP, SLOP])
        ci_len = record.info.get("CILEN", [-SLOP, SLOP])

        # Treat break ends like long deletions (different ends separately)
        if record.info.get("SVTYPE", None) == "BND":
            m = re.match(r"[a-zA-Z]*[][]([^:]+):(\d+)[][][a-zA-Z]*", record.alts[0])
            if m is None:
                raise ValueError(
                    "Weird thing with alt of %s : %s" % (record.id, record.alts[0])
                )
            else:
                end_contig, end_pos = m.groups()
                end_pos = int(end_pos) - 1  # Zero based
        else:
            end_contig, end_pos = record.contig, record.stop

        if end_pos - start_pos < 1000 and start_contig == end_contig:
            min_start = max(0, start_pos + ci_pos[0])
            max_end = min(
                bam.header.get_reference_length(end_contig) - 1,
                end_pos + ci_pos[1] + ci_len[1],
            )
            # reads = bam.fetch(start_contig, min_start, max_end)
            # mapping_qualities = [read.mapping_quality for read in reads]
            mapping_qualities = mapq_fetch.get_mapqs(min_start, max_end)
        else:
            # For long variants (on reference) count the values for reads around breakpoints.
            logging.info(
                "Long variant: %s at %s:%d -  %s:%d",
                record.id,
                start_contig,
                start_pos,
                end_contig,
                end_pos,
            )

            # reads = bam.fetch(
            #     start_contig,
            #     max(0, start_pos + ci_pos[0]),
            #     min(contig_len - 1, start_pos + ci_pos[1]),
            # )
            # mapping_qualities = {
            #     read.query_name: read.mapping_quality for read in reads
            # }
            min_start, max_start = max(0, start_pos + ci_pos[0]), min(
                contig_len - 1, start_pos + ci_pos[1]
            )
            mapping_qualities = mapq_fetch.get_mapqs(min_start, max_start)

            # SV end mapping qualities
            min_end, max_end = max(0, end_pos + ci_pos[0] + ci_len[0]), min(
                bam.header.get_reference_length(end_contig) - 1,
                end_pos + ci_pos[1] + ci_len[1],
            )
            if (
                end_contig != csplit.CHROM
                or min_end < csplit.START
                or max_end > csplit.END
            ):
                reads = bam.fetch(
                    end_contig,
                    min_end,
                    max_end,
                )

                mapping_qualities.update(
                    {
                        BamMapqFetcher.read_hasher(read): read.mapping_quality
                        for read in reads
                    }
                )

                # mapping_qualities = list(mapping_qualities.values())
            else:
                mapping_qualities.update(mapq_fetch.get_mapqs(min_end, max_end))

        mapping_qualities = list(mapping_qualities.values())
        # Calculate Mean Mapping Quality (MQ) for each read
        mean_mapping_quality = sum(mapping_qualities) / len(mapping_qualities)

        # Calculate RMS Mapping Quality
        rms_mapping_quality = (
            sum(mq**2 for mq in mapping_qualities) / len(mapping_qualities)
        ) ** 0.5

        mapq_zero_reads = sum(1 for q in mapping_qualities if q == 0)
        # Add INFO tags to the record
        record.info["MQ"] = rms_mapping_quality
        record.info["DP"] = len(mapping_qualities)
        record.info["MQ0"] = mapq_zero_reads
        out_records.append(record)
    bam.close()
    vcf_reader.close()
    t = time.time() - start_time
    logging.debug("Done %s:%d-%d in %gs", *csplit, t)
    return out_records


if __name__ == "__main__":
    sys.exit(main())

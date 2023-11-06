import argparse
import sys
import logging
from typing import Dict, Iterable, List, Tuple, Union

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata


VERSION = metadata.version("cuddlySV")

from pathlib import Path


class WorkDir:
    def __init__(self,temporary_dir:str):
        self.temporary_dir = Path(temporary_dir)
        if not self.temporary_dir.is_dir:
          raise FileNotFoundError("[Errno 2] No such directory: '%s'" % temporary_dir)
        
    @property
    def path(self):
        return self.temporary_dir
        

    @property
    def idx(self)->Dict[str,Dict[Union[str,Tuple[str,str]],int]]:
        """Return starting file positions of chr or (chr1,chr2) in signature files. First level of
        indexing is the signature type

        Returns:
            Dict[svtype,Dict[Union[chrom,Tuple[chrom1,chrom2]],int]]: _description_
        """
        if not hasattr(self,"_index"):
            self._index = self.index_temp_files()
        return self._index
        
    def temp_dir_empty(self) -> bool:
        """Check if temporary directory is finished or needs update

        Args:
            temporary_dir (Path): Used work dir

        Returns:
            bool: True if any necessary file in the workdir is empty or missing
        """
        for kind in ["DEL", "DUP", "INS", "INV", "TRA", "reads"]:
            f = self.temporary_dir / f"{kind}.sigs"
            if not f.exists():
                f=f.parent / (f.name+".gz")
            try:
                if f.stat().st_size == 0:
                    return True
            except FileNotFoundError:
                return True
        return False
    
    def load_valuable_chr(self) -> Dict[str, Union[List[str], Dict[str, List[str]]]]:
        """Load a dictionary of signature types to list (or dict of lists) containing the chromosome
    names with that type of signatures.

    Args:
        path (str): The temporary work path for the signature files

    Returns:
        Dict[str,Union[List[str],Dict[str,List[str]]]]: From signature type to list of chromosome names.
        """
        v={}
        for svtype,chrd in self.idx.items():
            svtype =svtype.split(".")[0]
            if svtype=="TRA":
                v[svtype] = {}
                for chr1,chr2 in chrd.keys():
                    v[svtype].setdefault(chr1,list()).append(chr2)
                for _,chrom_list in v[svtype].items():
                    chrom_list.sort()
            else:
                v[svtype] = list(chrd.keys())
                
        return v
    
    def lines(self,svtype:str,chrom:str,chrom2=None)->Iterable[str]:
        """Generate lines from svtype signature collection for chrom (and chrom2 for TRA )

        Args:
            svtype (str): _description_
            chrom (str): _description_
            chrom2 (_type_, optional): _description_. Defaults to None.

        Returns:
            Iterable[str]: _description_

        Yields:
            Iterator[Iterable[str]]: _description_
        """
        fpath = self.path / f"{svtype}.sigs"
        with fpath.open("rt") as f:
            end_pos = sys.maxsize
            k = chrom if chrom2 is None else (chrom,chrom2)
            try:
                pos = self.idx[fpath.name][k]
                greater_end_pos = [x for x in self.idx[fpath.name].values() if x > pos]
                if len(greater_end_pos)>0:
                    end_pos = min(greater_end_pos)
                f.seek(pos)
            except KeyError as e:
                logging.warning("Couldn't find %s from %s",str(k),str(fpath))
            
            line = f.readline()
            while line!='' and f.tell()<=end_pos:
                yield line
                line = f.readline()
                
    
    def index_temp_files(self):
        """Find the start positions in the sig files

        Raises:
            IOError: _description_
        """
        idxs={}
        for sigfile in self.temporary_dir.glob("*.sigs"):
            if sigfile.name == "reads.sigs": continue

            idxs[sigfile.name] = {}
            prev_chrom = None
                
            logging.info("Starting to index %s",str(sigfile))
            with sigfile.open("rt") as inf:                    
                while True:                        
                    f_line = inf.readline()
                    if f_line=='':
                        break
                    
                    if sigfile.name in ("INS.sigs","DUP.sigs","DEL.sigs","INV.sigs"):
                        cur_chrom = f_line.split('\t')[1]
                    elif sigfile.name in ("TRA.sigs"):
                        parts = f_line.split('\t')
                        cur_chrom = parts[1],parts[4]
                        
                        
                    if prev_chrom != cur_chrom:
                        c_pos = inf.tell()- len(f_line)
                        idxs[sigfile.name][cur_chrom] = c_pos
                        prev_chrom = cur_chrom 
                            
        logging.info("Indexing done!")
        return idxs
    

class cuddlySVdp(object):
    """
    Detailed descriptions of cuddlySV version and its parameters.
    """

    USAGE = """\
		
	Current version: v%s
	Author: Tao Jiang, Kimmo Palin
	Contact: tjiang@hit.edu.cn (or google for Kimmo Palin)

	If you use cuteSV in your work, please cite:
		Jiang T et al. Long-read-based human genomic structural variation detection with cuteSV. 
		Genome Biol 21,189(2020). https://doi.org/10.1186/s13059-020-02107-y


	Suggestions:

	For PacBio CLR data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	200
		--diff_ratio_merging_DEL	0.5

	For PacBio CCS(HIFI) data:
		--max_cluster_bias_INS		1000
		--diff_ratio_merging_INS	0.9
		--max_cluster_bias_DEL	1000
		--diff_ratio_merging_DEL	0.5

	For ONT data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	100
		--diff_ratio_merging_DEL	0.3


	""" % (
        VERSION
    )

    # MinSizeDel = 'For current version of cuteSV, it can detect deletions larger than this size.'


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        prog="cuddlySV",
        description=cuddlySVdp.USAGE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version="%(prog)s {version}".format(version=VERSION),
    )

    # **************Parameters of input******************
    parser.add_argument(
        "input",
        metavar="[BAM]",
        type=str,
        help="Sorted .bam file from NGMLR or Minimap2.",
    )
    parser.add_argument(
        "reference", type=str, help="The reference genome in fasta format."
    )
    parser.add_argument("output", type=str, help="Output VCF format file.")
    parser.add_argument(
        "work_dir", type=str, help="Work-directory for distributed jobs"
    )

    # ************** Other Parameters******************
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use.[%(default)s]",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--verbose",
        help="Verbose output.[%(default)s]",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-b",
        "--batches",
        help="Batch of genome segmentation interval.[%(default)s]",
        default=10000000,
        type=int,
    )
    # The description of batches needs to improve.
    parser.add_argument(
        "-S", "--sample", help="Sample name/id", default="NULL", type=str
    )

    parser.add_argument(
        "--retain_work_dir",
        help="Enable to retain temporary folder and files.",
        action="store_true",
    )

    parser.add_argument(
        "--report_readid",
        help="Enable to report supporting read ids for each SV.",
        action="store_true",
    )
    parser.add_argument(
        "--report_readgroup",
        help="Enable to report supporting readgroup ids for each SV. Only usefull with report_readid.",
        action="store_true",
    )
    parser.add_argument(
        "--max_ref_allele",
        help="Maximum length of reported reference allele for deletions. Longer variants use ALT allele <DEL>.[%(default)s]",
        default=100,
        type=int,
    )
    # **************Parameters in signatures collection******************
    GroupSignaturesCollect = parser.add_argument_group("Collection of SV signatures")
    GroupSignaturesCollect.add_argument(
        "-p",
        "--max_split_parts",
        help="Maximum number of split segments a read may be aligned before it is ignored. All split segments are considered when using -1. \
			(Recommand -1 when applying assembly-based alignment.)[%(default)s]",
        default=7,
        type=int,
    )
    GroupSignaturesCollect.add_argument(
        "-q",
        "--min_mapq",
        help="Minimum mapping quality value of alignment to be taken into account.[%(default)s]",
        default=20,
        type=int,
    )
    GroupSignaturesCollect.add_argument(
        "-r",
        "--min_read_len",
        help="Ignores reads that only report alignments with not longer than bp.[%(default)s]",
        default=500,
        type=int,
    )
    GroupSignaturesCollect.add_argument(
        "-md",
        "--merge_del_threshold",
        help="Maximum distance of deletion signals to be merged. In our paper, I used -md 500 to process HG002 real human sample data.[%(default)s]",
        default=0,
        type=int,
    )
    GroupSignaturesCollect.add_argument(
        "-mi",
        "--merge_ins_threshold",
        help="Maximum distance of insertion signals to be merged. In our paper, I used -mi 500 to process HG002 real human sample data.[%(default)s]",
        default=100,
        type=int,
    )
    GroupSignaturesCollect.add_argument(
        "-include_bed",
        help="Optional given bed file. Only detect SVs in regions in the BED file. [NULL]",
        default=None,
        type=str,
    )
    # The min_read_len in last version is 2000.
    # signatures with overlap need to be filtered

    # **************Parameters in clustering******************
    GroupSVCluster = parser.add_argument_group("Generation of SV clusters")
    GroupSVCluster.add_argument(
        "-s",
        "--min_support",
        help="Minimum number of reads that support a SV to be reported.[%(default)s]",
        default=10,
        type=int,
    )
    GroupSVCluster.add_argument(
        "-l",
        "--min_size",
        help="Minimum size of SV to be reported.[%(default)s]",
        default=30,
        type=int,
    )
    GroupSVCluster.add_argument(
        "-L",
        "--max_size",
        help="Maximum size of SV to be reported. All SVs are reported when using -1. [%(default)s]",
        default=-1,
        type=int,
    )
    GroupSVCluster.add_argument(
        "-sl",
        "--min_siglength",
        help="Minimum length of SV signal to be extracted.[%(default)s]",
        default=10,
        type=int,
    )

    # **************Parameters in genotyping******************
    GroupGenotype = parser.add_argument_group("Computing genotypes")
    GroupGenotype.add_argument(
        "--genotype", help="Enable to generate genotypes.", action="store_true"
    )
    GroupGenotype.add_argument(
        "--gt_round",
        help="Maximum round of iteration for alignments searching if perform genotyping.[%(default)s]",
        default=500,
        type=int,
    )
    # GroupGenotype.add_argument('--hom',
    # 	help = "Threshold on allele frequency for homozygous.[%(default)s]",
    # 	default = 0.8,
    # 	type = float)
    # GroupGenotype.add_argument('--het',
    # 	help = "Threshold on allele frequency for heterozygous.[%(default)s].",
    # 	default = 0.2,
    # 	type = float)

    # Just a parameter for debug.
    # Will be removed in future.
    # GroupSVCluster.add_argument('--preset',
    # 	help = "Parameter presets for different sequencing technologies (pbclr/pbccs/ont).[%(default)s]",
    # 	default = "pbccs",
    # 	type = str)

    # **************Parameters in force calling******************
    GroupGenotype = parser.add_argument_group("Force calling")
    GroupGenotype.add_argument(
        "-Ivcf",  #'--MERGED_VCF',
        help="Optional given vcf file. Enable to perform force calling. [NULL]",
        default=None,
        type=str,
    )

    # **************Advanced Parameters******************
    GroupAdvanced = parser.add_argument_group("Advanced")

    # ++++++INS++++++
    GroupAdvanced.add_argument(
        "--max_cluster_bias_INS",
        help="Maximum distance to cluster read together for insertion.[%(default)s]",
        default=100,
        type=int,
    )
    GroupAdvanced.add_argument(
        "--diff_ratio_merging_INS",
        help="Do not merge breakpoints with basepair identity more than [%(default)s] for insertion.",
        default=0.3,
        type=float,
    )
    # GroupAdvanced.add_argument('--diff_ratio_filtering_INS',
    # 	help = "Filter breakpoints with basepair identity less than [%(default)s] for insertion.",
    # 	default = 0.6,
    # 	type = float)

    # ++++++DEL++++++
    GroupAdvanced.add_argument(
        "--max_cluster_bias_DEL",
        help="Maximum distance to cluster read together for deletion.[%(default)s]",
        default=200,
        type=int,
    )
    GroupAdvanced.add_argument(
        "--diff_ratio_merging_DEL",
        help="Do not merge breakpoints with basepair identity more than [%(default)s] for deletion.",
        default=0.5,
        type=float,
    )
    # GroupAdvanced.add_argument('--diff_ratio_filtering_DEL',
    # 	help = "Filter breakpoints with basepair identity less than [%(default)s] for deletion.",
    # 	default = 0.7,
    # 	type = float)

    # ++++++INV++++++
    GroupAdvanced.add_argument(
        "--max_cluster_bias_INV",
        help="Maximum distance to cluster read together for inversion.[%(default)s]",
        default=500,
        type=int,
    )

    # ++++++DUP++++++
    GroupAdvanced.add_argument(
        "--max_cluster_bias_DUP",
        help="Maximum distance to cluster read together for duplication.[%(default)s]",
        default=500,
        type=int,
    )

    # ++++++TRA++++++
    GroupAdvanced.add_argument(
        "--max_cluster_bias_TRA",
        help="Maximum distance to cluster read together for translocation.[%(default)s]",
        default=50,
        type=int,
    )
    GroupAdvanced.add_argument(
        "--diff_ratio_filtering_TRA",
        help="Filter breakpoints with basepair identity less than [%(default)s] for translocation.",
        default=0.6,
        type=float,
    )

    GroupAdvanced.add_argument(
        "--remain_reads_ratio",
        help="The ratio of reads remained in cluster. Set lower when the alignment data have high quality but recommand over 0.5.[%(default)s]",
        default=1.0,
        type=float,
    )

    # parser.add_argument('-d', '--max_distance',
    # 	help = "Maximum distance to group SV together..[%(default)s]",
    # 	default = 1000, type = int)

    # These parameters are drawn lessons from pbsv v2.2.0
    # parser.add_argument('--min_del_size',
    # 	help = "Minimum size of a deletion.[%(default)s]",
    # 	default = 20, type = int)

    args = parser.parse_args(argv)
    return args


def Generation_VCF_header(file, contiginfo, sample, argv):
    # General header
    file.write("##fileformat=VCFv4.2\n")
    file.write("##source=cuddlySV-%s\n" % (VERSION))
    import time

    file.write(
        "##fileDate=%s\n" % (time.strftime("%Y-%m-%d %H:%M:%S %w-%Z", time.localtime()))
    )
    for i in contiginfo:
        file.write("##contig=<ID=%s,length=%d>\n" % (i[0], i[1]))

    # Specific header
    # ALT
    file.write(
        '##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">\n'
    )
    file.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">\n')
    file.write(
        '##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">\n'
    )
    file.write('##ALT=<ID=INV,Description="Inversion of reference sequence">\n')
    file.write('##ALT=<ID=BND,Description="Breakend of translocation">\n')

    # INFO
    file.write(
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">\n'
    )
    file.write(
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">\n'
    )
    file.write(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    )
    file.write(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
    )
    file.write(
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n'
    )
    file.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
    )
    file.write(
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Min,Max interval around POS for imprecise variants">\n'
    )
    file.write(
        '##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Min,Max interval around inserted/deleted material between breakends">\n'
    )
    # file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    file.write(
        '##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">\n'
    )
    file.write(
        '##INFO=<ID=STRAND,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">\n'
    )
    file.write(
        '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names of SVs (comma separated)">\n'
    )
    file.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">\n')
    file.write('##FILTER=<ID=q5,Description="Quality below 5">\n')
    # FORMAT
    # file.write("\n")
    file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    file.write(
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# High-quality reference reads">\n'
    )
    file.write(
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">\n'
    )
    file.write(
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="# Phred-scaled genotype likelihoods rounded to the closest integer">\n'
    )
    file.write(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality">\n'
    )

    file.write('##CommandLine="cuddlySV %s"\n' % (" ".join(argv)))


def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s:%(process)d] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    # logging.info("Running %s" % " ".join(sys.argv))
    logging.debug("Since debug is %s, using logLevel %s", str(debug), str(logLevel))

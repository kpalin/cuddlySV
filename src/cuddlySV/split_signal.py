from collections import namedtuple
from typing import List
import cigar

SplitRead = namedtuple(
    "SplitRead", ["read_start", "read_end", "ref_start", "ref_end", "chrom", "strand"]
)


def acquire_clip_pos(deal_cigar):
    seq = list(cigar.Cigar(deal_cigar).items())
    if seq[0][1] == "S":
        first_pos = seq[0][0]
    else:
        first_pos = 0
    if seq[-1][1] == "S":
        last_pos = seq[-1][0]
    else:
        last_pos = 0

    ref_span = 0
    for i in seq:
        if i[1] == "M" or i[1] == "D" or i[1] == "=" or i[1] == "X":
            ref_span += i[0]
    return [first_pos, last_pos, ref_span]


def organize_split_signal(
    primary_info,
    Supplementary_info,
    total_L,
    SV_size,
    min_mapq,
    max_split_parts,
    read_name,
    candidate,
    MaxSize,
    query,
):
    split_read: List[SplitRead] = list()
    if len(primary_info) > 0:
        split_read.append(SplitRead(*primary_info))
        min_mapq = 0
    for i in Supplementary_info:
        seq = i.split(",")
        local_chr = seq[0]
        local_start = int(seq[1])
        local_cigar = seq[3]
        local_strand = seq[2]
        local_mapq = int(seq[4])
        if local_mapq >= min_mapq:
            # if local_mapq >= 0:
            first_pos, last_pos, ref_span = acquire_clip_pos(local_cigar)
            if local_strand == "+":
                split_read.append(
                    SplitRead(
                        first_pos,
                        total_L - last_pos,
                        local_start,
                        local_start + ref_span - 1,
                        local_chr,
                        local_strand,
                    )
                )
            else:
                split_read.append(
                    SplitRead(
                        last_pos,
                        total_L - first_pos,
                        local_start,
                        local_start + ref_span - 1,
                        local_chr,
                        local_strand,
                    )
                )
    if len(split_read) <= max_split_parts or max_split_parts == -1:
        analysis_split_read(
            split_read, SV_size, total_L, read_name, candidate, MaxSize, query
        )


def analysis_bnd(ele_1: SplitRead, ele_2: SplitRead, read_name: str, candidate: List):
    """
    *********Description*********
    *	TYPE A:		N[chr:pos[	*
    *	TYPE B:		N]chr:pos]	*
    *	TYPE C:		[chr:pos[N	*
    *	TYPE D:		]chr:pos]N	*
    *****************************
    """
    if ele_2.read_start - ele_1.read_end <= 100:
        if ele_1.strand == "+":
            if ele_2.strand == "+":
                # +&+
                if ele_1.chrom < ele_2.chrom:
                    candidate.append(
                        [
                            "A",
                            ele_1.ref_end,
                            ele_2.chrom,
                            ele_2.ref_start,
                            read_name,
                            "TRA",
                            ele_1.chrom,
                        ]
                    )
                    # N[chr:pos[
                else:
                    candidate.append(
                        [
                            "D",
                            ele_2.ref_start,
                            ele_1.chrom,
                            ele_1.ref_end,
                            read_name,
                            "TRA",
                            ele_2.chrom,
                        ]
                    )
                    # ]chr:pos]N
            else:
                # +&-
                if ele_1.chrom < ele_2.chrom:
                    candidate.append(
                        [
                            "B",
                            ele_1.ref_end,
                            ele_2.chrom,
                            ele_2.ref_end,
                            read_name,
                            "TRA",
                            ele_1.chrom,
                        ]
                    )
                    # N]chr:pos]
                else:
                    candidate.append(
                        [
                            "B",
                            ele_2.ref_end,
                            ele_1.chrom,
                            ele_1.ref_end,
                            read_name,
                            "TRA",
                            ele_2.chrom,
                        ]
                    )
                    # N]chr:pos]
        else:
            if ele_2.strand == "+":
                # -&+
                if ele_1.chrom < ele_2.chrom:
                    candidate.append(
                        [
                            "C",
                            ele_1.ref_start,
                            ele_2.chrom,
                            ele_2.ref_start,
                            read_name,
                            "TRA",
                            ele_1.chrom,
                        ]
                    )
                    # [chr:pos[N
                else:
                    candidate.append(
                        [
                            "C",
                            ele_2.ref_start,
                            ele_1.chrom,
                            ele_1.ref_start,
                            read_name,
                            "TRA",
                            ele_2.chrom,
                        ]
                    )
                    # [chr:pos[N
            else:
                # -&-
                if ele_1.chrom < ele_2.chrom:
                    candidate.append(
                        [
                            "D",
                            ele_1.ref_start,
                            ele_2.chrom,
                            ele_2.ref_end,
                            read_name,
                            "TRA",
                            ele_1.chrom,
                        ]
                    )
                    # ]chr:pos]N
                else:
                    candidate.append(
                        [
                            "A",
                            ele_2.ref_end,
                            ele_1.chrom,
                            ele_1.ref_start,
                            read_name,
                            "TRA",
                            ele_2.chrom,
                        ]
                    )
                    # N[chr:pos[


def is_1d2_pair(ele_1: SplitRead, ele_2: SplitRead) -> bool:
    """Check if the split reads form merged 1d2 molecule

    Args:
        ele_1 (SplitRead): One split
        ele_2 (SplitRead): Next split

    Returns:
        bool: True if the input splits are same chromosome and opposite strand
    """
    if ele_1.chrom == ele_2.chrom and ele_1.strand != ele_2.strand:
        overlap_len = min(ele_1.ref_end, ele_2.ref_end) - max(
            ele_1.ref_start, ele_2.ref_start
        )
        shorter_ele_len = min(
            ele_1.ref_end - ele_1.ref_start, ele_2.ref_end - ele_2.ref_start
        )
        return overlap_len / shorter_ele_len > 0.95
    else:
        return False


def analysis_inv(
    ele_1: SplitRead, ele_2: SplitRead, read_name: str, candidate: List, SV_size: int
):
    if is_1d2_pair(ele_1, ele_2):
        return
    if ele_1.strand == "+":
        # +-
        if ele_1.ref_end - ele_2.ref_end >= SV_size:
            if (
                ele_2.read_start + 0.5 * (ele_1.ref_end - ele_2.ref_end)
                >= ele_1.read_end
            ):
                candidate.append(
                    ["++", ele_2.ref_end, ele_1.ref_end, read_name, "INV", ele_1.chrom]
                )
                # head-to-head
                # 5'->5'
        if ele_2.ref_end - ele_1.ref_end >= SV_size:
            if (
                ele_2.read_start + 0.5 * (ele_2.ref_end - ele_1.ref_end)
                >= ele_1.read_end
            ):
                candidate.append(
                    ["++", ele_1.ref_end, ele_2.ref_end, read_name, "INV", ele_1.chrom]
                )
                # head-to-head
                # 5'->5'
    else:
        # -+
        if ele_2.ref_start - ele_1.ref_start >= SV_size:
            if (
                ele_2.read_start + 0.5 * (ele_2.ref_start - ele_1.ref_start)
                >= ele_1.read_end
            ):
                candidate.append(
                    [
                        "--",
                        ele_1.ref_start,
                        ele_2.ref_start,
                        read_name,
                        "INV",
                        ele_1.chrom,
                    ]
                )
                # tail-to-tail
                # 3'->3'
        if ele_1.ref_start - ele_2.ref_start >= SV_size:
            if (
                ele_2.read_start + 0.5 * (ele_1.ref_start - ele_2.ref_start)
                >= ele_1.read_end
            ):
                candidate.append(
                    [
                        "--",
                        ele_2.ref_start,
                        ele_1.ref_start,
                        read_name,
                        "INV",
                        ele_1.chrom,
                    ]
                )
                # tail-to-tail
                # 3'->3'


def analysis_split_read(
    split_read: List[SplitRead],
    SV_size: int,
    RLength: int,
    read_name: str,
    candidate: List,
    MaxSize: int,
    query: str,
) -> None:
    """
    read_start	read_end	ref_start	ref_end	chr	strand
    #0			#1			#2			#3		#4	#5
    """
    SP_list = sorted(split_read, key=lambda x: x.read_start)
    # print(read_name)
    # for i in SP_list:
    # 	print(i)

    # Store Strands of INV

    for sp_idx in range(len(SP_list) - 1):
        ele_1, ele_2 = SP_list[sp_idx : (sp_idx + 2)]

        if ele_1.chrom == ele_2.chrom:
            analyse_ins_dup_del(
                SV_size, RLength, read_name, candidate, MaxSize, query, ele_1, ele_2
            )
        else:
            analysis_bnd(ele_1, ele_2, read_name, candidate)

    if len(SP_list) >= 3:
        # detect INS involoved in a translocation
        trigger_INS_TRA = 0
        # over three splits
        for a in range(len(SP_list[1:-1])):
            query, trigger_INS_TRA = analyse_over_three_splits(
                SV_size, RLength, read_name, candidate, MaxSize, query, SP_list, a
            )

        if trigger_INS_TRA == 1:
            analyse_insertion_in_translocation(
                SV_size, RLength, read_name, candidate, MaxSize, query, SP_list
            )


def analyse_ins_dup_del(
    SV_size,
    RLength,
    read_name,
    candidate: List,
    MaxSize,
    query,
    ele_1: SplitRead,
    ele_2: SplitRead,
):
    if ele_1.strand != ele_2.strand:
        # Same chromosome but opposite strand
        analysis_inv(ele_1, ele_2, read_name, candidate, SV_size)

    else:
        # Same chromosome and strand
        # dup & ins & del
        if ele_1.strand == "-":
            ele_1, ele_2 = (
                SplitRead(
                    RLength - ele_2.read_end, RLength - ele_2.read_start, *ele_2[2:]
                ),
                SplitRead(
                    RLength - ele_1.read_end, RLength - ele_1.read_start, *ele_1[2:]
                ),
            )
            query = query[::-1]

        if ele_1.ref_end - ele_2.ref_start >= SV_size:
            # if ele_2.read_end - ele_1.read_end >= ele_1.ref_end - ele_2.ref_start:
            if ele_2.read_start - ele_1.read_end >= ele_1.ref_end - ele_2.ref_start:
                candidate.append(
                    [
                        (ele_1.ref_end + ele_2.ref_start) / 2,
                        ele_2.read_start
                        + ele_1.ref_end
                        - ele_2.ref_start
                        - ele_1.read_end,
                        read_name,
                        str(
                            query[
                                ele_1.read_end
                                + int(
                                    (ele_1.ref_end - ele_2.ref_start) / 2
                                ) : ele_2.read_start
                                - int((ele_1.ref_end - ele_2.ref_start) / 2)
                            ]
                        ),
                        "INS",
                        ele_2.chrom,
                    ]
                )
            else:
                candidate.append(
                    [
                        ele_2.ref_start,
                        ele_1.ref_end,
                        read_name,
                        "DUP",
                        ele_2.chrom,
                    ]
                )

        delta_length = (
            ele_2.read_start + ele_1.ref_end - ele_2.ref_start - ele_1.read_end
        )
        if (
            ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
            and delta_length >= SV_size
        ):
            if ele_2.ref_start - ele_1.ref_end <= max(100, delta_length / 5) and (
                delta_length <= MaxSize or MaxSize == -1
            ):
                candidate.append(
                    [
                        (ele_2.ref_start + ele_1.ref_end) / 2,
                        delta_length,
                        read_name,
                        str(
                            query[
                                ele_1.read_end
                                + int(
                                    (ele_2.ref_start - ele_1.ref_end) / 2
                                ) : ele_2.read_start
                                - int((ele_2.ref_start - ele_1.ref_end) / 2)
                            ]
                        ),
                        "INS",
                        ele_2.chrom,
                    ]
                )
        delta_length = (
            ele_2.ref_start - ele_2.read_start + ele_1.read_end - ele_1.ref_end
        )
        if (
            ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
            and delta_length >= SV_size
        ):
            if ele_2.read_start - ele_1.read_end <= max(100, delta_length / 5) and (
                delta_length <= MaxSize or MaxSize == -1
            ):
                candidate.append(
                    [ele_1.ref_end, delta_length, read_name, "DEL", ele_2.chrom]
                )

    return ele_1, ele_2


def analyse_insertion_in_translocation(
    SV_size, RLength, read_name, candidate, MaxSize, query, SP_list
):
    if SP_list[0].chrom == SP_list[-1].chrom:
        # print(SP_list[0])
        # print(SP_list[-1])
        if SP_list[0].strand != SP_list[-1].strand:
            pass
        else:
            if SP_list[0].strand == "+":
                ele_1 = SP_list[0]
                ele_2 = SP_list[-1]
            else:
                ele_1 = SplitRead(
                    RLength - SP_list[-1].read_end,
                    RLength - SP_list[-1].read_start,
                    *SP_list[-1][2:],
                )
                ele_2 = SplitRead(
                    RLength - SP_list[0].read_end,
                    RLength - SP_list[0].read_start,
                    *SP_list[0][2:],
                )
                query = query[::-1]
                # print(ele_1)
                # print(ele_2)
            dis_ref = ele_2.ref_start - ele_1.ref_end
            dis_read = ele_2.read_start - ele_1.read_end
            if (
                dis_ref < 100
                and dis_read - dis_ref >= SV_size
                and (dis_read - dis_ref <= MaxSize or MaxSize == -1)
            ):
                # print(min(ele_2.ref_start, ele_1.ref_end), dis_read - dis_ref, read_name)
                candidate.append(
                    [
                        min(ele_2.ref_start, ele_1.ref_end),
                        dis_read - dis_ref,
                        read_name,
                        str(
                            query[
                                ele_1.read_end + int(dis_ref / 2) : ele_2.read_start
                                - int(dis_ref / 2)
                            ]
                        ),
                        "INS",
                        ele_2.chrom,
                    ]
                )

            if dis_ref <= -SV_size:
                candidate.append(
                    [ele_2.ref_start, ele_1.ref_end, read_name, "DUP", ele_2.chrom]
                )


def analyse_over_three_splits(
    SV_size,
    RLength,
    read_name,
    candidate: List,
    MaxSize,
    query,
    SP_list: List[SplitRead],
    sp_idx: int,
):
    ele_1 = SP_list[sp_idx]
    ele_2 = SP_list[sp_idx + 1]
    ele_3 = SP_list[sp_idx + 2]
    trigger_INS_TRA = 0

    if ele_1.chrom == ele_2.chrom:
        if ele_2.chrom == ele_3.chrom:
            if ele_1.strand == ele_3.strand and ele_1.strand != ele_2.strand:
                if ele_2.strand == "-":
                    # +-+
                    if (
                        ele_2.read_start + 0.5 * (ele_3.ref_start - ele_1.ref_end)
                        >= ele_1.read_end
                        and ele_3.read_start + 0.5 * (ele_3.ref_start - ele_1.ref_end)
                        >= ele_2.read_end
                    ):
                        # No overlaps in split reads
                        if (
                            ele_2.ref_start >= ele_1.ref_end
                            and ele_3.ref_start >= ele_2.ref_end
                        ):
                            candidate.append(
                                [
                                    "++",
                                    ele_1.ref_end,
                                    ele_2.ref_end,
                                    read_name,
                                    "INV",
                                    ele_1.chrom,
                                ]
                            )
                            # head-to-head
                            # 5'->5'
                            candidate.append(
                                [
                                    "--",
                                    ele_2.ref_start,
                                    ele_3.ref_start,
                                    read_name,
                                    "INV",
                                    ele_1.chrom,
                                ]
                            )
                            # tail-to-tail
                            # 3'->3'
                else:
                    # -+-
                    if (
                        ele_1.read_end
                        <= ele_2.read_start + 0.5 * (ele_1.ref_start - ele_3.ref_end)
                        and ele_3.read_start + 0.5 * (ele_1.ref_start - ele_3.ref_end)
                        >= ele_2.read_end
                    ):
                        # No overlaps in split reads
                        if (
                            ele_2.ref_start - ele_3.ref_end >= -50
                            and ele_1.ref_start - ele_2.ref_end >= -50
                        ):
                            candidate.append(
                                [
                                    "++",
                                    ele_3.ref_end,
                                    ele_2.ref_end,
                                    read_name,
                                    "INV",
                                    ele_1.chrom,
                                ]
                            )
                            # head-to-head
                            # 5'->5'
                            candidate.append(
                                [
                                    "--",
                                    ele_2.ref_start,
                                    ele_1.ref_start,
                                    read_name,
                                    "INV",
                                    ele_1.chrom,
                                ]
                            )
                            # tail-to-tail
                            # 3'->3'

            if len(SP_list) - 3 == sp_idx:
                if ele_1.strand != ele_3.strand:
                    if ele_2.strand == ele_1.strand:
                        # ++-/--+
                        analysis_inv(ele_2, ele_3, read_name, candidate, SV_size)
                    else:
                        # +--/-++
                        analysis_inv(ele_1, ele_2, read_name, candidate, SV_size)

            if ele_1.strand == ele_3.strand and ele_1.strand == ele_2.strand:
                # dup & ins & del
                if ele_1.strand == "-":
                    ele_1 = SplitRead(
                        RLength - SP_list[sp_idx + 2].read_end,
                        RLength - SP_list[sp_idx + 2].read_start,
                        *SP_list[sp_idx + 2][2:],
                    )
                    ele_2 = SplitRead(
                        RLength - SP_list[sp_idx + 1].read_end,
                        RLength - SP_list[sp_idx + 1].read_start,
                        *SP_list[sp_idx + 1][2:],
                    )
                    ele_3 = SplitRead(
                        RLength - SP_list[sp_idx].read_end,
                        RLength - SP_list[sp_idx].read_start,
                        *SP_list[sp_idx][2:],
                    )
                    query = query[::-1]

                if (
                    ele_2.ref_end - ele_3.ref_start >= SV_size
                    and ele_2.ref_start < ele_3.ref_end
                ):
                    candidate.append(
                        [
                            ele_3.ref_start,
                            ele_2.ref_end,
                            read_name,
                            "DUP",
                            ele_2.chrom,
                        ]
                    )

                if sp_idx == 0:
                    if ele_1.ref_end - ele_2.ref_start >= SV_size:
                        candidate.append(
                            [
                                ele_2.ref_start,
                                ele_1.ref_end,
                                read_name,
                                "DUP",
                                ele_2.chrom,
                            ]
                        )

                delta_length = (
                    ele_2.read_start + ele_1.ref_end - ele_2.ref_start - ele_1.read_end
                )
                if (
                    ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2.ref_start - ele_1.ref_end <= max(
                        100, delta_length / 5
                    ) and (delta_length <= MaxSize or MaxSize == -1):
                        if ele_3.ref_start >= ele_2.ref_end:
                            candidate.append(
                                [
                                    (ele_2.ref_start + ele_1.ref_end) / 2,
                                    delta_length,
                                    read_name,
                                    str(
                                        query[
                                            ele_1.read_end
                                            + int(
                                                (ele_2.ref_start - ele_1.ref_end) / 2
                                            ) : ele_2.read_start
                                            - int((ele_2.ref_start - ele_1.ref_end) / 2)
                                        ]
                                    ),
                                    "INS",
                                    ele_2.chrom,
                                ]
                            )
                delta_length = (
                    ele_2.ref_start - ele_2.read_start + ele_1.read_end - ele_1.ref_end
                )
                if (
                    ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2.read_start - ele_1.read_end <= max(
                        100, delta_length / 5
                    ) and (delta_length <= MaxSize or MaxSize == -1):
                        if ele_3.ref_start >= ele_2.ref_end:
                            candidate.append(
                                [
                                    ele_1.ref_end,
                                    delta_length,
                                    read_name,
                                    "DEL",
                                    ele_2.chrom,
                                ]
                            )

                if len(SP_list) - 3 == sp_idx:
                    ele_1 = ele_2
                    ele_2 = ele_3

                    delta_length = (
                        ele_2.read_start
                        + ele_1.ref_end
                        - ele_2.ref_start
                        - ele_1.read_end
                    )
                    if (
                        ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                        and delta_length >= SV_size
                    ):
                        if ele_2.ref_start - ele_1.ref_end <= max(
                            100, delta_length / 5
                        ) and (delta_length <= MaxSize or MaxSize == -1):
                            candidate.append(
                                [
                                    (ele_2.ref_start + ele_1.ref_end) / 2,
                                    delta_length,
                                    read_name,
                                    str(
                                        query[
                                            ele_1.read_end
                                            + int(
                                                (ele_2.ref_start - ele_1.ref_end) / 2
                                            ) : ele_2.read_start
                                            - int((ele_2.ref_start - ele_1.ref_end) / 2)
                                        ]
                                    ),
                                    "INS",
                                    ele_2.chrom,
                                ]
                            )

                    delta_length = (
                        ele_2.ref_start
                        - ele_2.read_start
                        + ele_1.read_end
                        - ele_1.ref_end
                    )
                    if (
                        ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                        and ele_2.ref_start
                        - ele_2.read_start
                        + ele_1.read_end
                        - ele_1.ref_end
                        >= SV_size
                    ):
                        if ele_2.read_start - ele_1.read_end <= max(
                            100, delta_length / 5
                        ) and (delta_length <= MaxSize or MaxSize == -1):
                            candidate.append(
                                [
                                    ele_1.ref_end,
                                    delta_length,
                                    read_name,
                                    "DEL",
                                    ele_2.chrom,
                                ]
                            )

            if (
                len(SP_list) - 3 == sp_idx
                and ele_1.strand != ele_2.strand
                and ele_2.strand == ele_3.strand
            ):
                ele_1 = ele_2
                ele_2 = ele_3
                ele_3 = None
            if ele_3 is None or (
                ele_1.strand == ele_2.strand and ele_2.strand != ele_3.strand
            ):
                if ele_1.strand == "-":
                    ele_1 = SplitRead(
                        RLength - SP_list[sp_idx + 2].read_end,
                        RLength - SP_list[sp_idx + 2].read_start,
                        *SP_list[sp_idx + 2][2:],
                    )
                    ele_2 = SplitRead(
                        RLength - SP_list[sp_idx + 1].read_end,
                        RLength - SP_list[sp_idx + 1].read_start,
                        *SP_list[sp_idx + 1][2:],
                    )
                    query = query[::-1]
                delta_length = (
                    ele_2.read_start + ele_1.ref_end - ele_2.ref_start - ele_1.read_end
                )
                if (
                    ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2.ref_start - ele_1.ref_end <= max(
                        100, delta_length / 5
                    ) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append(
                            [
                                (ele_2.ref_start + ele_1.ref_end) / 2,
                                delta_length,
                                read_name,
                                str(
                                    query[
                                        ele_1.read_end
                                        + int(
                                            (ele_2.ref_start - ele_1.ref_end) / 2
                                        ) : ele_2.read_start
                                        - int((ele_2.ref_start - ele_1.ref_end) / 2)
                                    ]
                                ),
                                "INS",
                                ele_2.chrom,
                            ]
                        )

                delta_length = (
                    ele_2.ref_start - ele_2.read_start + ele_1.read_end - ele_1.ref_end
                )
                if (
                    ele_1.ref_end - ele_2.ref_start < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2.read_start - ele_1.read_end <= max(
                        100, delta_length / 5
                    ) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append(
                            [
                                ele_1.ref_end,
                                delta_length,
                                read_name,
                                "DEL",
                                ele_2.chrom,
                            ]
                        )

    else:
        trigger_INS_TRA = 1
        analysis_bnd(ele_1, ele_2, read_name, candidate)

        if len(SP_list) - 3 == sp_idx:
            if ele_2.chrom != ele_3.chrom:
                analysis_bnd(ele_2, ele_3, read_name, candidate)
    return query, trigger_INS_TRA

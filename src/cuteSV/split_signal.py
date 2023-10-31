from collections import namedtuple
from typing import List
from cuteSV.cuddlySV import acquire_clip_pos

SplitRead = namedtuple(
    "SplitRead", ["read_start", "read_end", "ref_start", "ref_end", "chrom", "strand"]
)


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
            local_set = acquire_clip_pos(local_cigar)
            if local_strand == "+":
                split_read.append(
                    SplitRead(
                        local_set[0],
                        total_L - local_set[1],
                        local_start,
                        local_start + local_set[2],
                        local_chr,
                        local_strand,
                    )
                )
            else:
                try:
                    split_read.append(
                        SplitRead(
                            local_set[1],
                            total_L - local_set[0],
                            local_start,
                            local_start + local_set[2],
                            local_chr,
                            local_strand,
                        )
                    )
                except:
                    pass
    if len(split_read) <= max_split_parts or max_split_parts == -1:
        analysis_split_read(
            split_read, SV_size, total_L, read_name, candidate, MaxSize, query
        )


def analysis_bnd(ele_1, ele_2, read_name, candidate):
    """
    *********Description*********
    *	TYPE A:		N[chr:pos[	*
    *	TYPE B:		N]chr:pos]	*
    *	TYPE C:		[chr:pos[N	*
    *	TYPE D:		]chr:pos]N	*
    *****************************
    """
    if ele_2[0] - ele_1[1] <= 100:
        if ele_1[5] == "+":
            if ele_2[5] == "+":
                # +&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(
                        ["A", ele_1[3], ele_2[4], ele_2[2], read_name, "TRA", ele_1[4]]
                    )
                    # N[chr:pos[
                else:
                    candidate.append(
                        ["D", ele_2[2], ele_1[4], ele_1[3], read_name, "TRA", ele_2[4]]
                    )
                    # ]chr:pos]N
            else:
                # +&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(
                        ["B", ele_1[3], ele_2[4], ele_2[3], read_name, "TRA", ele_1[4]]
                    )
                    # N]chr:pos]
                else:
                    candidate.append(
                        ["B", ele_2[3], ele_1[4], ele_1[3], read_name, "TRA", ele_2[4]]
                    )
                    # N]chr:pos]
        else:
            if ele_2[5] == "+":
                # -&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(
                        ["C", ele_1[2], ele_2[4], ele_2[2], read_name, "TRA", ele_1[4]]
                    )
                    # [chr:pos[N
                else:
                    candidate.append(
                        ["C", ele_2[2], ele_1[4], ele_1[2], read_name, "TRA", ele_2[4]]
                    )
                    # [chr:pos[N
            else:
                # -&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(
                        ["D", ele_1[2], ele_2[4], ele_2[3], read_name, "TRA", ele_1[4]]
                    )
                    # ]chr:pos]N
                else:
                    candidate.append(
                        ["A", ele_2[3], ele_1[4], ele_1[2], read_name, "TRA", ele_2[4]]
                    )
                    # N[chr:pos[


def analysis_inv(ele_1, ele_2, read_name, candidate, SV_size):
    if ele_1[5] == "+":
        # +-
        if ele_1[3] - ele_2[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[3] - ele_2[3]) >= ele_1[1]:
                candidate.append(["++", ele_2[3], ele_1[3], read_name, "INV", ele_1[4]])
                # head-to-head
                # 5'->5'
        if ele_2[3] - ele_1[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[3] - ele_1[3]) >= ele_1[1]:
                candidate.append(["++", ele_1[3], ele_2[3], read_name, "INV", ele_1[4]])
                # head-to-head
                # 5'->5'
    else:
        # -+
        if ele_2[2] - ele_1[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[2] - ele_1[2]) >= ele_1[1]:
                candidate.append(["--", ele_1[2], ele_2[2], read_name, "INV", ele_1[4]])
                # tail-to-tail
                # 3'->3'
        if ele_1[2] - ele_2[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[2] - ele_2[2]) >= ele_1[1]:
                candidate.append(["--", ele_2[2], ele_1[2], read_name, "INV", ele_1[4]])
                # tail-to-tail
                # 3'->3'


# TODO:  This function has a bug. Read aligned in 3 parts, large deletion followed by false inversion is not classified as deletion signal.
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

    # detect INS involoved in a translocation
    trigger_INS_TRA = 0

    # Store Strands of INV

    if len(SP_list) == 2:
        ele_1 = SP_list[0]
        ele_2 = SP_list[1]
        if ele_1[4] == ele_2[4]:
            if ele_1[5] != ele_2[5]:
                analysis_inv(ele_1, ele_2, read_name, candidate, SV_size)

            else:
                # dup & ins & del
                a = 0
                if ele_1[5] == "-":
                    ele_1 = SplitRead(
                        RLength - SP_list[a + 1][1],
                        RLength - SP_list[a + 1][0],
                        *SP_list[a + 1][2:]
                    )
                    ele_2 = SplitRead(
                        RLength - SP_list[a][1],
                        RLength - SP_list[a][0],
                        *SP_list[a][2:]
                    )
                    query = query[::-1]

                if ele_1[3] - ele_2[2] >= SV_size:
                    # if ele_2[1] - ele_1[1] >= ele_1[3] - ele_2[2]:
                    if ele_2[0] - ele_1[1] >= ele_1[3] - ele_2[2]:
                        candidate.append(
                            [
                                (ele_1[3] + ele_2[2]) / 2,
                                ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1],
                                read_name,
                                str(
                                    query[
                                        ele_1[1]
                                        + int((ele_1[3] - ele_2[2]) / 2) : ele_2[0]
                                        - int((ele_1[3] - ele_2[2]) / 2)
                                    ]
                                ),
                                "INS",
                                ele_2[4],
                            ]
                        )
                    else:
                        candidate.append(
                            [ele_2[2], ele_1[3], read_name, "DUP", ele_2[4]]
                        )

                delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                if (
                    ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2[2] - ele_1[3] <= max(100, delta_length / 5) and (
                        delta_length <= MaxSize or MaxSize == -1
                    ):
                        candidate.append(
                            [
                                (ele_2[2] + ele_1[3]) / 2,
                                delta_length,
                                read_name,
                                str(
                                    query[
                                        ele_1[1]
                                        + int((ele_2[2] - ele_1[3]) / 2) : ele_2[0]
                                        - int((ele_2[2] - ele_1[3]) / 2)
                                    ]
                                ),
                                "INS",
                                ele_2[4],
                            ]
                        )
                delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                if (
                    ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                    and delta_length >= SV_size
                ):
                    if ele_2[0] - ele_1[1] <= max(100, delta_length / 5) and (
                        delta_length <= MaxSize or MaxSize == -1
                    ):
                        candidate.append(
                            [ele_1[3], delta_length, read_name, "DEL", ele_2[4]]
                        )
        else:
            trigger_INS_TRA = 1
            analysis_bnd(ele_1, ele_2, read_name, candidate)

    else:
        # over three splits
        for a in range(len(SP_list[1:-1])):
            ele_1 = SP_list[a]
            ele_2 = SP_list[a + 1]
            ele_3 = SP_list[a + 2]

            if ele_1[4] == ele_2[4]:
                if ele_2[4] == ele_3[4]:
                    if ele_1[5] == ele_3[5] and ele_1[5] != ele_2[5]:
                        if ele_2[5] == "-":
                            # +-+
                            if (
                                ele_2[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_1[1]
                                and ele_3[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_2[1]
                            ):
                                # No overlaps in split reads

                                if ele_2[2] >= ele_1[3] and ele_3[2] >= ele_2[3]:
                                    candidate.append(
                                        [
                                            "++",
                                            ele_1[3],
                                            ele_2[3],
                                            read_name,
                                            "INV",
                                            ele_1[4],
                                        ]
                                    )
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(
                                        [
                                            "--",
                                            ele_2[2],
                                            ele_3[2],
                                            read_name,
                                            "INV",
                                            ele_1[4],
                                        ]
                                    )
                                    # tail-to-tail
                                    # 3'->3'
                        else:
                            # -+-
                            if (
                                ele_1[1] <= ele_2[0] + 0.5 * (ele_1[2] - ele_3[3])
                                and ele_3[0] + 0.5 * (ele_1[2] - ele_3[3]) >= ele_2[1]
                            ):
                                # No overlaps in split reads

                                if (
                                    ele_2[2] - ele_3[3] >= -50
                                    and ele_1[2] - ele_2[3] >= -50
                                ):
                                    candidate.append(
                                        [
                                            "++",
                                            ele_3[3],
                                            ele_2[3],
                                            read_name,
                                            "INV",
                                            ele_1[4],
                                        ]
                                    )
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(
                                        [
                                            "--",
                                            ele_2[2],
                                            ele_1[2],
                                            read_name,
                                            "INV",
                                            ele_1[4],
                                        ]
                                    )
                                    # tail-to-tail
                                    # 3'->3'

                    if len(SP_list) - 3 == a:
                        if ele_1[5] != ele_3[5]:
                            if ele_2[5] == ele_1[5]:
                                # ++-/--+
                                analysis_inv(
                                    ele_2, ele_3, read_name, candidate, SV_size
                                )
                            else:
                                # +--/-++
                                analysis_inv(
                                    ele_1, ele_2, read_name, candidate, SV_size
                                )

                    if ele_1[5] == ele_3[5] and ele_1[5] == ele_2[5]:
                        # dup & ins & del
                        if ele_1[5] == "-":
                            ele_1 = SplitRead(
                                RLength - SP_list[a + 2][1],
                                RLength - SP_list[a + 2][0],
                                *SP_list[a + 2][2:]
                            )
                            ele_2 = SplitRead(
                                RLength - SP_list[a + 1][1],
                                RLength - SP_list[a + 1][0],
                                *SP_list[a + 1][2:]
                            )
                            ele_3 = SplitRead(
                                RLength - SP_list[a][1],
                                RLength - SP_list[a][0],
                                *SP_list[a][2:]
                            )
                            query = query[::-1]

                        if ele_2[3] - ele_3[2] >= SV_size and ele_2[2] < ele_3[3]:
                            candidate.append(
                                [ele_3[2], ele_2[3], read_name, "DUP", ele_2[4]]
                            )

                        if a == 0:
                            if ele_1[3] - ele_2[2] >= SV_size:
                                candidate.append(
                                    [ele_2[2], ele_1[3], read_name, "DUP", ele_2[4]]
                                )

                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if (
                            ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                            and delta_length >= SV_size
                        ):
                            if ele_2[2] - ele_1[3] <= max(100, delta_length / 5) and (
                                delta_length <= MaxSize or MaxSize == -1
                            ):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append(
                                        [
                                            (ele_2[2] + ele_1[3]) / 2,
                                            delta_length,
                                            read_name,
                                            str(
                                                query[
                                                    ele_1[1]
                                                    + int(
                                                        (ele_2[2] - ele_1[3]) / 2
                                                    ) : ele_2[0]
                                                    - int((ele_2[2] - ele_1[3]) / 2)
                                                ]
                                            ),
                                            "INS",
                                            ele_2[4],
                                        ]
                                    )
                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if (
                            ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                            and delta_length >= SV_size
                        ):
                            if ele_2[0] - ele_1[1] <= max(100, delta_length / 5) and (
                                delta_length <= MaxSize or MaxSize == -1
                            ):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append(
                                        [
                                            ele_1[3],
                                            delta_length,
                                            read_name,
                                            "DEL",
                                            ele_2[4],
                                        ]
                                    )

                        if len(SP_list) - 3 == a:
                            ele_1 = ele_2
                            ele_2 = ele_3

                            delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                            if (
                                ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                                and delta_length >= SV_size
                            ):
                                if ele_2[2] - ele_1[3] <= max(
                                    100, delta_length / 5
                                ) and (delta_length <= MaxSize or MaxSize == -1):
                                    candidate.append(
                                        [
                                            (ele_2[2] + ele_1[3]) / 2,
                                            delta_length,
                                            read_name,
                                            str(
                                                query[
                                                    ele_1[1]
                                                    + int(
                                                        (ele_2[2] - ele_1[3]) / 2
                                                    ) : ele_2[0]
                                                    - int((ele_2[2] - ele_1[3]) / 2)
                                                ]
                                            ),
                                            "INS",
                                            ele_2[4],
                                        ]
                                    )

                            delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                            if (
                                ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                                and ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size
                            ):
                                if ele_2[0] - ele_1[1] <= max(
                                    100, delta_length / 5
                                ) and (delta_length <= MaxSize or MaxSize == -1):
                                    candidate.append(
                                        [
                                            ele_1[3],
                                            delta_length,
                                            read_name,
                                            "DEL",
                                            ele_2[4],
                                        ]
                                    )

                    if (
                        len(SP_list) - 3 == a
                        and ele_1[5] != ele_2[5]
                        and ele_2[5] == ele_3[5]
                    ):
                        ele_1 = ele_2
                        ele_2 = ele_3
                        ele_3 = None
                    if ele_3 == None or (ele_1[5] == ele_2[5] and ele_2[5] != ele_3[5]):
                        if ele_1[5] == "-":
                            ele_1 = SplitRead(
                                RLength - SP_list[a + 2][1],
                                RLength - SP_list[a + 2][0],
                                *SP_list[a + 2][2:]
                            )
                            ele_2 = SplitRead(
                                RLength - SP_list[a + 1][1],
                                RLength - SP_list[a + 1][0],
                                *SP_list[a + 1][2:]
                            )
                            query = query[::-1]
                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if (
                            ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                            and delta_length >= SV_size
                        ):
                            if ele_2[2] - ele_1[3] <= max(100, delta_length / 5) and (
                                delta_length <= MaxSize or MaxSize == -1
                            ):
                                candidate.append(
                                    [
                                        (ele_2[2] + ele_1[3]) / 2,
                                        delta_length,
                                        read_name,
                                        str(
                                            query[
                                                ele_1[1]
                                                + int(
                                                    (ele_2[2] - ele_1[3]) / 2
                                                ) : ele_2[0]
                                                - int((ele_2[2] - ele_1[3]) / 2)
                                            ]
                                        ),
                                        "INS",
                                        ele_2[4],
                                    ]
                                )

                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if (
                            ele_1[3] - ele_2[2] < max(SV_size, delta_length / 5)
                            and delta_length >= SV_size
                        ):
                            if ele_2[0] - ele_1[1] <= max(100, delta_length / 5) and (
                                delta_length <= MaxSize or MaxSize == -1
                            ):
                                candidate.append(
                                    [ele_1[3], delta_length, read_name, "DEL", ele_2[4]]
                                )

            else:
                trigger_INS_TRA = 1
                analysis_bnd(ele_1, ele_2, read_name, candidate)

                if len(SP_list) - 3 == a:
                    if ele_2[4] != ele_3[4]:
                        analysis_bnd(ele_2, ele_3, read_name, candidate)

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        if SP_list[0][4] == SP_list[-1][4]:
            # print(SP_list[0])
            # print(SP_list[-1])
            if SP_list[0][5] != SP_list[-1][5]:
                pass
            else:
                if SP_list[0][5] == "+":
                    ele_1 = SP_list[0]
                    ele_2 = SP_list[-1]
                else:
                    ele_1 = SplitRead(
                        RLength - SP_list[-1][1],
                        RLength - SP_list[-1][0],
                        *SP_list[-1][2:]
                    )
                    ele_2 = SplitRead(
                        RLength - SP_list[0][1],
                        RLength - SP_list[0][0],
                        *SP_list[0][2:]
                    )
                    query = query[::-1]
                # print(ele_1)
                # print(ele_2)
                dis_ref = ele_2[2] - ele_1[3]
                dis_read = ele_2[0] - ele_1[1]
                if (
                    dis_ref < 100
                    and dis_read - dis_ref >= SV_size
                    and (dis_read - dis_ref <= MaxSize or MaxSize == -1)
                ):
                    # print(min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name)
                    candidate.append(
                        [
                            min(ele_2[2], ele_1[3]),
                            dis_read - dis_ref,
                            read_name,
                            str(
                                query[
                                    ele_1[1]
                                    + int(dis_ref / 2) : ele_2[0]
                                    - int(dis_ref / 2)
                                ]
                            ),
                            "INS",
                            ele_2[4],
                        ]
                    )

                if dis_ref <= -SV_size:
                    candidate.append([ele_2[2], ele_1[3], read_name, "DUP", ele_2[4]])

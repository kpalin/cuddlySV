# This test code was written by the `hypothesis.extra.ghostwriter` module
# and is provided under the Creative Commons Zero public domain dedication.

from pytest import mark
import cuteSV.cuddlySV
import cuteSV.split_signal
from hypothesis import given, strategies as st
from hypothesis import settings, Verbosity

# TODO: replace st.nothing() with appropriate strategies

chroms = {f"chr{i+1:d}" for i in range(22)}
chroms.update({"chrX", "chrY", "chrM"})
chroms = list(sorted(chroms))


@st.composite
def supplementary_info(draw: st.DrawFn):
    "Generate supplementary alignment info like 'chr1,30178190,-,108S6025M120D6314S,60,780'"
    chrom = draw(st.sampled_from(chroms))
    pos = draw(st.integers(min_value=0))
    strand = draw(st.sampled_from("+-"))
    mapq = draw(st.integers(min_value=0, max_value=60))
    just_a_number = draw(st.integers(min_value=0))
    cigar_lens = draw(st.lists(st.integers(min_value=1), min_size=3, max_size=100))
    if len(cigar_lens) % 2 == 0:
        # Need odd number of lengths
        cigar_lens.append(draw(st.integers(min_value=1)))
    n_cigar_ops = len(cigar_lens)

    clip_type = draw(st.sampled_from("SH"))
    cigar_ops = [f"{cigar_lens[0]}{clip_type}", f"{cigar_lens[1]}M"]

    for idx in range(2, n_cigar_ops - 1, 2):
        cig_op = draw(st.sampled_from("DI"))
        cigar_ops.append(f"{cigar_lens[idx]}{cig_op}")
        cigar_ops.append(f"{cigar_lens[idx+1]}M")

    cigar_ops.append(f"{cigar_lens[-1]}{clip_type}")
    cigar = "".join(cigar_ops)
    return f"{chrom},{pos},{strand},{cigar},{mapq},{just_a_number}"


# @given(
#     MaxSize=st.one_of(st.just(-1), st.integers(min_value=10)),
#     SV_size=st.integers(min_value=2),
#     Supplementary_info=st.lists(supplementary_info(), min_size=1, max_size=15),
#     max_split_parts=st.integers(min_value=1, max_value=10),
#     min_mapq=st.integers(min_value=0, max_value=60),
#     primary_info=st.tuples(
#         st.integers(min_value=0, max_value=100000),
#         st.integers(min_value=1, max_value=100000),
#         st.integers(min_value=0, max_value=int(1e9)),
#         st.integers(min_value=1, max_value=int(1e8)),
#         st.sampled_from(chroms),
#         st.sampled_from("+-"),
#     ).map(lambda x: (x[0], x[0] + x[1], x[2], x[2] + x[3], x[4], x[5])),
#     query=st.text(alphabet=set("ACGT"), min_size=10),
#     read_name=st.text(min_size=1, max_size=50),
#     total_L=st.integers(min_value=1),
# )
@given(
    MaxSize=st.just(-1),
    SV_size=st.integers(min_value=2),
    Supplementary_info=st.lists(supplementary_info(), min_size=1, max_size=15),
    max_split_parts=st.integers(min_value=1, max_value=10),
    min_mapq=st.integers(min_value=0, max_value=60),
    primary_info=st.tuples(
        st.integers(min_value=0, max_value=100000),
        st.integers(min_value=1, max_value=100000),
        st.integers(min_value=0, max_value=int(1e9)),
        st.integers(min_value=1, max_value=int(1e8)),
        st.sampled_from(chroms),
        st.sampled_from("+-"),
    ).map(lambda x: (x[0], x[0] + x[1], x[2], x[2] + x[3], x[4], x[5])),
    read_name=st.just("read-name"),
    total_L=st.integers(min_value=1),
)
# @settings(verbosity=Verbosity.verbose)
def test_equivalent_organize_split_signal_organize_split_signal_cuddly(
    MaxSize,
    SV_size,
    Supplementary_info,
    max_split_parts,
    min_mapq,
    primary_info,
    read_name,
    total_L,
):
    query = "N" * 100
    candidate = list()
    candidate_cuddly = list()
    result_organize_split_signal = cuteSV.cuddlySV.organize_split_signal(
        primary_info=list(primary_info),
        Supplementary_info=Supplementary_info.copy(),
        total_L=total_L,
        SV_size=SV_size,
        min_mapq=min_mapq,
        max_split_parts=max_split_parts,
        read_name=read_name,
        candidate=candidate,
        MaxSize=MaxSize,
        query=query,
    )
    result_organize_split_signal_cuddly = cuteSV.split_signal.organize_split_signal(
        primary_info=list(primary_info),
        Supplementary_info=Supplementary_info.copy(),
        total_L=total_L,
        SV_size=SV_size,
        min_mapq=min_mapq,
        max_split_parts=max_split_parts,
        read_name=read_name,
        candidate=candidate_cuddly,
        MaxSize=MaxSize,
        query=query,
    )
    assert result_organize_split_signal == result_organize_split_signal_cuddly, (
        result_organize_split_signal,
        result_organize_split_signal_cuddly,
    )
    candidate = {tuple(x) for x in candidate}
    candidate_cuddly = {tuple(x) for x in candidate_cuddly}
    assert candidate.issubset(candidate_cuddly), (candidate_cuddly, candidate)
    # assert len(candidate_cuddly - candidate) == 0, candidate_cuddly - candidate
    if not len(candidate_cuddly - candidate) == 0:
        print(candidate_cuddly - candidate)


@mark.skip("Not ready for this")
@given(
    MaxSize=st.one_of(st.just(-1), st.integers(min_value=10)),
    SV_size=st.integers(min_value=2),
    Supplementary_info=st.lists(supplementary_info(), min_size=1, max_size=15),
    max_split_parts=st.integers(min_value=1, max_value=10),
    min_mapq=st.integers(min_value=0, max_value=60),
    primary_info=st.tuples(
        st.integers(min_value=0),
        st.integers(min_value=1),
        st.integers(min_value=0),
        st.integers(min_value=1),
        st.sampled_from(chroms),
        st.sampled_from("+-"),
    ).map(lambda x: (x[0], x[0] + x[1], x[2], x[2] + x[3], x[4], x[5])),
    query=st.text(alphabet=set("ACGT"), min_size=10),
    read_name=st.text(min_size=1, max_size=50),
    total_L=st.integers(min_value=1),
)
@settings(verbosity=Verbosity.verbose)
def test_cuddly_candidate_out(
    MaxSize,
    SV_size,
    Supplementary_info,
    max_split_parts,
    min_mapq,
    primary_info,
    query,
    read_name,
    total_L,
):
    candidate_cuddly = list()
    cuteSV.split_signal.organize_split_signal(
        primary_info=list(primary_info),
        Supplementary_info=Supplementary_info.copy(),
        total_L=total_L,
        SV_size=SV_size,
        min_mapq=min_mapq,
        max_split_parts=max_split_parts,
        read_name=read_name,
        candidate=candidate_cuddly,
        MaxSize=MaxSize,
        query=query,
    )

    for i, cand in enumerate(candidate_cuddly):
        assert cand[3] in ("DUP", "INS", "DEL", "INV"), (i, cand)
        assert cand[4] in chroms, cand

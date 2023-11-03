#!/bin/bash

#  Created on Thursday, 19 October  2023

# For debugging
#set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

OUTPATH="panel_of_normals_wrk/"
usage() {
    echo -e "usage:
$0 -o OUTPATH/ -i input.tsv

-o OUTPATH/
-i input.tsv   Tab separated list of normal sample \"workdir/\\tcalls_with_RNAMES.vcf\"
-h          Show this message and exit." >&2
    exit 1

}

FIFTHofMEM=1G
while getopts "o:hn:i:" flag; do
    case "$flag" in
    o)
        OUTPATH="$OPTARG"
        ;;
    i)
        INPUTFILE="$OPTARG"
        ;;
    h | *)
        usage
        ;;
    esac
done
shift $((OPTIND - 1))
OPTIND=1
OUTPATH=$(readlink -f "$OUTPATH")
test ! -e "${OUTPATH}" || (
    echo "File ${OUTPATH} exist. Will not do anything!"
    exit 1
)
# Make temporary directory which will be cleaned at script exit
TEMPDIR=$(mktemp --directory -p "$(dirname "$OUTPATH")")
function _cleanup {
    rm -r "$TEMPDIR"
}
trap _cleanup EXIT

cat "$INPUTFILE" | while read -r NORMALPATH NORMALVCF; do
    echo Running $NORMALPATH $NORMALVCF
    bcftools query -f "%RNAMES\n" "${NORMALVCF}" | tr "," '\n' | sort -u >"${TEMPDIR}/_reads.lst"
    test -s "${TEMPDIR}/_reads.lst" || (
        echo "Couldn't find RNAMES from ${NORMALVCF}!"
        exit 1
    )
    for SIG_TYPE in DEL INS INV TRA DUP; do
        grep -Ff "${TEMPDIR}/_reads.lst" "${NORMALPATH}/${SIG_TYPE}.sigs" >>"${TEMPDIR}/${SIG_TYPE}.sigs" &
    done
    wait
done
mkdir "${OUTPATH}"

# Merge (Can't remove duplicates because CuteSV)

sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} "${TEMPDIR}/DEL.sigs" >"${OUTPATH}/DEL.sigs" &

sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} "${TEMPDIR}/INS.sigs" >"${OUTPATH}/INS.sigs" &

sort -k 2,2 -k 3,3 -k 4,5n --unique -S ${FIFTHofMEM} "${TEMPDIR}/INV.sigs" >"${OUTPATH}/INV.sigs" &

sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n -k 6,6n --unique -S ${FIFTHofMEM} "${TEMPDIR}/TRA.sigs" >"${OUTPATH}/TRA.sigs" &

sort -k 1,1r -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} "${TEMPDIR}/DUP.sigs" >"${OUTPATH}/DUP.sigs" &

wait

echo "Verifying output!"
test -s "${OUTPATH}/DUP.sigs"
test -s "${OUTPATH}/TRA.sigs"
test -s "${OUTPATH}/INV.sigs"
test -s "${OUTPATH}/INS.sigs"
test -s "${OUTPATH}/DEL.sigs"
echo "Looks fine!."

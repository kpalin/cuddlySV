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

usage() {
    echo -e "usage:
$0 -o OUTPATH/ -i INPATH/ 


Subset OUTPATH panel of normals to INPATH selecting only samples listed in samples.tsv

-o OUTPATH/
-i INPATH/
-s samples.tsv samples to subset 
-h          Show this message and exit." >&2
    exit 1

}

FIFTHofMEM=1G
while getopts "o:hi:s:" flag; do
    case "$flag" in
    o)
        OUTPATH="$OPTARG"
        ;;
    i)
        INPATH="$OPTARG"
        ;;
    s)
        SAMPLES="$OPTARG"
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

test -s "${SAMPLES}" || (
    echo "File ${SAMPLES} does not exist or is empty. Will not do anything!"
    exit 1
)

mkdir "${OUTPATH}"

# Merge. Remove duplicates for PoN

grep -Fwf ${SAMPLES} "${INPATH}/DEL.sigs" | sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${OUTPATH}/DEL.sigs" &

grep -Fwf ${SAMPLES} "${INPATH}/INS.sigs" | sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${OUTPATH}/INS.sigs" &

grep -Fwf ${SAMPLES} "${INPATH}/INV.sigs" | sort -k 2,2 -k 3,3 -k 4,5n --unique -S ${FIFTHofMEM} >"${OUTPATH}/INV.sigs" &

grep -Fwf ${SAMPLES} "${INPATH}/TRA.sigs" | sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n -k 6,6n --unique -S ${FIFTHofMEM} >"${OUTPATH}/TRA.sigs" &

grep -Fwf ${SAMPLES} "${INPATH}/DUP.sigs" | sort -k 1,1r -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${OUTPATH}/DUP.sigs" &

wait

echo "Verifying output!"
test -s "${OUTPATH}/DUP.sigs"
test -s "${OUTPATH}/TRA.sigs"
test -s "${OUTPATH}/INV.sigs"
test -s "${OUTPATH}/INS.sigs"
test -s "${OUTPATH}/DEL.sigs"
echo "Looks fine!."

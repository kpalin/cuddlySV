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

OUTPATH="merged_work/"

usage() {
    echo -e "usage:
$0 -o OUTPATH/ -t TUMORINPATH1/ -n PONINPATH/ [MORE_IN_PATH..]

-o OUTPATH/
-t tumor_work_dir
-n panel_of_normals_dir
-h          Show this message and exit." >&2
    exit 1

}
FIFTHofMEM=1G
while getopts "o:ht:n:" flag; do
    case "$flag" in
    o)
        OUTPATH="$OPTARG"
        ;;
    t)
        TUMORPATH="$OPTARG"
        ;;
    n)
        NORMALPATHS=("$OPTARG" "${@:$OPTIND}")
        ;;
    h | *)
        usage
        ;;
    esac
done
shift $((OPTIND - 1))
OPTIND=1
OUTPATH=$(readlink -f "$OUTPATH")
test ! -e "${OUTPATH}"
# Make temporary directory which will be cleaned at script exit
TEMPDIR=$(mktemp --directory -p "$(dirname "$OUTPATH")")
function _cleanup {
    rm -r "$TEMPDIR"
}
trap _cleanup EXIT

# Can not remove duplicate reads since we need to keep the 'normal' information

echo Will merge following paths: "${TUMORPATH}" "${NORMALPATHS[@]}"

sort -k 2,2 -k 3,4n -S ${FIFTHofMEM} "${TUMORPATH}/DEL.sigs" "${NORMALPATHS[@]/%/\/DEL.sigs}" >"${TEMPDIR}/DEL.sigs" &

sort -k 2,2 -k 3,4n -S ${FIFTHofMEM} "${TUMORPATH}/INS.sigs" "${NORMALPATHS[@]/%/\/INS.sigs}" >"${TEMPDIR}/INS.sigs" &

sort -k 2,2 -k 3,3 -k 4,5n -S ${FIFTHofMEM} "${TUMORPATH}/INV.sigs" "${NORMALPATHS[@]/%/\/INV.sigs}" >"${TEMPDIR}/INV.sigs" &

sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n -k 6,6n -S ${FIFTHofMEM} "${TUMORPATH}/TRA.sigs" "${NORMALPATHS[@]/%/\/TRA.sigs}" >"${TEMPDIR}/TRA.sigs" &

sort -k 1,1r -k 2,2 -k 3,4n -S ${FIFTHofMEM} "${TUMORPATH}/DUP.sigs" "${NORMALPATHS[@]/%/\/DUP.sigs}" >"${TEMPDIR}/DUP.sigs" &

wait
ln -s "${TUMORPATH}/reads.sigs" "${TEMPDIR}/"

test -s "${TEMPDIR}/DUP.sigs"
test -s "${TEMPDIR}/TRA.sigs"
test -s "${TEMPDIR}/INV.sigs"
test -s "${TEMPDIR}/INS.sigs"
test -s "${TEMPDIR}/DEL.sigs"

mv "${TEMPDIR}" "${OUTPATH}"
trap EXIT

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
$0 -o OUTPATH/  -n INPATH2/ [INPATHx/..]

-o OUTPATH/
-n normal_work_dir normal2_work_dir ..
-h          Show this message and exit." >&2
    exit 1

}
declare -a NORMALPATHS
FIFTHofMEM=1G
while getopts "o:hn:" flag; do
    case "$flag" in
    o)
        OUTPATH="$OPTARG"
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

# Merge, but remove duplicates
find "${NORMALPATHS[@]}" -maxdepth 1 -name DEL.sigs -print0 |
    xargs -0 sed -e 's/\b[a-zA-Z0-9_-]\+:/PanelOfNormals:/' |
    sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${TEMPDIR}/DEL.sigs" &

find "${NORMALPATHS[@]}" -maxdepth 1 -name INS.sigs -print0 |
    xargs -0 sed -e 's/\b[a-zA-Z0-9_-]\+:/PanelOfNormals:/' |
    sort -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${TEMPDIR}/INS.sigs" &

find "${NORMALPATHS[@]}" -maxdepth 1 -name INV.sigs -print0 |
    xargs -0 sed -e 's/\b[a-zA-Z0-9_-]\+:/PanelOfNormals:/' |
    sort -k 2,2 -k 3,3 -k 4,5n --unique -S ${FIFTHofMEM} >"${TEMPDIR}/INV.sigs" &

find "${NORMALPATHS[@]}" -maxdepth 1 -name TRA.sigs -print0 |
    xargs -0 sed -e 's/\b[a-zA-Z0-9_-]\+:/PanelOfNormals:/' |
    sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n -k 6,6n --unique -S ${FIFTHofMEM} >"${TEMPDIR}/TRA.sigs" &

find "${NORMALPATHS[@]}" -maxdepth 1 -name DUP.sigs -print0 |
    xargs -0 sed -e 's/\b[a-zA-Z0-9_-]\+:/PanelOfNormals:/' |
    sort -k 1,1r -k 2,2 -k 3,4n --unique -S ${FIFTHofMEM} >"${TEMPDIR}/DUP.sigs" &

wait

echo "Verifying output!"
test -s "${TEMPDIR}/DUP.sigs"
test -s "${TEMPDIR}/TRA.sigs"
test -s "${TEMPDIR}/INV.sigs"
test -s "${TEMPDIR}/INS.sigs"
test -s "${TEMPDIR}/DEL.sigs"

mv "${TEMPDIR}" "${OUTPATH}"
trap EXIT

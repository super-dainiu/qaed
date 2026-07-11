#!/bin/bash
# Download QTFM 3.4 (Quaternion Toolbox for MATLAB, Sangwine & Le Bihan) into
# qtfm/ and apply our patch (general non-Hermitian support in @quaternion/eig
# via eigq). Run once from the repository root.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

[ -d qtfm ] && { echo "qtfm/ already exists; remove it first."; exit 1; }

curl -L -o qtfm.zip "https://downloads.sourceforge.net/project/qtfm/qtfm/3.4/qtfm_3_4.zip"
unzip -q qtfm.zip
rm -f qtfm.zip
[ -d qtfm ] || mv qtfm_* qtfm

patch -p1 < patches/qtfm-nonhermitian-eig.patch
echo "qtfm 3.4 installed and patched."

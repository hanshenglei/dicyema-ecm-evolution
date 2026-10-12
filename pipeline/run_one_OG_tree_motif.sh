#!/usr/bin/env bash
set -euo pipefail

OG_DIR=$1
THREADS=${2:-10}
MEME_NMOTIFS=${3:-10}
MEME_MAXSIZE=${4:-20000000}
RUN_MEME=${5:-yes}

og=$(basename "${OG_DIR}")
faa="${OG_DIR}/${og}.faa"
log="${OG_DIR}/${og}.pipeline.log"

{
echo "[$(date '+%F %T')] Start ${og}"
echo "OG_DIR=${OG_DIR}"
echo "THREADS=${THREADS}"
echo "RUN_MEME=${RUN_MEME}"

if [ ! -s "${faa}" ]; then
    echo "ERROR: missing fasta ${faa}"
    exit 1
fi

nseq=$(grep -c "^>" "${faa}" || true)
echo "nseq=${nseq}"

if [ "${nseq}" -lt 3 ]; then
    echo "Skip ${og}: less than 3 sequences"
    exit 0
fi

# 1. MAFFT
if [ ! -s "${OG_DIR}/${og}.mafft.fa" ]; then
    echo "[$(date '+%F %T')] MAFFT ${og}"
    mafft --auto --thread "${THREADS}" "${faa}" > "${OG_DIR}/${og}.mafft.fa"
else
    echo "MAFFT exists, skip"
fi

# 2. trimAl
if [ ! -s "${OG_DIR}/${og}.mafft.trim.fa" ]; then
    echo "[$(date '+%F %T')] trimAl ${og}"
    trimal \
      -in "${OG_DIR}/${og}.mafft.fa" \
      -out "${OG_DIR}/${og}.mafft.trim.fa" \
      -automated1
else
    echo "trimAl exists, skip"
fi
# 如果 trimAl 没有生成有效结果，则使用原始 MAFFT 比对
if [ ! -s "${OG_DIR}/${og}.mafft.trim.fa" ]; then
    echo "WARNING: trimAl did not produce output for ${og}; use raw MAFFT alignment instead"
    cp "${OG_DIR}/${og}.mafft.fa" "${OG_DIR}/${og}.mafft.trim.fa"
    echo "untrimmed_mafft_used" > "${OG_DIR}/alignment_trim.note.txt"
fi

# 3. IQ-TREE
if [ ! -s "${OG_DIR}/${og}.mafft.trim.fa.treefile" ]; then
    echo "[$(date '+%F %T')] IQ-TREE ${og}"
    iqtree2 \
      -s "${OG_DIR}/${og}.mafft.trim.fa" \
      -m MFP \
      -B 1000 \
      -alrt 1000 \
      -T "${THREADS}" \
      --prefix "${OG_DIR}/${og}.mafft.trim.fa"
else
    echo "IQ-TREE treefile exists, skip"
fi

# 4. MEME
if [ "${RUN_MEME}" = "yes" ]; then
    if [ ! -s "${OG_DIR}/meme/meme.xml" ]; then
        echo "[$(date '+%F %T')] MEME ${og}"
        meme "${faa}" \
          -oc "${OG_DIR}/meme" \
          -protein \
          -mod zoops \
          -nmotifs "${MEME_NMOTIFS}" \
          -minw 6 \
          -maxw 50 \
          -maxsize "${MEME_MAXSIZE}" \
          -p "${THREADS}"
    else
        echo "MEME exists, skip"
    fi
else
    echo "Skip MEME by RUN_MEME=${RUN_MEME}"
fi

echo "[$(date '+%F %T')] Done ${og}"

} > "${log}" 2>&1

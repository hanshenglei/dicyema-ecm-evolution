#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=${1:-10.family_analysis/selected_OGs}
THREADS=${2:-16}
MEME_NMOTIFS=${3:-10}
MEME_MAXSIZE=${4:-20000000}
MAX_MEME_SEQ=${5:-100}

echo "ROOT_DIR=${ROOT_DIR}"
echo "THREADS=${THREADS}"
echo "MEME_NMOTIFS=${MEME_NMOTIFS}"
echo "MEME_MAXSIZE=${MEME_MAXSIZE}"
echo "MAX_MEME_SEQ=${MAX_MEME_SEQ}"

for faa in ${ROOT_DIR}/OG*/OG*.faa
do
    og=$(basename "$faa" .faa)
    dir=$(dirname "$faa")

    nseq=$(grep -c "^>" "$faa" || true)

    echo "[$(date '+%F %T')] Processing ${og}, nseq=${nseq}"

    if [ "$nseq" -lt 3 ]; then
        echo "Skip ${og}: less than 3 sequences"
        continue
    fi

    # 1. MAFFT
    if [ ! -s "${dir}/${og}.mafft.fa" ]; then
        mafft --auto --thread ${THREADS} "$faa" > "${dir}/${og}.mafft.fa"
    fi

    # 2. trimAl
    if [ ! -s "${dir}/${og}.mafft.trim.fa" ]; then
        trimal \
          -in "${dir}/${og}.mafft.fa" \
          -out "${dir}/${og}.mafft.trim.fa" \
          -automated1
    fi

    # 3. IQ-TREE
    if [ ! -s "${dir}/${og}.mafft.trim.fa.treefile" ]; then
        iqtree2 \
          -s "${dir}/${og}.mafft.trim.fa" \
          -m MFP \
          -B 1000 \
          -alrt 1000 \
          -T AUTO \
          --prefix "${dir}/${og}.mafft.trim.fa"
    fi

    # 4. MEME motif
    # 大 OG 不强行跑 MEME，避免卡死
    if [ "$nseq" -le "$MAX_MEME_SEQ" ]; then
        if [ ! -s "${dir}/meme/meme.xml" ]; then
            meme "$faa" \
              -oc "${dir}/meme" \
              -protein \
              -mod zoops \
              -nmotifs ${MEME_NMOTIFS} \
              -minw 6 \
              -maxw 50 \
              -maxsize ${MEME_MAXSIZE}
        fi
    else
        echo "Skip MEME for ${og}: nseq=${nseq} > MAX_MEME_SEQ=${MAX_MEME_SEQ}" \
          > "${dir}/meme.skipped.txt"
    fi

done

echo "Done."

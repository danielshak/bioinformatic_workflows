RSEMind=/home/dshak/cse185_sp19_team33/RSEM/RSEM
BAMDIR=/home/dshak/cse185_sp19_team33/RNA-seq/STAR_aligned_out/txBams/ERR315477Aligned.toTranscriptome.out.bam

rsem-calculate-expression \
    -p 8 \
    --paired-end \
    --fragment-length-mean -1 \
    --seed-length 25 \
    --bam ${BAMDIR} \
    ${RSEMind} \
    /home/dshak/cse185_sp19_team33/RNA-seq
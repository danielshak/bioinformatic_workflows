OUTDIR=/home/dshak/cse185_sp19_team33/RNA-seq/STAR_aligned_out
STARDIR=/home/dshak/cse185_sp19_team33/STAR
FASTQ=/home/dshak/cse185_sp19_team33/RNA-seq/fastqs/brain/ERR315477_from29
FQ1=ERR315477_1.fastq.gz
FQ2=ERR315477_2.fastq.gz

STAROPTS="--outSAMattributes NH HI AS NM MD \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --limitBAMsortRAM 50000000000"

STAR \
--runThreadN 24 \
--genomeDir ${STARDIR} \
--readFilesIn ${FASTQ}/${FQ1} ${FASTQ}/${FQ2} \
--readFilesCommand gunzip -c \
--outFileNamePrefix ${OUTDIR}/ERR315477 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM ${STAROPTS}

# Reorganize the output files
#mv ${OUTDIR}/${f}Aligned.toTranscriptome.out.bam ${OUTDIR}/txBams/
#mv ${OUTDIR}/${f}Aligned.sortedByCoord.out.bam ${OUTDIR}/genomeBams/
#samtools index ${OUTDIR}/genomeBams/${f}Aligned.sortedByCoord.out.bam
    

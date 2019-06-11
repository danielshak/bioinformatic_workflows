REFFA=/home/dshak/cse185_sp19_team33/Genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
GTF=/home/dshak/cse185_sp19_team33/Genomes/Homo_sapiens.GRCh38.83.gtf
RSEMind=/home/dshak/cse185_sp19_team33/RSEM/RSEM

rsem-prepare-reference \
    --gtf ${GTF} ${REFFA} ${RSEMind}


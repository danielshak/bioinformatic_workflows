REF=/home/dshak/cse185_sp19_team33/Genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
GTF=/home/dshak/cse185_sp19_team33/Genomes/Homo_sapiens.GRCh38.83.gtf
OUT=/home/dshak/cse185_sp19_team33/STAR

STAR \
--genomeSAsparseD 2 \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir ${OUT} \
--genomeFastaFiles ${REF} \
--sjdbGTFfile ${GTF} \
--sjdbOverhang 100
for file in ERR315432 ERR315455 ERR315477

do
~/.aspera/connect/bin/ascp -TQ -l200m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR315/${file} /datasets/home/96/596/dshak/ascp_out/rna_fastqs/brain
done
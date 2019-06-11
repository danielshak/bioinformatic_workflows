#!/bin/bash

set -e

TARGETVCF=$1
REFVCF=$2

# Fix 23andme VCF file to make alleles uppercase
echo "Fixing 23andMe VCF file to have uppercase alleles"
prefix=$(basename $TARGETVCF .vcf.gz)
zcat ${TARGETVCF} | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/g/G/g' | sed 's/c/C/g' | bgzip -c > ~/week3/${prefix}_fixcase.vcf.gz 
tabix -p vcf ~/week3/${prefix}_fixcase.vcf.gz 

# Get list of variants we will extract using plink
echo "Getting list of variants to extract"
zcat ${REFVCF} | grep -v "^#" | cut -f 1-5 | uniq > ~/week3/gtdata_alleles.tab
zcat ${REFVCF} | grep -v "^#" | cut -f 3 | uniq > ~/week3/var_ids.txt

# Run plink to align reference alleles across the two VCFs. Keep only variants present in our reference
echo "Align reference alleles and extract variants"
plink --vcf ~/week3/${prefix}_fixcase.vcf.gz  --recode vcf bgz --a2-allele ~/week3/gtdata_alleles.tab 4 3 '#' --real-ref-alleles --out ~/week3/${prefix}_corr  --extract ~/week3/var_ids.txt
tabix -p vcf ~/week3/${prefix}_corr.vcf.gz

# Finally, merge them! (and remove missing genotypes while we're at it)
echo "Merge VCFs to ~/week3/gtdata_merged.vcf.gz"
bcftools merge ~/week3/${prefix}_corr.vcf.gz ${REFVCF} | grep -v "\./\." | bgzip -c > ~/week3/gtdata_merged.vcf.gz  
tabix -p vcf ~/week3/gtdata_merged.vcf.gz 

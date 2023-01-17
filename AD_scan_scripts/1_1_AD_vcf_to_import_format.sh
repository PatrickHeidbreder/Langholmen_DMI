

#####
##### 0. Load Modules and set up directories 
#####

module load biokit


#####
##### 1. Extract AIM sites and only males from vcf 
#####

# AIMs were previously identified by Pierre Nouhaud. Where AIM = biallilic SNPs with an allele frequency difference of >90% between species.

# Male samples are used 

gunzip /scratch/project_2003480/allele_counts_all_parents_noFin_AllScaffolds_biall_125percNA.tsv.gz | cut -f 1-2 | sed 's/CHR/CHROM/' > AIMs_2col_format.table
gzip AIMs_2col_format.table


VCF=/scratch/project_2003480/patrick/data/vcf/filtered/final_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.MminDP3FminDP6.indMiss45.siteMiss40.vcf.gz
male_list=
bcftools view -R AIMs_2col_format.table.gz -S $VCF > 




#####
##### 2. Converting Reference allels to ancestral alleles ()
#####

# Filter based on last column (abs.diff.al) > 90% or 95% ( I think 95 is used in Li et al.)





#####
##### 1. 
#####
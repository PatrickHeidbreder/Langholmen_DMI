

cd /scratch/project_2003480/patrick/data/gvfc/raw
# Combine VCFs
for file in *noDiploidMales_filtered.vcf.gz ; do tabix $file ; done

ls *noDiploidMales_filtered.vcf.gz > all_filtered_vcfs.list

bcftools concat \
--allow-overlaps \
-f all_filtered_vcfs.list \
-Oz > ../concat_vcf/allnoDiploidMales_filtered.vcf.gz

tabix ../concat_vcf/allnoDiploidMales_filtered.vcf.gz



###
### PIXY COMPUTATION ------------------------------------------------------------------------------
###


# create population list
cd /scratch/project_2003480/patrick/analyses/pixy
awk '$6 == "f" {print $0}' ../../../sample_table_R_ordered_LanFix.tab | cut -f1,15 | grep -v 'FA07_5q' > pop.list


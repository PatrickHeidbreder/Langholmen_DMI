# additional filtering to be done.


# now that we have the final AIMs in the genome,
# we need to filter sites by LD to make sure sites are spread across the genome (following methods of Li et al 2022)
# then filter by LD (Li et al.)

module load plink/1.90

cd /scratch/project_2003480/patrick/analyses/ADscan/output

thresh=0.9
echo $thresh > name.tmp
thresh_name=$(sed 's/\.//' name.tmp)
rm name.tmp

PED_ORIG=painted_genotypes_males_95percAFD_wholegenome.ped
MAP_ORIG=painted_genotypes_males_95percAFD_wholegenome.map

plink --ped $PED_ORIG --map $MAP_ORIG --make-founders --indep-pairwise 10kb 1 $thresh --allow-extra-chr --out painted_genotypes_males_95percAFD_wholegenome_${thresh_name}LDthinned


# Load modules
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2003480/tmp" >> ~/.Renviron


# arg order:
# dir must match prune.in file dir
# ped full path # set above
# map full path # set above
# prune.in
DIR=$(pwd)
PRUNE=painted_genotypes_males_95percAFD_wholegenome_${thresh_name}LDthinned.prune.in

# Run the R script

RSCRIPT=/scratch/project_2003480/patrick/Langholmen_DMI/AD_scan_scripts/pruneinsites_ped_map.R

Rscript --no-save $RSCRIPT $DIR $PED_ORIG $MAP_ORIG $PRUNE 

### END ---


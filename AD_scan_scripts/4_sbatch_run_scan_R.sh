#!/bin/bash -l
#SBATCH --job-name=r_multithread_run_AD
#SBATCH --account=project_2003480
#SBATCH --output=/scratch/project_2003480/patrick/analyses/ADscan/logs/AD_scan_%j.out
#SBATCH --error=/scratch/project_2003480/patrick/analyses/ADscan/logs/AD_scan_%j.err
#SBATCH --partition=small
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=400M
#SBATCH --array=1

###
#### Step 0. Load modules, set directories ----
###
echo Loading modules, setting R environment ...

cd /scratch/project_2003480/patrick/analyses/ADscan/output

# Load modules
module load r-env

module load plink/1.90

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2003480/tmp" >> ~/.Renviron

# Match thread and core numbers
export APPTAINERENV_OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Thread affinity control
export APPTAINERENV_OMP_PLACES=cores
export APPTAINERENV_OMP_PROC_BIND=close


###
#### Step 1. Subsample sites ----
###

### thin sites to 700 making new .ped and .map files for input to AD scan.
# set input files
MAP_ORIG=/scratch/project_2003480/patrick/analyses/ADscan/output/painted_genotypes_males_95percAFD_wholegenome_09LDthinned.map
PED_ORIG=/scratch/project_2003480/patrick/analyses/ADscan/output/painted_genotypes_males_95percAFD_wholegenome_09LDthinned.ped

#SLURM_ARRAY_TASK_ID=1

cd subsample700/

# set number of sites to subsample
SUBSAMPLE=700

echo Preparing to subsample $SUBSAMPLE sites ...

# randomly remove positions until only 700 are left
plink --map $MAP_ORIG --ped $PED_ORIG --threads $SLURM_CPUS_PER_TASK --allow-extra-chr --make-bed --thin-count ${SUBSAMPLE} --out painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}

cut -f 2 painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.bim > painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.prune.in
rm cut -f 2 painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.bed
rm cut -f 2 painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.bim
rm cut -f 2 painted_genotypes_males_95percAFD_wholegenome_09LDthinned_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.fam


# use the files output above to extrat the random 700 positions fromt he original .ped and .map files.
DIR=$(pwd)
PRUNE=painted_genotypes_males_95percAFD_wholegenome_09LDthinned_700subsample_${SLURM_ARRAY_TASK_ID}.prune.in

## Run the R script to make the subsampled .ped
# arg order:
# dir must match prune.in file dir
# ped full path # set above
# map full path # set above
# prune.in

RSCRIPT=/scratch/project_2003480/patrick/Langholmen_DMI/AD_scan_scripts/pruneinsites_ped_map.R

srun apptainer_wrapper exec Rscript --no-save $RSCRIPT $DIR $PED_ORIG $MAP_ORIG $PRUNE 


###
#### Step 2. Run AD scan. ----
###

echo Running AD scan on bootstrap replicate: ${SLURM_ARRAY_TASK_ID}

# Set input files
MAP_SUB=${DIR}/painted_genotypes_males_95percAFD_wholegenome_09LDthinned_700subsample_${SLURM_ARRAY_TASK_ID}.map
PED_SUB=${DIR}/painted_genotypes_males_95percAFD_wholegenome_09LDthinned_700subsample_${SLURM_ARRAY_TASK_ID}.ped

# Run the R script
RSCRIPT=/scratch/project_2003480/patrick/Langholmen_DMI/AD_scan_scripts/fishStat.permutation.R

cd ../out/
srun apptainer_wrapper exec Rscript --no-save $RSCRIPT $MAP_SUB $PED_SUB AD_bootstrap_${SUBSAMPLE}subsample_${SLURM_ARRAY_TASK_ID}.out

# prefix=$output/5/TLMC_chr20
# Rscript $code/fishStat.permutation.R $prefix.map $prefix.ped $prefix.out

### END ---- 
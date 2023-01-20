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
#SBATCH --mem-per-cpu=50G
#SBATCH --array=2

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


###
#### Step 2. Run AD scan. ----
###

echo Running AD scan on bootstrap replicate: ${SLURM_ARRAY_TASK_ID}

# Set input files
MAP_SUB=/scratch/project_2003480/patrick/analyses/ADscan/output/painted_genotypes_males_95percAFD_wholegenome_09LDthinned.map
PED_SUB=/scratch/project_2003480/patrick/analyses/ADscan/output/painted_genotypes_males_95percAFD_wholegenome_09LDthinned.ped

# Run the R script
RSCRIPT=/scratch/project_2003480/patrick/Langholmen_DMI/AD_scan_scripts/fishStat.permutation.R

cd ../out/
srun apptainer_wrapper exec Rscript --no-save $RSCRIPT $MAP_SUB $PED_SUB AD_wholegenome_bootstrap_${SLURM_ARRAY_TASK_ID}.out

# prefix=$output/5/TLMC_chr20
# Rscript $code/fishStat.permutation.R $prefix.map $prefix.ped $prefix.out

### END ---- 
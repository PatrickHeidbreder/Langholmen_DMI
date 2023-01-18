#!/bin/bash -l
#SBATCH --job-name=r_multithread_run_AD
#SBATCH --account=project_2003480
#SBATCH --output=/scratch/project_2003480/patrick/analyses/ADscan/logs/AD_scan_%j.out
#SBATCH --error=/scratch/project_2003480/patrick/analyses/ADscan/logs/AD_scan_%j.err
#SBATCH --partition=hugemem
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=150G
#SBATCH --array=1

cd /scratch/project_2003480/patrick/analyses/ADscan/output/out

# Load r-env
module load r-env

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


# Set input files

MAP=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../map.list)
PED=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../ped.list)
SCAFFOLD=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2003480/patrick/reference_genome/scaffold.list)


# Run the R script

RSCRIPT=/scratch/project_2003480/patrick/Langholmen_DMI/AD_scan_scripts/fishStat.permutation.R

srun apptainer_wrapper exec Rscript --no-save $RSCRIPT $MAP $PED ${SCAFFOLD}_AD.out


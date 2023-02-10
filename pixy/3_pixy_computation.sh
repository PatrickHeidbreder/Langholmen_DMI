
#!/bin/bash -l
#SBATCH -J pixy
#SBATCH -o /scratch/project_2003480/patrick/analyses/pixy/logs/pixy2.out
#SBATCH -e /scratch/project_2003480/patrick/analyses/pixy/logs/pixy2.err
#SBATCH --account=project_2003480
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/pixy

pixy --stats pi fst dxy \
--vcf /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw/allFems_filtered.vcf.gz \
--populations pop.list \
--bed_file windows_20kb.bed \
--output_prefix allFems \
--n_cores 8


# sort output - FST
head -n1 allFems_fst.txt > tmp
grep -v 'pop' allFems_fst.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_fst_sort.txt

# sort output - Dxy
head -n1 allFems_dxy.txt > tmp
grep -v 'pop' allFems_dxy.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_dxy_sort.txt

# sort output - Pi
head -n1 allFems_pi.txt > tmp
grep -v 'pop' allFems_pi.txt | sort -k2,2 -k3,3n > tmp2
cat tmp tmp2 > allFems_pi_sort.txt

rm tmp tmp2
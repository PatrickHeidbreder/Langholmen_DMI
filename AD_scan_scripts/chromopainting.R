# chromosome_painting_genomewide.R
setwd("/scratch/project_2003480/patrick/analyses/ADscan")
options(stringsAsFactors = F)
library(stringr)

#### Read data ----

# metadata - you need it to determine the ploidy (sex) of your individuals 
# smpl = read.table("sample_table_R_ordered_LanFix.tab", h=T, comment.char = "", as.is = T)
smpl = read.table("/scratch/project_2003480/patrick/2022_batch_metadata_350IDfiltered_withLineages.tsv", h=T, comment.char = "", as.is = T)
# geno file, assuming the order of the inds in the geno file is the same as in the metadata table
geno = read.table("your.geno.gz", h = T, comment.char = "")
#geno = geno[,-which(names(geno) == "OUT")] # here I remove the outgroup
names(geno)[1] = "CHROM"

# counts and allele freq differences between the two species  
parent.cnts = read.table(paste0("/scratch/project_2003480/allele_counts_all_parents_noFin_AllScaffolds_biall_125percNA.tsv.gz"), h=T)


#### Extract aquilonia alleles ----

# We are only dealing with aquilonia alleles at AIMs
# Sites are strictly biallelic so if they are not aq, they must be pol - two sides of the same coin

# Here we define the allele frequency difference threshold used to define your AIMs - only sites >= to thresh will be used afterwards
thresh = .95

print(paste0("### Identifying AIMs with frequencies differing from at least ", thresh*100, "%..."))

# Create a vector scaffold_position for all AIMs
pos.aim = paste0(parent.cnts$CHR[which(round(parent.cnts$abs.diff.a1, 2) >= thresh)], 
                 "_",
                 parent.cnts$POS[which(round(parent.cnts$abs.diff.a1, 2) >= thresh)])

# keep only AIMs
aq.alleles.aim = parent.cnts[which(round(parent.cnts$abs.diff.a1, 2) >= thresh),]




# extract a vector containing aquilonia alleles for each site
cov = aq.alleles.aim[,5:6]
cov.max = apply(cov, 1, which.max)
al = aq.alleles.aim[,3:4]
aq.alleles.aim = NULL
for(i in 1:nrow(al)){
  aq.alleles.aim = c(aq.alleles.aim, al[i, cov.max[i]])
}

head(aq.alleles.aim)


#### Prep genotypes for hybrids & parents ----


print(paste0("### Prepare hybrid genotypes for painting..."))

# Subset genotypes at AIMs
tmp.pos = paste0(geno$CHROM, "_", geno$POS)
dat = geno[which(tmp.pos %in% pos.aim), ]

# subset aq AIMs in positions that re present in the data

tmp.pos = paste0(dat$CHROM, "_", dat$POS)
pos.aim <- pos.aim[which(pos.aim %in% tmp.pos)]

aq.alleles.aim <- aq.alleles.aim[which(pos.aim %in% tmp.pos)]

# Split genotypes and diploidize males (IIRC males are already diploid in your geno, in which case you can omit the if statement)
# hyb = dat[, -c(1:2)] # remove CHR & POS cols
# hyb.chr1 = hyb.chr2 = matrix(data = NA, ncol = ncol(hyb), nrow = nrow(hyb))
# for (i in 1:ncol(hyb)){  # we loop over individuals - cols in geno, rows in metadata
#   if(smpl_table$sex[which(smpl_table$id == names(hyb)[i])] == "f"){ # female samples
#     hyb.chr1[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,1]
#     hyb.chr2[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,2]
#   }else{ # male samples: the genotype is duplicated 
#     hyb.chr1[, i] = hyb.chr2[, i] = hyb[,i]
#   }
# }

# Split genotypes and diploidize males (IIRC males are already diploid in your geno, in which case you can omit the if statement)
hyb = dat[, -c(1:2)] # remove CHR & POS cols

# check that all sites inthe genotype matrix are homozygous

head(hyb)

hyb.chr1 = hyb.chr2 = matrix(data = NA, ncol = ncol(hyb), nrow = nrow(hyb))
for (i in 1:ncol(hyb)){  # we loop over individuals - cols in geno, rows in metadata
  hyb.chr1[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,1]
  hyb.chr2[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,2]
}

colnames(hyb.chr1) = colnames(hyb.chr2) = names(hyb)


hom.check <- hyb.chr1==hyb.chr2
sum(!hom.check)

# recode het sites as missing data

hyb.chr1[hyb.chr1!=hyb.chr2] <- "N"
hyb.chr2[hyb.chr1!=hyb.chr2] <- "N"

hom.check <- hyb.chr1==hyb.chr2
sum(!hom.check)


# Tag missing data (Ns in geno are missing)
hyb.chr1[hyb.chr1 == "N"] = NA
hyb.chr2[hyb.chr2 == "N"] = NA

# Male filtering (deprecated)
hyb.chr1.noNAmales = hyb.chr1 #[-which(rowSums(is.na(hyb.chr1))==40),]
hyb.chr2.noNAmales = hyb.chr2 #[-which(rowSums(is.na(hyb.chr1))==40),] # want to remove the same sites
aq.alleles.aim.noNAmales = aq.alleles.aim #[-which(rowSums(is.na(hyb.chr1))==40)]


pos.aim.noNAmales = pos.aim #[-which(rowS <- ums(is.na(hyb.chr1))==40)]


#### Chromosome painting ----

print(paste0("### Painting..."))

# real.pos.aim.noNAmales = as.numeric(matrix(ncol = 2, byrow = T, 
#                                            data = unlist(strsplit(x = pos.aim.noNAmales, split = "_")))[, 2]) # extract AIM positions

real.pos.aim.noNAmales = matrix(ncol = 2, byrow = T, 
                                           data = unlist(strsplit(x = pos.aim.noNAmales, split = "_"))) # extract AIM positions

tmp1 = hyb.chr1.noNAmales == aq.alleles.aim.noNAmales # this is for the first chromosome copy, gives TRUE if aq, FALSE otherwise (pol)
tmp2 = hyb.chr2.noNAmales == aq.alleles.aim.noNAmales # same for the second chromosome copy, gives TRUE if aq, FALSE otherwise (pol)
paint.dat = tmp1+tmp2 # we add the two copies: TRUE + TRUE = 2 (HOM aq), FALSE + FALSE = 0 (HOM pol), TRUE + FALSE = 1 (HET)

# save painting data:  2 = HOM aq, 0 = HOM pol, 1 = HET
writedat = cbind(real.pos.aim.noNAmales, paint.dat)
colnames(writedat)[1:2] = c("CHR", "POS")

setwd("/scratch/project_2003480/patrick/analyses/ADscan/output/painted_matrix_perscaff")
for(my.scaffold in unique(writedat[,1])){
write.table(x = writedat[which(writedat[,1] == my.scaffold),], 
            file = paste0("painted_genotypes_males_", thresh*100, "percAFD_", my.scaffold, ".tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")
}

table(writedat[,-c(1,2)])
# recode to follow input for ped file for AD method
# 0 = pol
# 1 = aq
# NA = missing data

setwd("/scratch/project_2003480/patrick/analyses/ADscan/output/painted_matrix_perscaff")
# this for loop changes each data matrix to ped format, by recoding alleles as above
# it then adds 6 columns in front of each row (.ped format documentation)
# it then duplicates each allele

makerow.ped <- function(x){
  # add mandatory columns, 1-6
  mand.col <- NULL
  ID <- paste0(x[1],"_",x[2])
  mand.col[1:6] <- c(ID, ID, -9, -9, 1, -9)
  
  # add dupliated alleles
  allele.single <- x[-c(1,2)]
  allele.dup <- rep(allele.single, each = 2)
  final.row <- c(mand.col,allele.dup)
  return(final.row)
}




setwd("/scratch/project_2003480/patrick/analyses/ADscan/output/ped")
for(my.scaffold in unique(writedat[,1])){
  # recode alleles
  forped <- writedat[which(writedat[,1] == my.scaffold),]
  forped[which(forped == 2)] <- 1
  
  
  # put into .ped format with makerow function
  ped <- apply(forped, 1, FUN = makerow.ped)
  ped <- t(ped)
  ped <- unname(ped)
  
  write.table(x = ped, 
              file = paste0("painted_genotypes_males_", thresh*100, "percAFD_", my.scaffold, ".ped"), 
              col.names = F, row.names = F, quote = F, sep = "\t")
}



makerow.map <- function(x){
  # add mandatory columns, 1-6
  mand.col <- NULL
  ID <- paste0(x[1],"_",x[2])
  mand.col[1:4] <- c(x[1], ID, -9, x[2])
}


setwd("/scratch/project_2003480/patrick/analyses/ADscan/output/map")
for(my.scaffold in unique(writedat[,1])){
  formap <- writedat[which(writedat[,1] == my.scaffold),]
  # put into .map format with makerow.map function
  map <- apply(formap, 1, FUN = makerow.map)
  map <- t(map)
  map <- unname(map)
  
  write.table(x = map, 
              file = paste0("painted_genotypes_males_", thresh*100, "percAFD_", my.scaffold, ".map"), 
              col.names = F, row.names = F, quote = F, sep = "\t")
}

#### END ----
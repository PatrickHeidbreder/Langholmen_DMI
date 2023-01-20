# chromosome_painting_genomewide.R

# This script takes a .geno formated file and outputs a painted genotype matrix based
# on the AIMs identified by Pierre from aq and pol samples.
#
# The painted genotype matrix is also converted to a .ped and .map file pair
# for input to PLINK.

setwd("/scratch/project_2003480/patrick/analyses/ADscan")
options(stringsAsFactors = F)
library(stringr)

###
#### Step 1. Read data ----
###
print(paste0("### Reading data ..."))

# metadata - you need it to determine the ploidy (sex) of your individuals 
# smpl = read.table("sample_table_R_ordered_LanFix.tab", h=T, comment.char = "", as.is = T)
smpl = read.table("/scratch/project_2003480/patrick/2022_batch_metadata_350IDfiltered_withLineages.tsv", h=T, comment.char = "", as.is = T)
# geno file, assuming the order of the inds in the geno file is the same as in the metadata table
geno = read.table("your.geno.gz", h = T, comment.char = "")
#geno = geno[,-which(names(geno) == "OUT")] # here I remove the outgroup
names(geno)[1] = "CHROM"

# counts and allele freq differences between the two species  
parent.cnts = read.table(paste0("/scratch/project_2003480/allele_counts_all_parents_noFin_AllScaffolds_biall_125percNA.tsv.gz"), h=T)

###
#### Step 2. Extract aquilonia alleles ----
###


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

###
#### Step 3. Prep genotypes for hybrids & parents ----
###

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


#### Step 4. Check for het sites ----

print(paste0("### Checking for het sites in gentype matrix ..."))
# check that all sites in the genotype matrix are homozygous

# head(hyb)

hyb.chr1 = hyb.chr2 = matrix(data = NA, ncol = ncol(hyb), nrow = nrow(hyb))
for (i in 1:ncol(hyb)){  # we loop over individuals - cols in geno, rows in metadata
  hyb.chr1[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,1]
  hyb.chr2[, i] = str_split_fixed(string = hyb[,i], pattern = "/", n = 2)[,2]
}

colnames(hyb.chr1) = colnames(hyb.chr2) = names(hyb)


hom.check <- hyb.chr1==hyb.chr2
print(paste0("Found ",sum(!hom.check), " sites with het genotypes."))


# recode het sites as missing data

hyb.chr1[hyb.chr1!=hyb.chr2] <- "N"
hyb.chr2[hyb.chr1!=hyb.chr2] <- "N"

hom.check <- hyb.chr1==hyb.chr2
print(paste0("Found ",sum(!hom.check), " het sites after recoding to missing."))


# Tag missing data (Ns in geno are missing)
hyb.chr1[hyb.chr1 == "N"] = NA
hyb.chr2[hyb.chr2 == "N"] = NA

# Male filtering (deprecated)
hyb.chr1.noNAmales = hyb.chr1 #[-which(rowSums(is.na(hyb.chr1))==40),]
hyb.chr2.noNAmales = hyb.chr2 #[-which(rowSums(is.na(hyb.chr1))==40),] # want to remove the same sites
aq.alleles.aim.noNAmales = aq.alleles.aim #[-which(rowSums(is.na(hyb.chr1))==40)]


pos.aim.noNAmales = pos.aim #[-which(rowS <- ums(is.na(hyb.chr1))==40)]

###
#### Step 5. Chromosome painting ----
###

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

setwd("/scratch/project_2003480/patrick/analyses/ADscan/output")
write.table(x = writedat, 
            file = paste0("painted_genotypes_males_", thresh*100, "percAFD_wholegenome.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")
###
#### Step 6. Output .ped and .map files ----
###

print(paste0("### Creating .ped and .map files..."))


# This script converts the painted genotype matrix into .ped and .map files
# for input to plink

setwd("/scratch/project_2003480/patrick/analyses/ADscan/output/")


###
##### Step 6.1. Convert to .ped format ----
###

# .ped files should have 6 meta data columns at the front
# metadata is followed by genotypes with two values each (one per chrom)
# columns 7+ should be coded as 012 with 0 for missing data
#
# .ped fil should be space delimited

# this for loop changes each data matrix row to ped format
# it adds 6 columns in front of each row (.ped format documentation)
# it then duplicates each allele

makerow.ped <- function(x){
  # add mandatory columns, 1-6
  mand.col <- NULL
  ID <- names(x)
  mand.col[1:6] <- c(ID, ID, -9, -9, 1, -9)
  
  # add dupliated alleles
  allele.single <- x[-c(1,2)]
  allele.dup <- rep(allele.single, each = 2)
  final.row <- c(mand.col,allele.dup)
  return(final.row)
}


# matrix coding for the painted genotype matrix:
# 0 = pol
# 1 = het # there are no het sites as they were removed
# 2 = aq
# N = missing

# recode for plink LD filtering:
# 0 = missing
# 1 = pol
# 2 = aq


# recode alleles int he painted matrix
forped <- writedat
# select matrix and recode to ped format
forped.values <- as.matrix(forped[1:nrow(forped),3:ncol(forped)])
forped.values[forped.values == 0] <- 1 # change pol code from 0 to 1
forped.values[is.na(forped.values)] <- 0 # change missing values from NA to 0

forped[1:nrow(forped),3:ncol(forped)] <- forped.values
forped <- forped[,-c(1,2)]

ped <- matrix(, nrow = 246, ncol = nrow(forped)*2+6)

# add mandatory columns, 1-6
for(i in 1:ncol(forped)){
  mand.col <- NULL
  ID <- colnames(forped)[i]
  mand.col[1:6] <- c(ID, ID, -9, -9, 1, -9)
  
  # add dupliated alleles
  allele.single <- forped[,i]
  allele.dup <- rep(allele.single, each = 2)
  final.row <- c(mand.col,allele.dup)
  ped[i,] <- final.row
}


# put patined matrix into .ped format with makerow function

ped <- unname(ped)

write.table(x = ped, 
            file = paste0("painted_genotypes_males_", thresh*100, "percAFD_wholegenome.ped"), 
            col.names = F, row.names = F, quote = F, sep = ' ')

###
##### Step 6.2. Convert to .map file ----
###

# .map files are companions to .ped files. They simply carry metadata about
# sites and samples

# there are 4 columns.
# chrom name , chrom name_site pos , -9 (missing), site pos
#
# columns are tab delimited

## old code that was to parallelize across scaffolds, this idea will nto work as
## it will not look for interchromosomal DMIs
# for(my.scaffold in unique(writedat[,1])){
#   formap <- writedat[which(writedat[,1] == my.scaffold),]
#   # put into .map format with makerow.map function
#   map <- apply(formap, 1, FUN = makerow.map)
#   map <- t(map)
#   map <- unname(map)
# 
#   write.table(x = map,
#               file = paste0("painted_genotypes_males_", thresh*100, "percAFD_", my.scaffold, ".map"),
#               col.names = F, row.names = F, quote = F, sep = "\t")
# }

## make the .map file rows
makerow.map <- function(x){
  # add mandatory columns, 1-6
  mand.col <- NULL
  ID <- paste0(str_replace_all(x[1], " ", ""),"_",str_replace_all(x[2], " ", ""))
  mand.col[1:4] <- c(x[1], ID, -9, x[2])
}

formap <- writedat
# put into .map format with makerow.map function
map <- apply(formap, 1, FUN = makerow.map)
map <- t(map)
map <- unname(map)

write.table(x = map,
            file = paste0("painted_genotypes_males_", thresh*100, "percAFD_wholegenome.map"),
            col.names = F, row.names = F, quote = F, sep = "\t", fileEncoding = "UTF-8")

### END ----

## This script filters sites .ped and .map files based on a 1 column CHR_pos text file (prune.in)
# for input to plink
args=commandArgs(trailingOnly = TRUE)

###
#### Step 1. Load in args ----
###

# arg order:
# dir must match prune.in file dir
# ped full path
# map full path
# prune.in


# load data

setwd(args[1])

ped <- read.table(args[2], header = F, sep = ' ')

map <- read.table(args[3], header = F, sep = "\t")

prune.in <- scan(args[4],
                 what = "character")

###
#### Step 2. Filter sites in .ped file ----
###

ped.filter <- ped[,c(rep(TRUE,6),rep(map[,2] %in% prune.in, each = 2))]

# extract new name for ped/map file name with regular expression (based on prune.in)
ped.map.name <- sub("\\..*", "", args[4])

write.table(x = ped.filter, 
            file = paste0(ped.map.name,".ped"), 
            col.names = F, row.names = F, quote = F, sep = ' ')

###
#### Step 3. Filter sites in the .map file ----
###

map.filter <- map[map[,2] %in% prune.in,]

write.table(x = map.filter,
            file = paste0(ped.map.name,".map"),
            col.names = F, row.names = F, quote = F, sep = "\t")


#### END ----
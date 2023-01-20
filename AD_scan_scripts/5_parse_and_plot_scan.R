library(stringr)

scan.data <- read.csv("/scratch/project_2003480/patrick/analyses/ADscan/output/out/AD_bootstrap_700subsample_1.out",
                        header = T,
                        sep = "\t")

 colnames(scan.data)[1] <- "Pos1"

# Parse positions into separate chr and position columns.
 Pos1_chr <- lapply(str_split(scan.data$Pos1, "_"), function(x){
   return(x[1])})
 scan.data$Pos1_chr <- unlist(Pos1_chr)
 
 Pos1_pos <- lapply(str_split(scan.data$Pos1, "_"), function(x){
   return(x[2])})
 scan.data$Pos1_pos <- as.integer(unlist(Pos1_pos))
 
 
 Pos2_chr <- lapply(str_split(scan.data$Pos2, "_"), function(x){
   return(x[1])})
 scan.data$Pos2_chr <- unlist(Pos2_chr) 
 
 Pos2_pos <- lapply(str_split(scan.data$Pos2, "_"), function(x){
   return(x[2])})
 scan.data$Pos2_pos <- as.integer(unlist(Pos2_pos)) 
 
 
 
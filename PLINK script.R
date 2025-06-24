# This script runs PLINK analysis on output files from Axiom Analysis Suite. 
# Prior to running PLINK, putative triploid samples are removed, and SNP chromosome locations are added

#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Genotyping.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Remove triploids from .ped file
ids_to_remove <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_remove <- ids_to_remove$V1
ped <- ped[!(ped$V1 %in% ids_to_remove), ]

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#Load .map file and BLAST results
map <- read.csv("JD_PFR_Genotyping.map", header = FALSE, sep ="\t")
BLAST <- read.csv("BLAST results.tsv", header = FALSE, sep = "\t")

#Match Marker IDs
match_ids <- match(map$V2, BLAST$V2)
matched <- !is.na(match_ids)
map[matched, ] <- BLAST[match_ids[matched], ]

#Remove header row, and set any chromosome numbers larger than 17 to 0
map = map[-1, ]
map$V1 <- as.numeric(as.character(map$V1))
map[map$V1 > 17, c("V1", "V4")] <- 0

#Save .map file (with locations)
write.table(map, "JD_PFR_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##RUNNING PLINK

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/PLINK Matching")

#Create subset of SNPs in PLINK
system("plink --file JD_PFR_PLINK --no-fid --no-parents --no-sex --no-pheno --missing-genotype 0 --geno 0.01 --maf 0.3 --indep 50 5 1.5")


#Run PLINK dup matching
system("plink --file JD_PFR_PLINK --no-fid --no-parents --no-sex --no-pheno --extract plink.prune.in --missing-genotype 0 --genome full")


#Convert plink.genome to a tab-delimited file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "plink_tab_delimited1.txt", sep = "\t", row.names = FALSE, quote = FALSE)




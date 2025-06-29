# This script runs PLINK analysis on output files from Axiom Analysis Suite. 
# Prior to running PLINK, putative triploid samples are removed, and SNP chromosome locations are added

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/PLINK Duplicate Identification/Inputs")


#.ped file curation

#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Genotyping.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Remove triploids from .ped file
ids_to_remove <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_remove <- ids_to_remove$V1
ped <- ped[!(ped$V1 %in% ids_to_remove), ]

#Add empty rows for PLINK formatting
ped <- add_column(ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
ped <- add_column(ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#.map file curation

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


##Running PLINK

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/PLINK Duplicate Identification/Inputs")


#Run PLINK
system("plink --file JD_PFR_PLINK --missing-genotype 0 --genome full")

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/PLINK Duplicate Identification/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Filter for PI_HAT >0.97 (duplicate threshold)
genome <- genome[!(genome$PI_HAT < 0.97), ]
genome <- subset(genome, select = c("IID1","IID2"))

##Grouping and graphing duplicates

#Group duplicates with igraph
graph <- graph_from_data_frame(genome, directed = FALSE)
components <- components(graph)

#Sort groupings by number of duplicates
group_sizes <- table(components$membership)
sorted_group_ids <- order(group_sizes)
new_ids <- match(components$membership, sorted_group_ids)
V(graph)$group <- new_ids
grouped_samples <- split(names(components$membership), new_ids)

#Pad group with length less than max length with NA's
max_len <- max(sapply(grouped_samples, length))
padded_list <- lapply(grouped_samples, function(x) {c(x, rep(" ", max_len - length(x)))})

#Write groupings to a dataframe
dd <- as.data.frame(do.call(rbind, padded_list))

#Add a number for each group
dd <- cbind(Group = seq_len(nrow(dd)), dd)

# Add the number of duplicates in each grouping
sample_counts <- rowSums(dd[, -1] != " ")
dd <- add_column(dd, SampleCount = sample_counts, .after = "Group")

#Rename columns
colnames(dd) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13")

#Save .csv of duplicate groupings
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/PLINK Duplicate Identification/Grouped_Duplicates.csv", row.names = FALSE)

#Graph duplicates

#Assign a group ID
V(graph)$group <- components$membership

#Assign a color to each group
group_colors <- rainbow(length(unique(V(graph)$group)))
V(graph)$color <- group_colors[V(graph)$group]

#Plot and save as .png
png(filename = "C:/Users/curly/Desktop/Apple Genotyping/Results/PLINK Duplicate Identification/PLINK Clonal Groups.png", width = 1000, height = 1000, pointsize = 20)
plot(graph, layout = layout_with_fr(graph), vertex.size = 3, vertex.label = NA, main = "PLINK Clonal Groups")

dev.off()




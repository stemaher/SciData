library(PhosR)
library(dplyr)
# Import MS data
PhosphoSites <- read.table(file = "Input/STY_79.9663.tsv", header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")

rownames(PhosphoSites) <- paste0(PhosphoSites$Protein,";",
                                           PhosphoSites$Gene,";",
                                           sapply(strsplit(PhosphoSites$Index, "_"), "[", 2) ,";",
                                           PhosphoSites$Peptide)

PhosphoSites <- PhosphoSites[,grepl("*MaxLFQ.Intensity*", colnames(PhosphoSites))]
PhosphoSites[PhosphoSites=="0"]<-NA
olnames(PhosphoSites) <- gsub("X220724kw_SMaher","",colnames(PhosphoSites))
colnames(PhosphoSites) <- gsub(".MaxLFQ.Intensity","",colnames(PhosphoSites))

# delete rows with all NAs  
PhosphoSites <- PhosphoSites[rowSums(is.na(PhosphoSites)) != ncol(PhosphoSites), ]


# Experimental design for the phospho data frame 
experimental_design <- read.table(file = "/Input/experimantal_design_phospho.csv", header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")
colnames(experimental_design) <- c("label","condition", "Cell_Line","Timepoint","Ligand","Biological.Replicate","Technical.Replicate")
experimental_design$label <- sub("-",".", experimental_design$label)
experimental_design$label <- sub(" ",".", experimental_design$label)
experimental_design$label <- sub("Intensity.","", experimental_design$label)

experimental_design$condition <- paste0(experimental_design$Cell_Line,"_",experimental_design$Timepoint)
experimental_design$condition <- sub("/",".",experimental_design$condition)
experimental_design$condition <- sub("TrK","Trk",experimental_design$condition)
rownames(experimental_design) <- experimental_design$label
experimental_design$condition2 <- sapply(strsplit(experimental_design$Cell_Line, "/"), "[", 1)
experimental_design$Trk <- sapply(strsplit(experimental_design$Cell_Line, "/"), "[", 2)

# checking whether the Intensity colnames are matched with the Exp des labels
intersect(colnames(PhosphoSites),row.names(experimental_design))
experimental_design <- experimental_design[colnames(PhosphoSites),]
# should be TRUE
all(colnames(PhosphoSites) == row.names(experimental_design))


################################
### Preprocessing Phosposites###
################################

temp <- PhosphoSites
temp <- temp[rowSums(is.na(temp)) != ncol(temp), ]
temp <- log2(temp)
#temp[is.na(temp)] <- 0
temp2 <- preprocessCore::normalize.quantiles(as.matrix(temp))
rownames(temp2) <- rownames(temp)
colnames(temp2) <- colnames(temp)
temp2[temp2=="0"]<-NA
temp <- selectGrps(as.matrix(temp2), experimental_design$condition, 0.5, n=1)
rm(temp2)
boxplot(temp)

temp <- PhosR::scImpute(temp, 0.5, experimental_design$condition)
temp <- PhosR::tImpute(temp)

phospho <- temp

saveRDS(phospho, file = "NormalizedData.Rds")
saveRDS(experimental_design, file = "experimental_design.Rds")



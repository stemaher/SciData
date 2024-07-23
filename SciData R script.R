library(PhosR)
library(dplyr)

#Phosphoproteomics ####
##Load Normalised Data & Annotation file ####
phospho <- readRDS("/mnt/Data/RStudio/smaher/PhosphoproteomicsNov2022/NormalizedData.Rds")
experimental_design <- readRDS("/mnt/Data/RStudio/smaher/PhosphoproteomicsNov2022/experimental_design.Rds")


##PCA ####
pca <- prcomp(t(as.matrix(phospho)), scale. = TRUE, center = TRUE)
library("FactoMineR")
library("factoextra")
fviz_pca_ind(pca, repel = TRUE, habillage = experimental_design[,8], invisible="quali", axes = c(1, 2),
             legend.title = "Parental Cell Line", pointsize = 3, show.legend = FALSE, geom="point",
             ggtheme=theme(axis.text=element_text(size=50)) )+
  labs(title ="") +
  #    scale_color_brewer(palette="Set2") +
  theme_minimal()



##Differently expressed Phosphosites (DEPSs) ####
library(limma)
library(data.table)
library(splines)
datalimma <- list()

groups <- as.factor(experimental_design$condition)
levels(groups)
design <-  model.matrix(~0 + groups)
colnames(design) <- sub("groups","",colnames(design))
colnames(design) <- gsub("/",".",colnames(design) )
colnames(design) <- gsub(":",".",colnames(design) )

contrast.matrix <- makeContrasts(
  SY5Y.TrkA.10 = SY5Y.TrkA_10 - SY5Y.TrkA_0,
  SY5Y.TrkA.45 = SY5Y.TrkA_45 - SY5Y.TrkA_0,
  SY5Y.TrkA.24 = SY5Y.TrkA_24 - SY5Y.TrkA_0,
  SY5Y.TrkB.10 = SY5Y.TrkB_10 - SY5Y.TrkB_0,
  SY5Y.TrkB.45 = SY5Y.TrkB_45 - SY5Y.TrkB_0,
  SY5Y.TrkB.24 = SY5Y.TrkB_24 - SY5Y.TrkB_0,
  SY5Y.TrkC.10 = SY5Y.TrkC_10 - SY5Y.TrkC_0,
  SY5Y.TrkC.45 = SY5Y.TrkC_45 - SY5Y.TrkC_0,
  SY5Y.TrkC.24 = SY5Y.TrkC_24 - SY5Y.TrkC_0,
  
  NLF.TrkA.10 = NLF.TrkA_10 - NLF.TrkA_0,
  NLF.TrkA.45 = NLF.TrkA_45 - NLF.TrkA_0,
  NLF.TrkA.24 = NLF.TrkA_24 - NLF.TrkA_0,
  NLF.TrkB.10 = NLF.TrkB_10 - NLF.TrkB_0,
  NLF.TrkB.45 = NLF.TrkB_45 - NLF.TrkB_0,
  NLF.TrkB.24 = NLF.TrkB_24 - NLF.TrkB_0,
  NLF.TrkC.10 = NLF.TrkC_10 - NLF.TrkC_0,
  NLF.TrkC.45 = NLF.TrkC_45 - NLF.TrkC_0,
  NLF.TrkC.24 = NLF.TrkC_24 - NLF.TrkC_0,
  
  NBLS.TrkA.10 = NBLS.TrkA_10 - NBLS.TrkA_0,
  NBLS.TrkA.45 = NBLS.TrkA_45 - NBLS.TrkA_0,
  NBLS.TrkA.24 = NBLS.TrkA_24 - NBLS.TrkA_0,
  NBLS.TrkB.10 = NBLS.TrkB_10 - NBLS.TrkB_0,
  NBLS.TrkB.45 = NBLS.TrkB_45 - NBLS.TrkB_0,
  NBLS.TrkB.24 = NBLS.TrkB_24 - NBLS.TrkB_0,
  NBLS.TrkC.10 = NBLS.TrkC_10 - NBLS.TrkC_0,
  NBLS.TrkC.45 = NBLS.TrkC_45 - NBLS.TrkC_0,
  NBLS.TrkC.24 = NBLS.TrkC_24 - NBLS.TrkC_0,
  levels=design)

colnames(contrast.matrix)
fit <- lmFit(phospho, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

for (s in 1:length(colnames(contrast.matrix))) {
  datalimma[[s]] <- topTable(fit2, coef=s, n=Inf)
  setDT(datalimma[[s]], keep.rownames = TRUE)[]
  datalimma[[s]]$type <- colnames(contrast.matrix)[s]
}
DEP <- do.call(rbind, datalimma)
colnames(DEP)[1] <- "site"

DEP$gene <- paste0(sapply(strsplit(DEP$site, ";"), "[", 2),"_",
                   sapply(strsplit(DEP$site, ";"), "[", 3))

DEP$gene2 <- paste0(sapply(strsplit(DEP$site, ";"), "[", 2))

names(datalimma) <- colnames(contrast.matrix)

selection <- lapply(datalimma, subset, adj.P.Val < 0.05 &  abs(logFC) > log2(1.5)) # for all conditions
DEP_list <- vector("list", length(selection))
# Iterate over the dataframes
for (i in seq_along(selection)) {
  # Extract the first column of the dataframe
  DEP_list[[i]] <- selection[[i]]$rn
}
names(DEP_list) <- names(selection)


##Venn Diagram Upregulated DEPs only ####

# i) Select upregulated DEPS
upreg_selection <- lapply(datalimma, subset, adj.P.Val < 0.05 & (logFC) > log2(1.5)) # for all conditions
upreg_DEP_list <- vector("list", length(upreg_selection))
# Iterate over the dataframes
for (i in seq_along(upreg_selection)) {
  # Extract the first column of the dataframe
  upreg_DEP_list[[i]] <- upreg_selection[[i]]$rn
}
names(upreg_DEP_list) <- names(upreg_selection)


# ii) Draw Venn diagram
library(VennDiagram)
library(RColorBrewer)
colors <- brewer.pal(n = 3, name = "Set1")
display.brewer.all() 
VD <- venn.diagram(upreg_DEP_list[c(1,19,10)], filename=NULL, fill = colors,
                     cex = 2, # size for values
                     fontfamily = "Arial",
                     #fontface = "bold",
                     #cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
                     cat.cex = 0.8, # font size for name of category
                     #cat.pos = 0.3,
                     # cat.dist = 0.07,
                     cat.fontfamily = "Arial", margin = 0.05)
grid.newpage()
grid.draw(VD)

# Total Proteomics ####

##Load Normalised Data & Annotation file ####

library(PhosR)
library(dplyr)

Intensities 
experimental_design
rownames(experimental_design) <- experimental_design$LABEL


##PCA ####
pca <- prcomp(t(as.matrix(Intensities)), scale. = TRUE, center = TRUE)
library("FactoMineR")
library("factoextra")
fviz_pca_ind(pca, repel = TRUE, habillage = experimental_design[,8]
             , invisible="quali", axes = c(1, 2),
             legend.title = "Parental Cell Line", pointsize = 3, show.legend = FALSE, geom="point",
             ggtheme=theme(axis.text=element_text(size=50)) )+
  labs(title ="") +
  #    scale_color_brewer(palette="Set2") +
  theme_minimal()

##Differently Expressed Proteins ####

library(limma)
library(data.table)
library(splines)

datalimma <- list()
groups <- as.factor(experimental_design$Condition)
levels(groups)
design <-  model.matrix(~0 + groups)
colnames(design) <- sub("groups","",colnames(design))
colnames(design) <- gsub("/",".",colnames(design) )
colnames(design) <- gsub(":",".",colnames(design) )

contrast.matrix <- makeContrasts(  
  SY5Y.TrkA.10 = SY5Y.TrKA_10 - SY5Y.TrKA_0,
  SY5Y.TrkA.45 = SY5Y.TrKA_45 - SY5Y.TrKA_0,
  SY5Y.TrkA.24 = SY5Y.TrKA_24 - SY5Y.TrKA_0,
  SY5Y.TrkB.10 = SY5Y.TrkB_10 - SY5Y.TrkB_0,
  SY5Y.TrkB.45 = SY5Y.TrkB_45 - SY5Y.TrkB_0,
  SY5Y.TrkB.24 = SY5Y.TrkB_24 - SY5Y.TrkB_0,
  SY5Y.TrkC.10 = SY5Y.TrkC_10 - SY5Y.TrkC_0,
  SY5Y.TrkC.45 = SY5Y.TrkC_45 - SY5Y.TrkC_0,
  SY5Y.TrkC.24 = SY5Y.TrkC_24 - SY5Y.TrkC_0,
  
  NBLS.TrkA.10 = NBLS.TrKA_10 - NBLS.TrKA_0,
  NBLS.TrkA.45 = NBLS.TrKA_45 - NBLS.TrKA_0,
  NBLS.TrkA.24 = NBLS.TrKA_24 - NBLS.TrKA_0,
  NBLS.TrkB.10 = NBLS.TrkB_10 - NBLS.TrkB_0,
  NBLS.TrkB.45 = NBLS.TrkB_45 - NBLS.TrkB_0,
  NBLS.TrkB.24 = NBLS.TrkB_24 - NBLS.TrkB_0,
  NBLS.TrkC.10 = NBLS.TrkC_10 - NBLS.TrkC_0,
  NBLS.TrkC.45 = NBLS.TrkC_45 - NBLS.TrkC_0,
  NBLS.TrkC.24 = NBLS.TrkC_24 - NBLS.TrkC_0,
  
  NLF.TrkA.10 = NLF.TrKA_10 - NLF.TrKA_0,
  NLF.TrkA.45 = NLF.TrKA_45 - NLF.TrKA_0,
  NLF.TrkA.24 = NLF.TrKA_24 - NLF.TrKA_0,
  NLF.TrkB.10 = NLF.TrkB_10 - NLF.TrkB_0,
  NLF.TrkB.45 = NLF.TrkB_45 - NLF.TrkB_0,
  NLF.TrkB.24 = NLF.TrkB_24 - NLF.TrkB_0,
  NLF.TrkC.10 = NLF.TrkC_10 - NLF.TrkC_0,
  NLF.TrkC.45 = NLF.TrkC_45 - NLF.TrkC_0,
  NLF.TrkC.24 = NLF.TrkC_24 - NLF.TrkC_0,
  levels=design)

colnames(contrast.matrix)
fit <- lmFit(Intensities, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


DEP$gene <- paste0(sapply(strsplit(DEP$protein, "_"), "[", 2))

names(datalimma) <- colnames(contrast.matrix)

DEP_list <- lapply(datalimma, subset, adj.P.Val < 0.05 & abs(logFC) > log2(1)) # adj.p.val =  Benjamini-Hochberg procedure FDR
DEP_list <- lapply(DEP_list, "[", , rn)

library(readxl)
library(editData)
library(openxlsx)

#Load annotation.xlsx
if (file.exists("out/annotation_copy.xlsx")) {
  print("Reading annotation_copy.xlsx")
  annotation <- read_excel("out/annotation_copy.xlsx")
} else {
  annotation <- read_excel("out/annotation.xlsx")
}

# Use DE for checking top genes of some clusters
DE_genes <- read.csv("temp/DE_genes.csv", row.names = 1)

#Check top 10 genes of a certain cluster
clust = 30
head(DE_genes[DE_genes$cluster == clust, ], 10)

# Change annotation
annotation <- editData(annotation)

write.xlsx(annotation, "out/annotation.xlsx")
if(file.exists("out/annotation_copy.xlsx")) {file.remove("out/annotation_copy.xlsx")}
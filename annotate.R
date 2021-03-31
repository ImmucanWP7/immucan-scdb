library(readxl)
library(editData)
library(openxlsx)

if (file.exists("out/annotation_copy.xlsx")) {
  print("Reading annotation_copy.xlsx")
  annotation <- read_excel("out/annotation_copy.xlsx")
} else {
  annotation <- read_excel("out/annotation.xlsx")
}

annotation <- editData(annotation)

# Use DE for checking top genes of some clusters
DE_genes <- read.csv("temp/DE_genes.csv", row.names = TRUE)

#Check top 10 genes of a certain cluster
clust = 19
head(DE_genes[DE_genes$cluster == clust, ], 10)

write.xlsx(annotation, "out/annotation.xlsx")
if(file.exists("out/annotation_copy.xlsx")) {file.remove("out/annotation_copy.xlsx")}
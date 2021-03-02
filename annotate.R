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

write.xlsx(annotation, "out/annotation.xlsx")

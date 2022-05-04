library(pheatmap)
library(Seurat)
data <- read.csv(file = "RPPA_data_Trk_normalized.csv", sep = "\t", stringsAsFactors=FALSE)
rownames(data) <- data$Antibodies
data$Antibodies <- NULL
data2 <- data[,grepl("DMSO",colnames(data))]
data2 <- log2(data2) 
data2  <- MinMax(data2, min = -2.1, max = 2.1)
png("HM11092020.png", width = 4, height = 7, units = 'in', res = 600)
pheatmap(data2,fontsize_row = 3.5, fontsize_col = 4)
dev.off()


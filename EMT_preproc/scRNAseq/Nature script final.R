#################
#### EMT ########
#################
library(Seurat)
# downloading data and DEG calculation

# function for DEG calculation
calculateDEGs <- function(UMImatrix,metadata,name) {
  data <- CreateSeuratObject(counts = UMImatrix, 
                             meta.data = metadata)
  data <- NormalizeData(data)
  stat <- as.data.frame(table(metadata$Drug))
  stat <- subset(stat, Freq > 70) #minimum number of cells per condition
  DEGlist <- list()
  for (i in unique(stat$Var1)){
    temp <- FindMarkers(data, ident.1 = i, ident.2 = paste0(name), group.by = 'Drug',  test.use = "wilcox")
    if (nrow(temp) > 0) {
      temp <- subset(temp, p_val_adj < 0.05)
      temp$gene <- rownames(temp)
      DEGlist[[paste0(deparse(substitute(metadata)),"_",i)]] <- as.data.frame(temp)
    }
  }
  return(DEGlist)
}

Func.Normalize <- function(matrix, metadata) {
  data <- CreateSeuratObject(counts = matrix)
  data <- NormalizeData(data)
  data2 <- as.matrix(data@assays$RNA@data)
  return(data2)
}

DEGlistALL <- list()
# A549 EGF
destfile <- "GSE147405_A549_EGF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_EGF_KinaseScreen_UMI_matrix.csv.gz", destfile)
A549_EGF_matrix <- read.csv("GSE147405_A549_EGF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_A549_EGF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_EGF_KinaseScreen_metadata.csv.gz", destfile)
A549_EGF_metadata <- read.csv("GSE147405_A549_EGF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- calculateDEGs(A549_EGF_matrix, A549_EGF_metadata,"Uninhibited_EGF")

#
A549_EGF.normalized <- Func.Normalize(A549_EGF_matrix)
write.table(A549_EGF.normalized, "/normalized/A549_EGF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(A549_EGF_metadata, "/normalized/A549_EGF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

# A549 TGFB1
destfile <- "GSE147405_A549_TGFB1_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_TGFB1_KinaseScreen_UMI_matrix.csv.gz", destfile)
A549_TGFB1_matrix <- read.csv("GSE147405_A549_TGFB1_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_A549_TGFB1_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_TGFB1_KinaseScreen_metadata.csv.gz", destfile)
A549_TGFB1_metadata <- read.csv("GSE147405_A549_TGFB1_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(A549_TGFB1_matrix, A549_TGFB1_metadata, "Uninhibited_TGFB1"))
#
A549_TGFB1.normalized <- Func.Normalize(A549_TGFB1_matrix)
write.table(A549_TGFB1.normalized, "/normalized/A549_TGFB1.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(A549_TGFB1_metadata, "/normalized/A549_TGFB1_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)



# A549 TNF
destfile <- "GSE147405_A549_TNF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_TNF_KinaseScreen_UMI_matrix.csv.gz", destfile)
A549_TNF_matrix <- read.csv("GSE147405_A549_TNF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_A549_TNF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_A549_TNF_KinaseScreen_metadata.csv.gz", destfile)
A549_TNF_metadata <- read.csv("GSE147405_A549_TNF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(A549_TNF_matrix, A549_TNF_metadata,"Uninhibited_TNF"))
#
A549_TNF.normalized <- Func.Normalize(A549_TNF_matrix)
write.table(A549_TNF.normalized, "/normalized/A549_TNF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(A549_TNF_metadata, "/normalized/A549_TNF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)



#DU145_EGF
destfile <- "GSE147405_DU145_EGF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_EGF_KinaseScreen_UMI_matrix.csv.gz", destfile)
DU145_EGF_matrix <- read.csv("GSE147405_DU145_EGF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_DU145_EGF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_EGF_KinaseScreen_metadata.csv.gz", destfile)
DU145_EGF_metadata <- read.csv("GSE147405_DU145_EGF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(DU145_EGF_matrix, DU145_EGF_metadata,"Uninhibited_EGF"))
#
DU145_EGF.normalized <- Func.Normalize(DU145_EGF_matrix)
write.table(DU145_EGF.normalized, "/normalized/DU145_EGF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(DU145_EGF_metadata, "/normalized/DU145_EGF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)


#DU145_TGFB1
destfile <- "GSE147405_DU145_TGFB1_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_TGFB1_KinaseScreen_UMI_matrix.csv.gz", destfile)
DU145_TGFB1_matrix <- read.csv("GSE147405_DU145_TGFB1_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_DU145_TGFB1_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_TGFB1_KinaseScreen_metadata.csv.gz", destfile)
DU145_TGFB1_metadata <- read.csv("GSE147405_DU145_TGFB1_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(DU145_TGFB1_matrix, DU145_TGFB1_metadata, "Uninhibited_TGFB1"))

#
DU145_TGFB1.normalized <- Func.Normalize(DU145_TGFB1_matrix)
write.table(DU145_TGFB1.normalized, "/normalized/DU145_TGFB1.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(DU145_TGFB1_metadata, "/normalized/DU145_TGFB1_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

#DU145_TNF
destfile <- "GSE147405_DU145_TNF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_TNF_KinaseScreen_UMI_matrix.csv.gz", destfile)
DU145_TNF_matrix <- read.csv("GSE147405_DU145_TNF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_DU145_TNF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_DU145_TNF_KinaseScreen_metadata.csv.gz", destfile)
DU145_TNF_metadata <- read.csv("GSE147405_DU145_TNF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(DU145_TNF_matrix, DU145_TNF_metadata,"Uninhibited_TNF"))

#
DU145_TNF.normalized <- Func.Normalize(DU145_TNF_matrix)
write.table(DU145_TNF.normalized, "/normalized/DU145_TNF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(DU145_TNF_metadata, "/normalized/DU145_TNF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)


#MCF7_EGF
destfile <- "GSE147405_MCF7_EGF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_EGF_KinaseScreen_UMI_matrix.csv.gz", destfile)
MCF7_EGF_matrix <- read.csv("GSE147405_MCF7_EGF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_MCF7_EGF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_EGF_KinaseScreen_metadata.csv.gz", destfile)
MCF7_EGF_metadata <- read.csv("GSE147405_MCF7_EGF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(MCF7_EGF_matrix, MCF7_EGF_metadata, "Uninhibited_EGF"))

MCF7_EGF.normalized <- Func.Normalize(MCF7_EGF_matrix)
write.table(MCF7_EGF.normalized, "/normalized/MCF7_EGF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(MCF7_EGF_metadata, "/normalized/MCF7_EGF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

#MCF7_TGFB1
destfile <- "GSE147405_MCF7_TGFB1_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_TGFB1_KinaseScreen_UMI_matrix.csv.gz", destfile)
MCF7_TGFB1_matrix <- read.csv("GSE147405_MCF7_TGFB1_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_MCF7_TGFB1_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_TGFB1_KinaseScreen_metadata.csv.gz", destfile)
MCF7_TGFB1_metadata <- read.csv("GSE147405_MCF7_TGFB1_KinaseScreen_metadata.csv", row.names = "X")
DEGlistALL <- append(DEGlistALL, calculateDEGs(MCF7_TGFB1_matrix, MCF7_TGFB1_metadata,"Uninhibited_TGFB1"))

MCF7_TGFB1.normalized <- Func.Normalize(MCF7_TGFB1_matrix)
write.table(MCF7_TGFB1.normalized, "/normalized/MCF7_TGFB1.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(MCF7_TGFB1_metadata, "/normalized/MCF7_TGFB1_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

#MCF7_TNF
destfile <- "GSE147405_MCF7_TNF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_TNF_KinaseScreen_UMI_matrix.csv.gz", destfile)
MCF7_TNF_matrix <- read.csv("GSE147405_MCF7_TNF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_MCF7_TNF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_MCF7_TNF_KinaseScreen_metadata.csv.gz", destfile)
MCF7_TNF_metadata <- read.csv("GSE147405_MCF7_TNF_KinaseScreen_metadata.csv", row.names = "X")
DEGlistALL <- append(DEGlistALL, calculateDEGs(MCF7_TNF_matrix, MCF7_TNF_metadata, "Uninhibited_TNF"))

MCF7_TNF.normalized <- Func.Normalize(MCF7_TNF_matrix)
write.table(MCF7_TNF.normalized, "/normalized/MCF7_TNF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(MCF7_TNF_metadata, "/normalized/MCF7_TNF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)


#OVCA420_EGF
destfile <- "GSE147405_OVCA420_EGF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_EGF_KinaseScreen_UMI_matrix.csv.gz", destfile)
OVCA420_EGF_matrix <- read.csv("GSE147405_OVCA420_EGF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_OVCA420_EGF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_EGF_KinaseScreen_metadata.csv.gz", destfile)
OVCA420_EGF_metadata <- read.csv("GSE147405_OVCA420_EGF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(OVCA420_EGF_matrix, OVCA420_EGF_metadata, "Uninhibited_EGF"))
OVCA420_EGF.normalized <- Func.Normalize(OVCA420_EGF_matrix)
write.table(OVCA420_EGF.normalized, "/normalized/OVCA420_EGF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(OVCA420_EGF_metadata, "/normalized/OVCA420_EGF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

#OVCA420_TGFB1
destfile <- "GSE147405_OVCA420_TGFB1_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TGFB1_KinaseScreen_UMI_matrix.csv.gz", destfile)
OVCA420_TGFB1_matrix <- read.csv("GSE147405_OVCA420_TGFB1_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_OVCA420_TGFB1_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TGFB1_KinaseScreen_metadata.csv.gz", destfile)
OVCA420_TGFB1_metadata <- read.csv("GSE147405_OVCA420_TGFB1_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(OVCA420_TGFB1_matrix,
                                               OVCA420_TGFB1_metadata,
                                               "Uninhibited_TGFB1"))

OVCA420_TGFB1.normalized <- Func.Normalize(OVCA420_TGFB1_matrix)
write.table(OVCA420_TGFB1.normalized, "/normalized/OVCA420_TGFB1.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(OVCA420_TGFB1_metadata, "/normalized/OVCA420_TGFB1_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

#OVCA420_TNF
destfile <- "GSE147405_OVCA420_TNF_KinaseScreen_UMI_matrix.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TNF_KinaseScreen_UMI_matrix.csv.gz", destfile)
OVCA420_TNF_matrix <- read.csv("GSE147405_OVCA420_TNF_KinaseScreen_UMI_matrix.csv", row.names = "X")
#
destfile <- "GSE147405_OVCA420_TNF_KinaseScreen_metadata.csv"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TNF_KinaseScreen_metadata.csv.gz", destfile)
OVCA420_TNF_metadata <- read.csv("GSE147405_OVCA420_TNF_KinaseScreen_metadata.csv", row.names = "X")

DEGlistALL <- append(DEGlistALL, calculateDEGs(OVCA420_TNF_UMI_matrix,
                                               OVCA420_TNF_metadata,
                                               "Uninhibited_TNF"))
#
OVCA420_TNF.normalized <- Func.Normalize(OVCA420_TNF_matrix)
write.table(OVCA420_TNF.normalized, "/normalized/OVCA420_TNF.normalized.txt", row.names = T, sep = "\t", quote = F, col.names = NA)
write.table(OVCA420_TNF_metadata, "/normalized/OVCA420_TNF_metadata.txt", row.names = T, sep = "\t", quote = F, col.names = NA)

# save results for the KEA analysis
saveRDS(DEGlistALL, file = "DEGlistALL.Rds")
DEGlist.long <- dplyr::bind_rows(DEGlistALL, .id = "condition")
DEGlist.long$condition <- sub("metadata_","",DEGlist.long$condition)
write.table(DEGlist.long, file = "DEGlist.long.txt", sep = "\t", quote = F, row.names = F)

stat2 <- as.data.frame(sapply(DEGlistALL, nrow))
stat2 <- subset(stat2, stat2[1] >= 20) # minimum number of DEGs
for (i in unique(rownames(stat2))){
  write.table(DEGlistALL[[i]]$gene, file = paste0("/KEAinput/",i,".txt"), sep = "\t",row.names = F, quote = F)
}

# import the KEA results
KEAresults <- list()
for (i in list.files("/KEAresults")) {
  KEAresults[[i]] <- read.csv(file = paste0(i), sep = "\t")
}

KEAresults.long <- dplyr::bind_rows(KEAresults, .id = "condition")
KEAresults.long$condition <- gsub("metadata_","",KEAresults.long$condition)
KEAresults.long$condition <- gsub(" ","",KEAresults.long$condition)
KEAresults.long$condition <- gsub(".txt","",KEAresults.long$condition)
KEAresults.wide <- reshape(KEAresults.long[,c(1,5,6)], timevar = "condition", idvar = "name", direction = "wide")
colnames(KEAresults.wide) <- sub("pvalue.","",colnames(KEAresults.wide))

#save the lists
write.table(KEAresults.wide, file = "KEAresults.wide.txt", sep = "\t", quote = F, row.names = F)
write.table(KEAresults.long[,c(1,5,6,2,3,4,7)], file = "KEAresults.long.txt", sep = "\t", quote = F, row.names = F)




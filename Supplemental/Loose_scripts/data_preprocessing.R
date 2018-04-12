sample.data <- function(df, n){
  names <- sample(df$Run, n)
  r.df <- df[df$Run %in% names,]
  return(r.df)
}

#subsetting 50 cells per platform
gse48968 <- read.table("~/FYP/fastq/GSE48968_dendritic.txt", sep = "\t", header = T)
gse48968$Run <- as.character(gse48968$Run)
gse48968 <- gse48968[,c("Run", "Instrument")]
hs25 <- gse48968[gse48968$Instrument == "Illumina HiSeq 2500",]
hs20 <- gse48968[gse48968$Instrument == "Illumina HiSeq 2000",]
hs25 <- sample.data(hs25, 50)
hs20 <- sample.data(hs20, 50)
gse48968 <- rbind(hs20, hs25)
acc <- gse48968$Run
write(acc, file = "~/FYP/fastq/GSE48968_accessions.txt", sep = "\n")
write.table(gse48968, file = "~/FYP/fastq/GSE48968_dendritic.txt", row.names = FALSE, sep = "\t")

#subset by a certain category, then subset further and sample from each group for equal representation
gse59114 <- read.table("~/FYP/fastq/GSE59114_HSCs.txt", sep = "\t", header = T)
gse59114 <- gse59114[,c("Run", "Instrument", "cell_type")]
gse59114$Run <- as.character(gse59114$Run)
hs25 <- gse59114[gse59114$Instrument == "Illumina HiSeq 2500",]
hs25.lt <- hs25[hs25$cell_type == "long term hematopoietic stem cell",]
hs25.mpp <- hs25[hs25$cell_type == "multipotent progenitor",]
hs25.st <- hs25[hs25$cell_type == "short term hematopoietic stem cell",]
hs20 <- gse59114[gse59114$Instrument == "Illumina HiSeq 2000",]
hs20.lt <- hs20[hs20$cell_type == "long term hematopoietic stem cell",]
hs20.mpp <- hs20[hs20$cell_type == "multipotent progenitor",]
hs20.st <- hs20[hs20$cell_type == "short term hematopoietic stem cell",]
hs25.lt <- sample.data(hs25.lt,20)
hs25.mpp <- sample.data(hs25.mpp,20)
hs25.st <- sample.data(hs25.st,20)
hs20.lt <- sample.data(hs20.lt,20)
hs20.mpp <- sample.data(hs20.mpp,20)
hs20.st <- sample.data(hs20.st,20)
hs20 <- rbind(hs20.lt, hs20.mpp, hs20.st)
hs25 <- rbind(hs25.lt, hs25.mpp, hs25.st)
gse59114 <- rbind(hs20, hs25)
acc <- gse59114$Run
write(acc, file = "~/FYP/fastq/GSE59114_accessions.txt", sep = "\n")
write.table(gse59114, file = "~/FYP/fastq/GSE59114_HSCs.txt", row.names = FALSE, sep = "\t")

##script for merging cell expression estimates into dataset level tables

lymphoid <- list.files(path = "FYP/OUTPUT/homo_lymphoid", pattern = "genes.results",full.names = TRUE)
cortex <- list.files(path = "FYP/OUTPUT/mus_cortex/", pattern = "genes.results",full.names = TRUE)
embryo <- list.files(path = "FYP/OUTPUT/mus_embryo/", pattern = "genes.results",full.names = TRUE)
categories <- c("length", "effective_length", "expected_count", "TPM", "FPKM")

#### metatable function ####
metatable <- function(filenames, field, group){
    cell.ids <- gsub(".genes.results", "" , basename(filenames))
    var.table <- list()
    for(i in seq_along(filenames)){
        ival <- read.table(filenames[i], header = TRUE, stringsAsFactors = FALSE, row.names = 1)
        icol <- ival[,field]
        names(icol) <- rownames(ival)
        var.table <- c(var.table, list(icol))
    }
    var.table <- as.data.frame(var.table)
    colnames(var.table) <- cell.ids
    var.table$names <- rownames(var.table)
    rownames(var.table) <- NULL
    var.table$names <- gsub("^[^_]+_", "", var.table$names)
    gene.names <- var.table$names
    merged <- aggregate(var.table[names(var.table) != "names"], by = list(gene.names), FUN = "mean")
    rownames(merged) <- merged$Group.1
    merged <- merged[,-1]
    write.table(merged, paste(group, field, "table.txt", sep = "_"), sep="\t")
    return(merged[1:5,1:10])
}
####  End of function ####

system.time(
    for(i in seq_along(categories)){
        metatable(lymphoid, categories[i], "Hlymphoid")
        metatable(cortex, categories[i], "Mcortex")
        metatable(embryo, categories[i], "Membryo")
    }
)


## getting transcript abundances from fastq using kallisto
GSE70580.quant.time <- system.time(
  GSE70580.quant <- getKallistoQuant(sampleFile = "~/FYP/GSE70580/GSE70580_sampFile.txt",
                                    index = "~/FYP/ref/human.idx", resultDir = "~/FYP/kallisto_out/GSE70580/",
                                    outputPrefix = "GSE70580", nCores = 2,
                                    fragment_length = 200, fragment_standard_deviation = 20)
)
GSE70580 <- readKallistoResults(GSE70580.quant)
saveRDS(GSE70580, file = "~/FYP/GSE70580.rds")

GSE59114.quant.time <- system.time(
  GSE59114.quant <- getKallistoQuant(sampleFile = "~/FYP/fastq/GSE59114/GSE59114_sampFile.txt",
                                 index = "~/FYP/mouse.idx", resultDir = "~/FYP/kallisto_out/",
                                 outputPrefix = "GSE59114", nCores = 2)
) 

GSE48968.quant.time <- system.time(
  GSE48968.quant <- getKallistoQuant(sampleFile = "~/FYP/fastq/GSE48968/GSE48968_sampFile.txt",
                                     index = "~/FYP/mouse.idx", resultDir = "~/FYP/kallisto_out/",
                                     outputPrefix = "GSE48968", nCores = 2)
)
gse48968 <- readKallistoResults(GSE48968.quant)
saveRDS(gse48968, file = "~/FYP/GSE48968.rds")

cort.quant.time <- system.time(
  cort.quant <- getKallistoQuant(sampleFile = "~/FYP/fastq/cortex/cortex_sampFile.txt",
                                 index = "~/FYP/mouse.idx", resultDir = "~/FYP/kallisto_out/",
                                 outputPrefix = "cortex", nCores = 2,
                                 fragment_length = 200, fragment_standard_deviation = 20)
) 
 
embryo.quant.time <- system.time(
  embryo.quant <- getKallistoQuant(sampleFile = "~/FYP/fastq/embryo/embryo_sampFile.txt",
                                   index = "~/FYP/mouse.idx", resultDir = "~/FYP/kallisto_out/",
                                   outputPrefix = "embryo", nCores = 2,
                                   fragment_length = 200, fragment_standard_deviation = 20)
) 

## converting kallisto transcript level results to gene level abundances
## for demonstration, only shown for GSE70580
add.on <- function(data, name){
    for(i in seq_along(colnames(data))){
        colnames(data)[i] <- paste(colnames(data)[i], name, sep = "_")
    }
    return(data)
}

GSE70580 <- readRDS("~/matt_fyp/kallisto/GSE70580.rds")
tpm.GSE70580 <- as.data.frame(tpm(GSE70580))
write.table(tpm.GSE70580, "~/matt_fyp/GSE70580_TPM.txt", sep = "\t")


tpm.GSE70580 <- tpm.GSE70580[sort(rownames(tpm.GSE70580)),]
raw.ids <- rownames(tpm.GSE70580)

mart <- read.table("~/matt_fyp/kallisto/mart_export.txt",header = T,sep = "\t")

use <- mart[,c(2,4)]
rownames(use) <- use[,2]
use <- use[,-2, drop = FALSE]
use <- use[rownames(use) %in% raw.ids, ,drop = FALSE]
use <- use[sort(rownames(use)),, drop = F]
raw.ids <- sort(raw.ids)
GSE70580.names <- cbind(use, tpm.GSE70580)

agg.GSE70580 <- aggregate(GSE70580.names[,2:101], by = list(genename=GSE70580.names$Gene.name), FUN = sum)
agg.GSE70580 <- data.frame(agg.GSE70580[,-1], row.names = agg.GSE70580[,1])

write.table(agg.GSE70580, "~/matt_fyp/GSE70580_TPM.txt", sep = "\t")
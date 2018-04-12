## generating the boxplots
load("~/matt_fyp/Seurat_results/FYP_mouse_RSEM_2018-01-16.Robj")
metadata <- sObj@meta.data
pdf(file = "~/matt_fyp/comparison_results/nGene_boxplot.pdf")
ngene.boxplot <- ggplot(metadata, mapping = aes(x = orig.ident, y = nGene))
ngene.boxplot + geom_boxplot()
dev.off()

pdf(file = "~/matt_fyp/comparison_results/mito_boxplot.pdf")
mito.boxplot <- ggplot(metadata, mapping = aes(x = orig.ident, y = percent.mito))
mito.boxplot + geom_boxplot()
dev.off()


## generating RSEM-kallisto comparison dotplots
names <- colnames(total)[1:200]
names <- sub("_.+", "", names)
for (i in seq_along(names)) {
    compare <- total[, grep(names[i], colnames(total))]
    colnames(compare) <- c("kallisto", "RSEM")
    compare <- log(1 + compare)
    comp.plot <- ggplot(compare, aes(kallisto, RSEM)) +
        geom_point(alpha = 0.1) + geom_abline(slope = 1, color = "red") +
        xlim(0, ceiling(max(compare))) + ylim(0, ceiling(max(compare))) +
        ggtitle(names[i])
    ggsave(filename = paste0(names[i], ".png"),
           plot = comp.plot,
           device = "png")
}

## Same plot as above, but for ILC dataset
kall <- read.table("~/matt_fyp/final/GSE70580_TPM.txt", header = T)[,1, drop = F]
rsem <- read.table("~/matt_fyp/final/RSEM/RSEM.genes.results", header = T)
rsem$gene_id <- sub("^[^_]+_", "", rsem$gene_id)
rsem$gene_id <- sub("PAR_Y_", "", rsem$gene_id)
rsem <- rsem[,c("gene_id", "TPM")]
tpm <- aggregate(rsem[,"TPM", drop = F], by = list(x=rsem$gene_id), FUN = sum)
tpm <- data.frame(tpm[,"TPM"], row.names = tpm$x)
colnames(tpm) <- "RSEM"
colnames(kall) <- "Kallisto"


compare <- merge.data.frame(tpm, kall, by = "row.names")
compare <- data.frame(compare[,-1], row.names = compare$Row.names)
compare <- log(1 + compare)
comp.plot <- ggplot(compare, aes(Kallisto, RSEM)) +
    geom_point(alpha = 0.1) + geom_abline(slope = 1, color = "red") +
    xlim(0, ceiling(max(compare))) + ylim(0, ceiling(max(compare))) +
    ggtitle(paste0("Correlation = ", round(cor(x = compare$RSEM, y = compare$Kallisto), digits = 2)))
ggsave(filename = "RSEMvsKallisto.png",
       plot = comp.plot,
       device = "png")

## generating venn diagrams
kvg <- read.table("~/matt_fyp/comparison_results/kal_varGenes.txt", stringsAsFactors = F)[,1]
rvg <- read.table("~/matt_fyp/comparison_results/rsem_varGenes.txt", stringsAsFactors = F)[,1]
vG <- list(kallisto = kvg,
           RSEM = rvg)
dge <- list(kallisto = read.csv("~/matt_fyp/comparison_results/kal_DEG_2018-02-07.csv", stringsAsFactors = F)[,1],
            RSEM = read.csv("~/matt_fyp/comparison_results/rsem_DEG_2018-02-07.csv", stringsAsFactors = F)[,1])
require("VennDiagram")
vG.plot <- venn.diagram(vG , NULL, fill=c("darkmagenta", "darkblue"),
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4,
                          category.names=c("kallisto", "RSEM"), main="Variable Genes")
dge.plot <- venn.diagram(dge , NULL, fill=c("darkmagenta", "darkblue"),
                        alpha=c(0.5,0.5), cex = 2, cat.fontface=4,
                        category.names=c("kallisto", "RSEM"), main="Differentially Expressed Genes")
pdf(file = "~/matt_fyp/comparison_results/varGene_vennDiagram.pdf")
grid.draw(vG.plot)
dev.off()
pdf(file = "~/matt_fyp/comparison_results/DEG_vennDiagram.pdf")
grid.draw(dge.plot)
dev.off()



# Installing the necessary packages
if (!require(BiocManager)) install.packages("BiocManager")  
installifnot <- function (pkg){
 if (!require(pkg, character.only=T)){
 BiocManager::install(pkg)
 }
}
installifnot("pd.hg.u133a")
installifnot("hgu133a.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("GOstats")

# Reading the data
cel <- list.files("dades", pattern="*.CEL")
print(cel)

# Construction of the targets object
library(Biobase)
targetsDF <-read.csv(file=file.path(dataDir,"targets.csv"), header =TRUE, sep=";")
#DEFINE SOME VARIABLES FOR PLOTS
sampleNames <- as.character(targetsDF$ShortName)
sampleColor <- as.character(targetsDF$Col)
# Creation of the AnnotatedDataFrame object
targets <- AnnotatedDataFrame(targetsDF)

# Creation of the rawData object of the CEL files for subsequent pre-processing
CELfiles <- targetsDF$Filenames
rawData <- read.celfiles(file.path("C:/Users/KAREN/Documents/R/PEC2 Datos_ómicos/dades",CELfiles), phenoData=targets)

# Pre-processing: quality control and data exploration
# Visual exploration of data using ad-hoc methods:
#BOXPLOT
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", cex.axis=0.6, col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", cex=0.7, hang=-1)

# PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="",
scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
 pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
 loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
 xlab<-c(paste("PC1",loads[1],"%"))
 ylab<-c(paste("PC2",loads[2],"%"))
 if (is.null(colors)) colors=1
 plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts,
 xlim=c(min(pcX$x[,1])-100000,
max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000,
max(pcX$x[,2])+100000))
 text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
 title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep="
"), cex=0.8)
}
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data",
colors=sampleColor,
 formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

# Standardization
# The data is normalized and saved in a file called NormData.txt
eset<-rma(rawData)
Background correcting
Normalizing
Calculating Expression
write.exprs(eset, file.path("C:/Users/KAREN/Documents/R/PEC2 Datos_ómicos/results", "NormData.txt"))

# Filtering
library(genefilter)
annotation(eset) <- "hgu133a.db"
eset_filtered <- nsFilter(eset, var.func=IQR, require.entrez = TRUE,
 var.cutoff=0.75, var.filter=TRUE,
 filterByQuantile=TRUE)

# NUMBER OF GENES REMOVED
print(eset_filtered)

$filter.log
$filter.log$numDupsRemoved

$filter.log$numLowVar

$filter.log$numRemoved.ENTREZID

$filter.log$feature.exclude


#NUMBER OF GENES IN
print(eset_filtered$eset)

# Data post-processing
# Gene selection
# To carry out gene selection based on a linear model, two matrices must be constructed: one for design and one for contrast:

# design matrix
library(limma)
treat <- pData(filteredEset)$Group
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design) <- levels(lev)
rownames(design) <- sampleNames

# Contrast Matrix:
cont.matrix1 <- makeContrasts(untreated.vs.treated = untreated-treated, levels = design)
comparisonName <- "Efecto del tratamiento"

# Model estimates:
fit1 <- lmFit(filteredData, design) #estimaciones del modelo
fit.main1 <- contrasts.fit(fit1, cont.matrix1) #contraste
fit.main1 <- eBayes(fit.main1) #transformación de la varianza

# With the toptable function you can take the previous object and extract the number of rows you want that meet the criteria we ask for below:
topTab <- topTable (fit.main1, number=nrow(fit.main1),
coef="untreated.vs.treated", adjust="fdr",lfc=3, p.value=0.05)
View(topTab)

# Recording of results
# The genes are then annotated with the help of the microarray annotation (hgu133a.db)
library(hgu133a.db)
keytypes(hgu133a.db)

anotaciones<- AnnotationDbi::select (hgu133a.db,
keys=rownames(filteredData), columns=c("ENTREZID", "SYMBOL"))

# Combining annotations with topTab
library(dplyr)
topTabAnotada <- topTab %>%
 mutate(PROBEID=rownames(topTab)) %>%
 left_join(anotaciones) %>%
 arrange(P.Value) %>%
 select(7,8,9, 1:6)
head(topTabAnotada)

# The result can be written to a text file or an html file:
write.csv2(topTabAnotada, file= file.path(resultsDir,"Genes
seleccionados.csv"))
print(xtable(topTab,align="lllllll"),type="html",html.table.attributes=""
,
 file=file.path(resultsDir,"Genes seleccionados.html"))

# Visualization of the data using a Volcano plot or using a heatmap:
# volcano plot:
genenames <- AnnotationDbi::select(hgu133a.db,
 rownames(fit.main1), c("SYMBOL"))$SYMBOL
volcanoplot(fit.main1, highlight=10, names=genenames,
 main = paste("Differentially expressed genes",
colnames(cont.matrix1), sep="\n"))
abline(v = c(-3, 3))

pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = genenames,
 main = paste("Differentially expressed genes",
colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()

#heatmap plot:
selectedRows <- rownames(filteredData) %in% rownames(topTab)
selectedData <- filteredData[selectedRows,]

my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
library(gplots)
heatmap.2(selectedData,
 Rowv=TRUE,
 Colv=TRUE,
 main="HeatMap Induced.vs.WT FC>=3",
 scale="row",
 col=my_palette,
 sepcolor="white",
 sepwidth=c(0.05,0.05),
 cexRow=0.5,
 cexCol=0.9,
 key=TRUE,
 keysize=1.5,
 density.info="histogram",
 ColSideColors=c(rep("red",4),rep("blue",4)),
 tracecol=NULL,
 srtCol=30)

# Biological significance analysis
# Gene overrepresentation analysis. The percentage of all genes in the annotation is compared to the percentage of genes in the list being analyzed:
ibrary(hgu133a.db)
probesUniverse <- rownames(filteredData)
entrezUniverse<- AnnotationDbi::select(hgu133a.db, probesUniverse,
"ENTREZID")$ENTREZID
topProbes <- rownames(selectedData)
entrezTop<- AnnotationDbi::select(hgu133a.db, topProbes,
"ENTREZID")$ENTREZID

# Elimination of possible duplicates:
topGenes <- entrezTop[!duplicated(entrezTop)]
entrezUniverse <- entrezUniverse[!duplicated(entrezUniverse)]

# The GOstats package is used for enrichment analysis:
library(GOstats)
GOparams = new("GOHyperGParams",
 geneIds=topGenes, universeGeneIds=entrezUniverse,
 annotation="hgu133a.db", ontology="BP",
 pvalueCutoff=0.01)
GOhyper = hyperGTest(GOparams)






















setwd("~/Desktop/wgcna/")

library(WGCNA)

options(stringsAsFactors = F)

#
myData = read.table("breast_wgcna_exprs.txt", sep = "\t", dec = ".", header = T)
annot = read.table("notch.txt", sep = "\t", dec = ".", header = F)
rownames(myData) = myData$Hybridization.REF
loc = which(row.names(myData) %in% annot$V1)
sel = myData[loc,]
write.table(sel, "breast_wgcna_notch_targets.txt", sep = "\t", dec = ".", row.names = F, col.names = T)
#

myData = read.table("breast_wgcna_notch_targets.txt", sep = "\t", dec = ".", header = T)
rownames(myData) = myData$Hybridization.REF
myData = myData[1:506]
datExpr0 = as.data.frame(t(myData[, -1]))
names(datExpr0) = myData$Hybridization.REF
rownames(datExpr0) = names(myData)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

traitData = read.table("breast_wgcna_clinical.txt", sep = "\t", dec = ".", header = T)
lungs = rownames(datExpr0)
traitRows = match(lungs, traitData$id)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits$type, signed = F)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits[2]), main = "Sample dendrogram and trait heatmap")

# enableWGCNAThreads()

# STEP-BY-STEP #
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

k = softConnectivity(datE = datExpr0, power = 4)
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(k)
scaleFreePlot(k)

softPower = 4
adjacency = adjacency(datExpr0, power = softPower)

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEsmoduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

plotTOM = dissTOM^4
diag(plotTOM) = NA
TOMplot(plotTOM, geneTree, as.character(moduleColors), main = "TOM heatmap plot, module genes", terrainColors = T)

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits$type, use = "p")
moduleTraitPValue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPValue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits[2]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

type = as.data.frame(datTraits$type)
names(type) = "type"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr0, type, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(type), sep="")
names(GSPvalue) = paste("p.GS.", names(type), sep="")

sizeGrWindow(8, 9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise"
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),nrgcols=300,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

sizeGrWindow(8,7)
which.module="turquoise"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
turquoise = datExpr0[,moduleColors==which.module ]
write.table(turquoise, "turquoise.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

MMturquoise = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
# write.table(MMblue, "MMblue.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

eset = as.data.frame(t(turquoise))
nBRCA = read.table("nBRCA_rnaseq.txt", sep = "\t", dec = ".", header = T)
rownames(nBRCA) = nBRCA$gene_id
nBRCA$gene_id = NULL

loc = which(row.names(nBRCA) %in% rownames(eset))
sel = nBRCA[loc,]

conc_turq = cbind(eset, sel)

eset = as.matrix(log2(conc_turq+0.0000001))
hr = hclust(as.dist(1-cor(t(eset), method="spearman")), method="complete")
hc = hclust(as.dist(1-cor(eset, method="spearman")), method="complete")

library(gplots)

lighten <- function(color, factor = 9){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = greenred(50)
col = sapply(col, lighten)

heatmap.2(eset, Rowv=as.dendrogram(hr), Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", cexRow = 0.5, labCol = F,
          main = "Breast-assoc genes", xlab = "patients", ylab = "genes",
          ColSideColors = c(rep("red", 97), rep("blue", 58), rep("yellow", 228), rep("purple", 122), rep("green", 113)))
coord = locator(1)
legend(350,2, legend = c("basal-like", "HER2-enriched", "luminal A", "luminal B", "normal"), 
       col = c("red", "blue", "yellow", "purple", "green"), lty= 1, lwd = 8, xpd = T, cex = 0.7)

GS1 = as.numeric(cor(datTraits$type, datExpr0, use = "p"))
GeneSignificance = abs(GS1)
ModuleSignificance = tapply(GeneSignificance, moduleColors, mean, na.rm = T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors)

status = as.data.frame(datTraits$type)
names(status) = "status"
MET = orderMEs(cbind(MEs, status))
sizeGrWindow(5, 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8,
                      xLabelsAngle = 90)

hubs = chooseTopHubInEachModule(datExpr0, colorh = moduleColors, omitColors = "grey", power = 6)
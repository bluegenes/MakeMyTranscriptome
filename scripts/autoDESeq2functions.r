#! /usr/bin/Rscript
##############################################################################
##############################################################################
# Helper Functions for DESeq2 and EdgeR Analysis on RNASeq datasets 
#  Tessa Pierce
#  5.5.2015
#  helpul resources:
#DESeq2:
    #http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
    #http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ 
    #https://benchtobioinformatics.wordpress.com/category/deseq2/
    #https://gist.github.com/stephenturner/f60c1934405c127f09a6
#EdgeR:
    #http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
##############################################################################
##############################################################################


# set up R environment
source("http://bioconductor.org/biocLite.R") # bioconductor
library(BiocInstaller)

#install packages if necessary. BiocLite for some, regular R install packages for rest
installasNecessary <- function(packageName){
    if(!require(packageName, character.only=T, quietly=T)){
        if(packageName %in% c("DESeq2", "vsn", "VennDiagram","edgeR", "matrixStats", "fission")){
            biocLite(packageName)
        } else {
            install.packages(packageName)
        library(packageName,character.only=T) ## should write in something to catch errors here + print useful "could not install" statement
        }
    }
}
# if everything is failing, you may need to recompile ALL R packages:
#source("http://bioconductor.org/biocLite.R")
#pkgs <- rownames(installed.packages())
#biocLite(pkgs, type="source")

installasNecessary("DESeq2")
plotMA.deseq2 <- plotMA # specify DESeq2 MA function so doesn't get masked
installasNecessary("vsn")
installasNecessary("VennDiagram")
installasNecessary("edgeR")
#installasNecessary("matrixStats")
installasNecessary("RColorBrewer")
installasNecessary("calibrate")
installasNecessary("gplots")
installasNecessary("ggplot2")
installasNecessary("genefilter")
installasNecessary("fission")

##########################################
           ## DESeq2 Functions ##
##########################################
runAllDESeq2 <- function(baseDir,treatmentInfo, deseqBasename, ddsData,PVAL=.05, filter=FALSE, filterBM=5){#, model="~condition"){
    #make new output directory for the basenames
    dirName = paste("deseq2",deseqBasename,sep="_")
    dir.create(file.path(baseDir,dirName), showWarnings = FALSE)
    setwd(file.path(baseDir, dirName))
    # Capture sessionInfo from R session 
    writeLines(capture.output(sessionInfo()), paste(dirName, "_sessionInfo.txt", sep="")) 
    deseqBaseDir <- getwd()
    # DESeq2 transformations (regularized log transform + variance stabilizing transform) for clustering & plotting
    mycols <- brewer.pal(8, "Dark2")[1:length(colnames(ddsData))] # get colors for conditions --> if have >8, need to switch color palette!!!
    #transformData(ddsData, deseqBasename, mycols) #this function transforms, plots, and outputs the transformed count data 
    resultsList <- getAllContrastResults(ddsData,treatmentInfo) #make all pairwise contrasts; store each results df as element of resultsList
    for (name in names(resultsList)){
        newDeseqBasename = paste(deseqBasename, '_', name,sep="")
        dir.create(file.path(deseqBaseDir,newDeseqBasename), showWarnings = FALSE) #create subfolder for this contrast results
        setwd(file.path(deseqBaseDir,newDeseqBasename))
        resultsD = resultsList[[name]]
        if (filter == TRUE){ #option to filter by base mean 
            newDeseqBasename = paste(newDeseqBasename, '_filterbm', filterBM, sep="")
            resultsD <- resultsD[resultsD$baseMean>filterBM, ]} #elim any gene with basemean < filterBM counts
        plotdeseq2indepFiltering(resultsD, newDeseqBasename)
        resData <- mergeResultsWithCountData(resultsD, ddsData, newDeseqBasename,PVAL) # merge DESeq2 results with count data + write out in tsv format 
        names(resData)[1] <- "Gene" #add column header for 'Gene' column
        plotDispersionsAndVolcano(ddsData, resData, newDeseqBasename, PVAL)
        #heatmap NOT WORKING #
        #plotDEHeatmap(ddsData,resData, newDeseqBasename, 50, nrow(treatmentInfo), mycols) # default top 50 genes for heatmap --> will this work??? Maybe need to plot differently
        deseqsigData <- subset(resData, padj < PVAL) #get signif genes
        write.table(as.data.frame(deseqsigData),file= paste(newDeseqBasename, "_p_", PVAL, ".tsv", sep=''), sep ='\t', quote=F, row.names=FALSE)
        deseqGenes <- deseqsigData$Gene #get gene names (mostly useful for edgeR comparsions --> maybe don't need here..)
        }
}




getAllContrastResults <- function(ddsData,treatmentInfo){
    resultsList = list()
    i = 1
    while (i <= ncol(treatmentInfo)){
        factor = colnames(treatmentInfo)[i]
        levels = levels(treatmentInfo[,i])
        pairwiseComb = combn(levels, 2) # get character vectors of each pair.
        j = 1
        while (j <= ncol(pairwiseComb)){ # number of comparisons
            resultsList[[paste(rev(pairwiseComb[,j]),collapse="_")]] <- runDESeqContrast(ddsData, factor, pairwiseComb[,j][2],pairwiseComb[,j][1]) # note order=reverse so earlier levels = denominators 
            j = j+1
            }
        i = i+1
        }
    resultsList 
    }

#function to run a single contrast & return the ordered dictionary
runDESeqContrast <- function(ddsData, factorName, numerator, denominator){ 
    resultsDict <- results(ddsData, contrast = c(factorName, numerator, denominator))
    resultsD <- resultsDict[order(resultsDict$padj),] 
    }
# just tired of writing these commands out. need some simple functions for evaluation work.
orderResults <- function(resultsDict){resultsDict[order(resultsDict$padj),]}
getSigResults <- function(resultsDict, PVAL){ subset(resultsDict, padj < PVAL)}

# Merge pvalue info with normalized count data
mergeResultsWithCountData <- function(resultsD, ddsData, basename,PVAL){
    resultsData <- merge(as.data.frame(resultsD), as.data.frame(counts(ddsData, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resultsData)[1] <- "Gene"
    write.table(resultsData, file= paste(basename, "_deseq2.tsv", sep=''), sep ='\t', quote=F, row.names=FALSE)
    write.table(table(resultsD$padj<PVAL ),file=paste(basename, "_padj.", PVAL, ".txt", sep=""), quote=F )
    return(resultsData)
}

plotDispersionsAndVolcano <- function(ddsData, resData, basename, PVAL){
    png(paste(basename, '_dispersonEstimates_deseq2.png', sep=''))
    plotDispEsts(ddsData, main="Dispersion plot")
    dev.off()
    png(paste(basename, '_tMA_deseq2.png', sep=''), 1500, 1000, pointsize=20) #turner MA plot
    maplot(resData, PVAL, FALSE, main="MA Plot")
    dev.off()
    png(paste(basename, '_Labeled_tMA_deseq2.png', sep=''), 1500, 1000, pointsize=20) #turner MA plot
    maplot(resData, PVAL, TRUE, main="MA Plot")
    dev.off()
    png(paste(basename, '_volcano_deseq2.png', sep=''), 1200, 1000, pointsize=20)
    volcanoplot(resData, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
    dev.off()
    png(paste(basename, '_pvalHist_deseq2.png', sep=''))
    hist(resData$pvalue, breaks=50,col="grey")
    dev.off()
}

plotdeseq2indepFiltering <- function(resultsD, basename){
    png(paste(basename, '_filteringHist_deseq2.png', sep=''))
    plot(attr(resultsD,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
    dev.off()
}

## need to improve heatmap function...
plotDEHeatmap <- function(ddsData,resdata, basename, numTopGenes, numSamples, COLORS){ #input = a merged dds 'results' dataset
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) #get some colors
    colsToSelect = 7 + numSamples # get index of last countData column
    selectDE <- as.matrix(resdata[2:numTopGenes, 8:colsToSelect]) # select top NUM DE genes + counts in all samples 
    row.names(selectDE) <- resdata[2:numTopGenes,1]
    png(paste(basename, '_heatmapDE_deseq2.png', sep=''))
    heatmap(selectDE, margin=c(10,6))
    dev.off()
    png(paste(basename, '_fancierHeatmapDE_deseq2.png', sep=''))
    heatmap.2(selectDE, col = hmcol, scale="none", dendrogram="both", trace="none", margin=c(10,6)) #need to change this command slightly
    dev.off()
    ## better heatmap --> plot largest difference from baseMean
    #r <- resdata[(resdata, padj<.05),]$Gene #get only signif genes
    #sigDDS <- ddsData[row.names(ddsData) %in% r,]
    #topVarGenes <- head(order(rowVars(assay(sigDDS)), decreasing=TRUE), numTopGenes)
    #png(paste(basename, '_topVarHeatmap_deseq2.png', sep=''), width=800, height=800)
    #heatmap.2(assay(sigDDS)[topVarGenes,], scale="row", trace="none", dendrogram="both", col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),ColSideColors=COLORS[colData(sigDDS)$condition],cexRow=1,margins=c(20,16))
    #dev.off()
}


transformData <- function(ddsData, basename, COLORS){
    #mycols <- brewer.pal(8, "Dark2")[1:length(colnames(ddsData))] # get plenty of colors 
    rldData <- rlogTransformation(ddsData, blind=FALSE)
    vsdData <- varianceStabilizingTransformation(ddsData,blind=FALSE)
    plotTransformedData(paste(basename, "_RLD", sep=""),rldData,ddsData, COLORS)
    plotTransformedData(paste(basename, "_VSD", sep=""),vsdData,ddsData, COLORS)
    #plotEffectofTransformation(rldData, vsdData, ddsData, basename)
}

plotTransformedData <- function(basename, transformedData, ddsData, COLORS, numGenes=50){
    #plot sample distances (heatmap + PCA)
    sampleDist <- as.matrix(dist(t(assay(transformedData))))  # sample distance heatmap (turner)
    #plot sample distance heatmap
    png(paste(basename, "_sampleDistHmap_deseq2.png",sep=""), w=1000, h=1000, pointsize=20) 
    heatmap.2(as.matrix(sampleDist), key=F, trace="none", col=colorpanel(100, "black", "white"),ColSideColors=COLORS[ddsData$condition],RowSideColors=COLORS[ddsData$condition],margin=c(15, 15), main=paste("Sample Distance Matrix: ", basename, sep=""))
    dev.off()
    #plot default Deseq2 PCA
    png(paste(basename,"_defaultPCA_deseq2.png",sep="")) #default PCA
    plotPCA(transformedData, intgroup = c("condition"))
    dev.off()
    #plot Turner PCA
    png(paste(basename,"_tPCA_deseq2.png",sep=""), 1000, 1000, pointsize=20) #turner pca
    rld_pca(transformedData, colors=COLORS, intgroup="condition", xlim=c(-75, 75), ylim=c(-40, 40))
    dev.off()
    #plot histogram of transformed data
    png(paste(basename,"hist.png", sep=""))
    hist(assay(transformedData)) #do a histogram, just to look at the data
    dev.off()
# from http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2#working-with-rlog-transformed-data
#The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s  average across all samples. Hence, we center and scale each genes’ values across samples, and plot a heatmap.
    topVarGenes <- head(order(rowVars(assay(transformedData)), decreasing=TRUE), numGenes)
    png(paste(basename, '_topVarHeatmap_deseq2.png', sep=''), width=800, height=800)
    heatmap.2(assay(transformedData)[topVarGenes,], scale="row", trace="none", dendrogram="both", col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),ColSideColors=COLORS[colData(transformedData)$condition],cexRow=1,margins=c(20,16))
    dev.off()
    #write table of transformed counts (in case want to use them later)
    write.table(as.data.frame(assay(transformedData)),file=paste(basename, "transformed_counts_deseq2.tsv",sep=""), quote=F, sep='\t') #write transformed counts
}

plotEffectofTransformation <- function(rld, vsd, ddsData, basename){
   #plot mean vs std dev for the original + both transformations (DESeq2)
    png(paste(basename, '_transformationEffect_deseq2.png', sep=''))
    par(mfrow=c(1,3))
    notAllZero <- (rowSums(counts(ddsData)) > 0)
    meanSdPlot(log2(counts(ddsData,normalized=TRUE)[notAllZero,] + 1))
    meanSdPlot(assay(rld[notAllZero,]))
    meanSdPlot(assay(vsd[notAllZero,]))
    dev.off()
}

##########################################
           ## EdgeR Functions ##
##########################################
runEdgeR <- function(edgeR_countData, edgeRbasename){
    dataset <- calcNormFactors(edgeR_countData)
    dataset <- estimateCommonDisp(dataset)
    dataset <- estimateTagwiseDisp(dataset)
    return(dataset)
}

# DESeq does independent filtering automatically --> edgeR recommends doing some. Here, we filter all below 2 cpm in at least two samples, which is the smallest      number of samples in group (as recommended in the edgeR user guide)
filterAndRunEdgeR <- function(edgeData, minNumReplicates){
    keep <- rowSums(cpm(edgeData)>1) >= minNumReplicates
    filteredData <- edgeData[keep,]
    edgeRFiltered <- runEdgeR(filteredData)
    return(edgeRFiltered)
}
getSigEdgeRGenes <- function(edgeRBasename, edgeR_ET, PVAL){
    numGood <- sum(p.adjust(edgeR_ET$table$PValue,method="BH") < PVAL)
    if( numGood == 0){ numGood =1}
    sigData = topTags(edgeR_ET, n=numGood)
    write.table(sigData,file=paste(edgeRBasename, "_", PVAL, '_edgeR.tsv', sep=""), sep ='\t', quote=F, row.names=TRUE) # write out edgeR sig genes 
    return(sigData)
}

plotEdgeR <- function(edgeRData, et, edgeRBasename, signifData){
    png(paste(edgeRBasename,"_BCV_edgeR.png", sep=""))    
    plotBCV(edgeRData)
    dev.off()
    plotMDSandVariance(edgeRData, edgeRBasename)
    plotEdgeRSmear(edgeRData, et, edgeRBasename, signifData)
}

plotMDSandVariance <- function(edgeData, eBasename){ #from edgeR tutorial: http://cgrlucb.wikispaces.com/file/view/edgeR_Tutorial.pdf
    png(paste(eBasename,"_meanVariance_edgeR.png", sep=""))    
    plotMeanVar(edgeData, show.raw.vars=TRUE, show.tagwise.vars=TRUE,show.binned.common.disp.vars=FALSE,show.ave.raw.vars=FALSE, NBline = TRUE, nbins = 100, pch = 16,xlab ="Mean Expression (Log10 Scale)", ylab = "Variance (Log10 Scale)", main = "Mean-Variance Plot" )
    dev.off()
    png(paste(eBasename,"_MDS_edgeR.png", sep=""))    
    plotMDS(edgeData, main = "MDS Plot for Count Data", labels = colnames(edgeData$counts))   
    dev.off()
}

plotEdgeRSmear <- function(edgeData, et, basename, signif){
    png(paste(basename,"_smear_edgeR.png", sep=""))   
    plotSmear(et, de.tags=rownames(signif))
    abline(h=c(-1, 1), col="blue")
    dev.off()
}


##########################################
# Plotting Functions from Stephen Turner  https://gist.github.com/stephenturner/f60c1934405c127f09a6
########################################## 
# Turner RLD PCA
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

## Turner MA Plot
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
## Turner Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
######## End Turner plotting functions #######


###### unused functions ######
plotDefaultMA <- function(ddsData, basename){
    png(paste(basename, '_MA_deseq2.png', sep=''))
    plotMA.deseq2(ddsData,ylim=c(-2,2),main="DESeq2 MA") # using built-in ma function
    dev.off()
}

#working, but inefficient; poor implementation. Rewritten above.
GETAllContrastResults <- function(ddsData, treatmentInfo){
    resultsList = list()
    i =1
    while (i <= ncol(treatmentInfo)){
        factor = colnames(treatmentInfo)[i]
        levels = levels(treatmentInfo[, i])
        j = length(levels)
        while (j > 1){
            resultsList[[paste(levels[j], '_', levels[j-1], sep = "")]] <- runDESeqContrast(ddsData, factor, levels[j],levels[j-1])
            if (j-2 >= 1){
                resultsList[[paste(levels[j], '_', levels[j-2], sep = "")]] <- runDESeqContrast(ddsData, factor, levels[j],levels[j-2])
            }
            if (j-3 >= 1){ 
                resultsList[[paste(levels[j], '_', levels[j-3], sep = "")]] <- runDESeqContrast(ddsData, factor, levels[j],levels[j-3])
            }
            j = j-1
        }
        i = i+1
    }
    resultsList
}

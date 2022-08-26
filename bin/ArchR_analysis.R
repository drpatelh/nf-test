#!/usr/bin/env Rscript
packages <- c("AnnotationHub","ArchR","argparse", "BiocManager","BSgenome","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Hsapiens.NCBI.GRCh38", "Cairo", "devtools","dplyr","GenomicFeatures","ggplot2","glue","magick","parallel","Seurat", "yaml")
suppressMessages(invisible(lapply(packages, library, character.only = TRUE, quietly=TRUE)))

# UTIL FUNCTIONS
# Splits a comma delimited string into a vector 
# if commas are present. Otherwise returns vector
# with string as single element 
splitIfList <- function(str) {
    if (grepl(",", str, fixed=TRUE)) {
        splitString <- strsplit(str, ",")
        return(unlist(splitString[[1]]))
    } else {
        return(c(str))
    }
}

formatAsScientificName <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

## Validation Functions 
validateFragmentsAndSampleNames <- function(fragments, sampleNames) {
    if (length(fragments) != length(sampleNames)) {
        stop(glue("The number of fragment file paths and the number of samples names provided must be equal. Currently {length(fragments)} and {length(sampleNames)}"))
    }
    for (path in fragments) {
        if (!file.exists(path)) {
            stop(glue("One of the specified fragment paths does not exist: {path}"))
        }
    } 
}

validateUsedAnnotationSet <- function(genomeAlias, gtf, BSGenome) {
    if (is.null(genomeAlias) && is.null(gtf)) {
        stop("An ArchR Genome Alias (-a) or a GTF (-g) & BSGenome (-b) must be specified")
    }
    if (!is.null(genomeAlias) && !is.null(gtf)) {
        message("Both an ArchR Genome Alias and GTF were provided. Using GTF to create a custom annotation set")
    } 
    else if (!is.null(genomeAlias)) {
        message(glue("Using built-in ArchR Genome {genomeAlias} for annotations"))
    }
    else if (!is.null(gtf) && is.null(BSGenome)) {
        stop("A BSGenome must be provided when creating a custom annotations from a GTF")
    } 
    else if (!is.null(gtf) && !is.null(BSGenome)) {
        message("Using passed GTF and BSGenome to create custom annotations")
    }
}

getGenome <- function(BSGenomeName) {
    if (!require(BSGenomeName, character.only=TRUE)) {
            print(glue("{BSGenomeName} not found. Installing now"))
            tryCatch(
                {
                    BiocManager::install(BSGenomeName)
                    library(packageName, character.only=TRUE)
                },
                error=function(cond){
                    message(glue("{BSGenomeName} could not be retrieved, ensure that BSGenome Package specified in genomes.json exists"))
                    message("Original Error message")
                    stop(cond)
                }
            )
    }
    genome <- getBSgenome(BSGenomeName)
    seqlevelsStyle(genome) <- "UCSC"
    return(genome)
}

# Gets OrgDb object from AnnotationHub for specified organism 
getOrgDb <- function(speciesName) {
    hub <- AnnotationHub()
    allAvailableSpecies <- hub$species
    if (speciesName %in% allAvailableSpecies) {
        orgDbSearch <- query(hub, c(speciesName,"OrgDb"))
        if (length(orgDbSearch) > 0) {
            print(glue("Utilizing {orgDbSearch$title[1]}"))
            orgdb <- orgDbSearch[[1]]
            return(orgdb)
        } else {
            stop(glue("No OrgDb could be found for {speciesName}"))
        }
    } else {
        stop(glue("There is no information on AnnotationHub for {speciesName}, Available species: {unique(sort(allAvailableSpecies))}"))
    }
}

# Makes TxDb object from gtf 
makeTxDb <- function(speciesName,gtf,referenceConvention) {
    txdb <- makeTxDbFromGFF(gtf, format="gtf", "Personal", speciesName)
    if (referenceConvention == "ENSEMBL") {
        seqlevelsStyle(txdb) <- "UCSC"
    }
    return(txdb)
}

# Add blacklist (means either more package installations or 
# interaction with AnnotationHub (download at runtime))
makeGenomeAnnotation <- function(genome) {
    chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
    genomeAnnotation <- createGenomeAnnotation(
        genome=genome,
        chromSizes=chromSizes
    )
}

# Creates custom set of gene annotations 
makeGeneAnnotation <- function(orgName, gtfPath, referenceConvention) {
    txdb <- makeTxDb(orgName, gtfPath, referenceConvention)
    orgdb <- getOrgDb(orgName)
    anStyle <- if (referenceConvention == "ENSEMBL") "ENSEMBL" else "SYMBOL"
    geneAnnot <- createGeneAnnotation(
        genome=NULL,
        TxDb=txdb,
        OrgDb=orgdb,
        # Column in orgDb that matches TxDb gene_id
        annoStyle=anStyle
    )
    return(geneAnnot)
}

# Gets the minimum number of non-zero features in the tile 
# Matrices of the passed ArrowFiles (used as an estimate of minimum variable features) 
getMinFeatureCounts <- function(arrowFiles) {
    featCounts <- c()
    for (file in arrowFiles) {
        tilemat <- getMatrixFromArrow(
        ArrowFile = file,
        useMatrix = "TileMatrix",
        useSeqnames = NULL,
        cellNames = NULL,
        ArchRProj = NULL,
        verbose = TRUE,
        binarize = TRUE
        )
        vec <- tabulate(assays(tilemat)[[1]]@i + 1)
        nzFeats = sum(vec != 0)
        featCounts <- c(featCounts, nzFeats)
    }
    return(min(featCounts))
}

# Creates a grid of plots for each marker in @markers, 
# coloring the given @project's UMAP embedding the geneScore 
plotMarkersOnEmbedding <- function(markers, project) {
    sampleCount <- length(getSampleNames(project))
    embed <- if (sampleCount==1) "UMAP" else "UMAPHarmony" 
    p <- plotEmbedding(
        ArchRProj = project, 
        colorBy = "GeneScoreMatrix", 
        name = markers, 
        embedding = embed,
        imputeWeights = getImputeWeights(project)
    )
    if (length(markers) == 1) {
        p <- list(p)
    }
    p2 <- lapply(p, function(x){
    x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 4.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
            )
        }
    )
    markersOnEmbedding <- do.call(cowplot::plot_grid, c(list(ncol = 2),p2))
    return(markersOnEmbedding)
}


# Writes TSSEnrichment and DoubletScores to metadata.tsv table 
# within analysis output directory 
writeProjectMetadata <- function(archRProj, outputDir) {
    Barcode <- archRProj$cellNames
    # sampleCount <- length(unique(archRProj$Sample))
    # if (sampleCount == 1) {
    #    Barcode <- lapply(Barcode, function(x) {
    #      return(sub(".*#", "", x))
    #    })
    # }
    TSSEnrichment <- archRProj$TSSEnrichment
    DoubletScore <- archRProj$DoubletScore
    class.df<- data.frame(Barcode, TSSEnrichment, DoubletScore)
    write.table(class.df, sep="\t", file=glue('{outputDir}/metadata.tsv'), quote=FALSE, row.names=FALSE)
}


# Writes gene score matrix, features.tsv and barcodes.tsv 
writeGeneScoreMatrix <- function(archrproj, outputDir) {
    dir.create(file.path(outputDir, "geneScore"), showWarnings = FALSE)
    # Features
    features <- getFeatures(
        ArchRProj = archrproj,
        useMatrix = "GeneScoreMatrix",
        select = NULL,
        ignoreCase = TRUE
    )
    write.table(features, file=glue('{outputDir}/geneScore/features.tsv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
    # Barcodes
    write.table(archrproj$cellNames, file=glue('{outputDir}/geneScore/barcodes.tsv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
    # Matrix 
    sparseMatrices <- assays(getMatrixFromProject(archrproj, "GeneScoreMatrix"))
    writeMM(sparseMatrices[[1]], file=glue('{outputDir}/geneScore/matrix.mtx'))
}


plotMarkerGenes <- function(top10Df, allMarkers, archrproj, chosenMarkerPath) {
    if (!is.null(chosenMarkerPath)) {
        listDF <- read.csv(chosenMarkerPath, header=FALSE)
        markersOfInterest <- unlist(listDF[[1]])
        chosenMarkerPlot <- plotMarkersOnEmbedding(markersOfInterest, archrproj)
        suppressWarnings(plotPDF(chosenMarkerPlot, name=glue("markerPlots.pdf"), width=9, height=9, ArchRProj=archrproj))
    }

    if (nrow(top10Df) > 0) {
    for (clusterVal in unique(top10Df$cluster)) {
        clusterMarkers <- subset(top10Df, subset=top10Df$cluster==clusterVal)
        if (nrow(clusterMarkers) > 0) {
            clusterMarkerPlot <- plotMarkersOnEmbedding(clusterMarkers$name, archrproj)
            suppressWarnings(plotPDF(clusterMarkerPlot, name=glue("markerPlot_{clusterVal}.pdf"), width=9, height=9, ArchRProj=archrproj))
        }
    }
        plotLog2FCBool = ncol(allMarkers)<= 2
        # Heatmap creation 
        heatmap <- plotMarkerHeatmap(
        seMarker = allMarkers, 
        cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
        transpose = TRUE,
        plotLog2FC = plotLog2FCBool
        )
    heatmapImage <- ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")    
    suppressWarnings(plotPDF(heatmapImage, name="markersHeatmap.pdf", width=9,height=9,ArchRProj=archrproj))
    }
}

# Writes the top 10 markers per cluster to a tsv for later use 
writeTop10MarkersPerCluster <- function(markers, outputDir) {
    markerList <- getMarkers(markers, cutOff="FDR <= 0.01 & Log2FC >= 1.5")
    for (cluster in sort(names(markerList))) {
        clusterMarkers <- markerList[[cluster]]
        topN <- min(c(nrow(clusterMarkers), 10))
        if (cluster == 'C1') {
            df <- setNames(data.frame(matrix(ncol = (length(colnames(clusterMarkers)) + 1), nrow = 0)), c(colnames(clusterMarkers),'cluster'))
        } 
        if (topN > 0) {
            addition <- head(clusterMarkers, topN)
            addition$cluster <- cluster
            df <- rbind(df, as.data.frame(addition))
        }
    }
    if (nrow(df) > 0) {
        write.table(df, file=glue("{outputDir}/Top10Markers.tsv"), sep="\t", quote=FALSE)
    } 
    return(df) 
}


# Returns markers as a summarized project object. 
calcMarkerGenes <- function(proj) {
    markers <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
    )
    return(markers)
}


# Plot UMAP projections for the passed project. 
plotClustering <- function(project, nameAppend, sampleCount, additionalFields) {
    toPlot <- list()
    p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters", embedding = "UMAP",  baseSize=14)
    toPlot <- append(toPlot, list(p1))
    p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP", baseSize=14)
    toPlot <- append(toPlot, list(p2))
    for (field in additionalFields) {
         plot <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = field, embedding = "UMAP", baseSize=14)
         toPlot <- append(toPlot, list(plot))
    }
    suppressWarnings(plotPDF(plotList=toPlot, name=glue("UMAP-{nameAppend}.pdf"),width=8,height=8,ArchRProj=project))
    if (sampleCount > 1) {
        toPlotHarmony <- list()
        p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony",  baseSize=14)
        toPlotHarmony <- append(toPlotHarmony, list(p3))
        p4 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony",  baseSize=14)
        toPlotHarmony <- append(toPlotHarmony, list(p4))
        for (field in additionalFields) {
        plot <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = field, embedding = "UMAPHarmony", baseSize=14)
        toPlotHarmony <- append(toPlotHarmony, list(plot))
        }
        suppressWarnings(plotPDF(toPlotHarmony, name=glue("UMAPHarmony-{nameAppend}.pdf"),width=8,height=8,ArchRProj=project))
    }
    project <- addImputeWeights(project)
    return(project)
}


# Performs dimensionality reduction 
clusterAndDimRed <- function(archrProj, sampleNames, varFeatureVal) {
    proj <- addIterativeLSI(ArchRProj = archrProj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations=quantParams$LSI_Iterations, varFeatures=varFeatureVal)
    proj <- addClusters(input = proj, method="Seurat", reducedDims = "IterativeLSI", maxClusters=quantParams$maxClusters)
    proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",  nNeighbors=quantParams$nNeighbors)

    # Batch effect correction if multiple samples 
    if (length(sampleNames) > 1) {
        proj <- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "Sample")
        proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors=quantParams$nNeighbors)
    }
    return(proj)
}

## DATA VISUALIZATION FUNCS 
# Creates a violin plot of any metadata field 
# in cellColData (specified using @field)
makeFieldPlot <- function(archrproj, field) {
    p1 <- plotGroups(
    ArchRProj = archrproj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = field,
    plotAs = "violin",
    alpha = 0.4, 
    addBoxPlot = TRUE,
    baseSize=9
    )
    suppressWarnings(plotPDF(
        p1,
        name=glue("{field}.pdf"),
        width=8,
        height=8,
        ArchRProj=archrproj
    ))
}


# Plots fragments size and TSS Enrichment in a single PDF for the given 
# @archrproj 
makeQCPlots <- function(archrproj) {
    p2 <-plotFragmentSizes(ArchRProj = archrproj)
    p3 <- plotTSSEnrichment(ArchRProj = archrproj)
    suppressWarnings(plotPDF(
        p2,
        p3,
        name="QC.pdf",
        width=8,
        height=8,
        ArchRProj=archrproj
    ))
}

# Builds ArchR project from the provided arguments. Loads one in if it already exists at @outputDir
# minTSS and minFrags temporarily set to non-zero values until fragment filtering logic added to pipeline 
# inputted fragments.tsv.gz pre-filtered by qc step    
buildArchrProject <- function(fragmentFiles, genomeAnnotation, geneAnnotation, outputDir, sampleNames, passingBarcodeVector) {
    if (file.exists(glue("{outputDir}/Save-ArchR-Project.rds"))) {
        proj <- loadArchRProject(outputDir)
    } else {
    arrowFiles <- createArrowFiles(
        inputFiles = fragmentFiles,
        sampleNames = sampleNames,
        validBarcodes = passingBarcodeVector,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        geneAnnotation = geneAnnotation,
        genomeAnnotation = genomeAnnotation,
        minFrags=quantParams$minimumFragments,
        minTSS=quantParams$minimumTSS
    )
    proj <- ArchRProject(
        ArrowFiles=arrowFiles, 
        outputDirectory=outputDir, 
        copyArrows=TRUE,
        showLogo=FALSE,
        geneAnnotation=geneAnnotation, 
        genomeAnnotation=genomeAnnotation
    )}
    return(proj)  
}

# Use ArchR's built in sets of annotations for your organism 
# instead of passing your own gtf 
useDefaultGenomeAnnots <- function(genomeAlias,referenceConvention) {
    # In genomes metadata file denote whether or not ArchR has built in support of genome
    addArchRGenome(genomeAlias, install=TRUE)
    includeChrPrefix <- !referenceConvention == "ENSEMBL"
    addArchRChrPrefix(chrPrefix = includeChrPrefix)
    return(list("geneAnnotation"=getGeneAnnotation(), "genomeAnnotation"=getGenomeAnnotation()))
}

#Creates a custom set of genome and gene  
createCustomGenomeAnnots <- function(orgName, referenceConvention, BSGenomeName, gtfPath) {
    genome <- getGenome(BSGenomeName)
    genomeAnnotation <- makeGenomeAnnotation(genome)
    geneAnnotation <- makeGeneAnnotation(orgName, gtfPath, referenceConvention)
    return(list("geneAnnotation"=geneAnnotation, "genomeAnnotation"=genomeAnnotation))
}

# Returns the list of 
extractPassingCells <- function(passingMetricsPath) {
    if (!is.null(passingMetricsPath) && file.exists(passingMetricsPath)) {
        allCellData <- read.table(file = passingMetricsPath, sep = '\t', header = TRUE)
        passingCellData <- allCellData[ which(allCellData$Filter=='Pass'), ]
        return(passingCellData)
    } else {
        return(NULL)
    }
}

# Assumption: Can currently only process one sample at a time, sampleNames 
# only contains one sample
addExternalMetadata <- function(archrProj, metadataDf, sampleNames) {
    if (!is.null(metadataDf)) {
        wantedCells <- unique(archrProj$cellNames)
        sample <- sampleNames[1]
        metadataDf$FullBarcode <- sprintf(glue("{sample}#%s"), metadataDf$Barcode)
        metadataForWantedCells <-  dplyr::filter(metadataDf, metadataDf$FullBarcode %in% wantedCells)
        metadataForWantedCells$FracMito <- metadataForWantedCells$MitoReads / metadataForWantedCells$TotalReads
        archrProj <- addCellColData(archrProj, data=as.numeric(metadataForWantedCells$TotalReads), name="TotalReads", cells=metadataForWantedCells$FullBarcode, force=TRUE)
        archrProj <- addCellColData(archrProj, data=as.numeric(metadataForWantedCells$UniqueReads), name="UniqueReads", cells=metadataForWantedCells$FullBarcode, force=TRUE)
        archrProj <- addCellColData(archrProj, data=as.numeric(metadataForWantedCells$FracMito), name="FracMito", cells=metadataForWantedCells$FullBarcode, force=TRUE)
        archrProj <- addCellColData(archrProj, data=as.numeric(metadataForWantedCells$FRiP), name="FRiP", cells=metadataForWantedCells$FullBarcode, force=TRUE)
        archrProj <- addCellColData(archrProj, data=as.numeric(metadataForWantedCells$FRiT), name="FRiT", cells=metadataForWantedCells$FullBarcode, force=TRUE)
    }
    return(archrProj)   
}

calcAndFilterDoublets <- function(archrProj, variableFeatures) {
    archrProj <- addDoubletScores(
        input=archrProj, 
        k=10,
        knnMethod="UMAP", 
        LSIMethod=1,
        LSIParams=list(varFeatures=variableFeatures)
    )
    return(filterDoublets(archrProj)) 
}

removeUnclusteredCells <- function(archrProj) {
    hasUnclusteredCells <- NA %in% names(table(archrProj$Clusters,useNA='ifany'))
    if (hasUnclusteredCells) {
        totalCells <- nCells(archrProj)
        clusteredCells <- dim(archrProj@reducedDims$IterativeLSI$matSVD)[1]
        unclusteredCells <- totalCells - clusteredCells
        warning(glue("WARNING: {unclusteredCells} out of {totalCells} Cells were not clustered by iterativeLSI.\n This suggests the chosen QC-parameters were too lenient. NA cells are being filtered out for now (consider rerunning with stricter parameters)"))
        archrProj <- archrProj[!is.na(archrProj$Clusters), ]
    }
    return(archrProj)
}


# Runs basic ArchR analysis, generating clustering and  
# QC images for metrics calcualted by ArchR
runArchRAnalysis <- function(genomeAlias, orgName, fragments, sampleNames, markerGenesPath, gtfPath, referenceConvention, BSGenomeName,outputDir, passingMetrics) {
    addArchRVerbose(verbose = FALSE)
    if (!file.exists("ArchR")) {
        dir.create("ArchR")
    }
    annotations <- if ((is.null(gtfPath) & !is.null(genomeAlias))) useDefaultGenomeAnnots(genomeAlias, referenceConvention) else createCustomGenomeAnnots(orgName,referenceConvention,BSGenomeName,gtfPath)
    passingCells <- extractPassingCells(passingMetrics)
    passingBarcodes <- if(!is.null(passingCells)) passingCells$Barcode else NULL
    proj <- buildArchrProject(fragments, annotations[['genomeAnnotation']], annotations[['geneAnnotation']], outputDir, sampleNames, passingBarcodes)
    proj <- addExternalMetadata(proj, passingCells, sampleNames)
    arrowFiles <- getArrowFiles(proj)
    variableFeatureCount <- min(c(getMinFeatureCounts(arrowFiles),25000))

    # QC Visualization
    makeQCPlots(proj)
    makeFieldPlot(proj, "TSSEnrichment")

    saveArchRProject(proj, load=FALSE)

    # Atleast 1000 variable features required for iterative LSI to run 
    if (variableFeatureCount >= 1000) {
        proj <- calcAndFilterDoublets(proj, variableFeatureCount)
        proj <- clusterAndDimRed(proj, sampleNames, variableFeatureCount)
        saveArchRProject(proj, load=FALSE)

        proj <- removeUnclusteredCells(proj)

        # Data Visualization 
        if (quantParams$plotUMAP) {
            proj <- plotClustering(proj, "pipeline-run", length(sampleNames), c("TSSEnrichment","FRiP", "FRiT", "FracMito", "UniqueReads"))
        }
        if (quantParams$getMarkerGenes) {
            markers <- calcMarkerGenes(proj)
            top10Df <- writeTop10MarkersPerCluster(markers, outputDir)
            plotMarkerGenes(top10Df, markers, proj, markerGenesPath)
        }
        writeGeneScoreMatrix(proj, outputDir)
        writeProjectMetadata(proj, outputDir)
    } else {
        stop("ONE OR MORE OF THE SPECIFIED SAMPLE HAS <1000 VARIABLE REGIONS DETECTED BETWEEN CELLS (1000 NEEDED FOR DIMENSIONALITY REDUCTION).\n  THIS MAY BE THE RESULT OF EITHER POOR DATA QUALITY OR POOR QC FILTERING.\nTRY THE FOLLOWING:\n- Increase the value for 'minimumTSS' in the user specified qcAndArchRParams file (see 'Quality Control and ArchR Analysis Parameters' section of the README for guidance)")
    }
    saveArchRProject(proj,load=FALSE)
}

parser <- ArgumentParser()
# Required 
parser$add_argument("-s", "--species",type="character",required=TRUE, help="Full scientific species name (i.e-Homo sapiens)")
parser$add_argument('-f', "--fragments", type="character", required=TRUE, help="Comma delimited list of paths to fragments.tsv.gz files to included in analysis (i.e- path1,path2,path3)")
parser$add_argument('-sN', "--sampleNames", type="character", required=TRUE, help="The names of the samples corresponding to the elements in the fragments argument")
parser$add_argument("-o", "--outputDir", type="character", required=TRUE, help="Name of output directory to be created")
parser$add_argument("-r", "--referenceConvention", choices=c("ENSEMBL", "UCSC"),type="character", default=NULL, required=TRUE, help="Contig naming format convention followed by gtf and reference genome used for alignment")
parser$add_argument("-q", "--quantitativeParams", type="character", required=TRUE, help="Path to .yml containing parameters for various analysis steps")
# Optional (Depending on annotation derivation)
parser$add_argument("-b", "--BSGenome", required=FALSE, default=NULL, help="Name of the known BSGenome to use in order to derive a custom genome annotation")
parser$add_argument("-a", "--archrGenomeAlias", choices=c("hg19","hg38","mm9", "mm10"), type="character", default=NULL, required=FALSE,help="If using ArchR's built in annotations, specify which set to use")
parser$add_argument("-g", "--gtf", type="character", default=NULL, required=FALSE, help=" (optional) Path to GTF containing gene annotations of interest")
parser$add_argument("-qc", "--qcMetrics", type="character", default=NULL, required=FALSE, help="Path to tsv denoting which cells passed qc filtering")
# Optional 
parser$add_argument("-th", "--threads", type="integer", required=FALSE, default=3, help="Number of threads to use to run the analysis")
parser$add_argument("-m", "--markerGenes",type="character", required=FALSE, default=NULL,help=" (optional) A line delimited list of Marker genes to create plots for")
parser$add_argument("-se", "--seed", type="integer", default=NULL, required=FALSE, help=" (optional) Seed to use in order to maintain reproducibility")

args <- parser$parse_args()

fragments <- splitIfList(args$fragments)
sampleNames <- splitIfList(args$sampleNames)
validateFragmentsAndSampleNames(fragments, sampleNames)
validateUsedAnnotationSet(args$archrGenomeAlias, args$gtf, args$BSGenome)
outputDir <- glue("ArchR/{args$outputDir}")

if (!is.null(args$seed)) {
    set.seed(args$seed)
}
if (args$threads  > 1) {
    # Done in order to allow for random seeding when 
    # multiple threads are utilized 
    RNGkind("L'Ecuyer-CMRG")
}

addArchRThreads(threads = args$threads)
quantParams <- yaml.load_file(args$quantitativeParams)
runArchRAnalysis(args$archrGenomeAlias, formatAsScientificName(args$species), fragments, sampleNames, args$markerGenes, args$gtf, args$referenceConvention, args$BSGenome, outputDir, args$qcMetrics)

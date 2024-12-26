#load required packages
required.packages <- c("tidyverse",
                       "data.table",
                       "edgeR",
                       "limma",
                       "DESeq2",
                       "ashr",
                       "combinat",
                       "BiocParallel")
installed.packages <- rownames(installed.packages())
if (all(required.packages %in% installed.packages)) {
  for (package in required.packages) {
    suppressMessages(library(package, character.only=TRUE))
  }
} else {
  required.packages[!required.packages %in% installed.packages]
  stop(paste("Install these packages: ", paste(required.packages[!required.packages %in% installed.packages], collapse=",")))
}


#generate variables based on command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Expected arguments: 1) dge object 2) dea method 3) model equation")
}
dge_file <- args[1]
method <- args[2]
if (grepl("shrink", method, ignore.case=TRUE)) { 
  shrink=TRUE
} else {
  shrink=FALSE
}
model <- args[3]
model <- ifelse(startsWith(model, "~"), model, paste0("~", model))
model_notilde <- gsub("~", "", model)
lfc_cutoff <- 1
padj_cutoff <- 0.05


#create output files
outdir <- file.path(dirname(dge_file), "differentialExpression")
dir.create(outdir, showWarnings=F, recursive=T)


#read in DGE object
dgel <- readRDS(dge_file)
dge <- dgel[["All"]]


#remove any reference samples
keep <- is.na(dge$samples$Reference)
dge$samples <- dge$samples[keep,]
dge$counts <- dge$counts[,keep]
dge$length <- dge$length[,keep]
dge$log2nx <- dge$log2nx[,keep]



runDESeq2 <- function(dge, design, contrast, shrink=FALSE) {
  # Replace counts with high Cook's distance 
  # only necessary if design includes continuous predictors b/c won't be done automatically
  # requires one column named "group"
  if (length(setdiff(unique(as.numeric(design)), c(0,1))) > 0) {
    dds <- DESeqDataSetFromMatrix(round(dge$counts), colData=dge$samples, design=as.formula("~group"))
    dds <- DESeq(dds, minReplicatesForReplace=3)
    counts <- assays(dds)[["replaceCounts"]]
  } else {
    counts <- round(dge$counts)
  }
  dds <- DESeqDataSetFromMatrix(counts, colData=dge$samples, design=design)
  assays(dds)[[2]] <- dge$length
  names(assays(dds))[2] <- "avgTxLength"
  
  # Pre-filter genes with no length info
  keep <- !is.na(rowSums(assays(dds)[[2]]))
  dds <- dds[keep,]
  dds <- estimateSizeFactors(dds, type="poscount")
  
  # Pre-filter genes with low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Run deseq, includes:
  # 1) calculation of size factors (if not already done)
  # 2) estimation of dispersion
  # 3) negative binomial GLM
  dds <- DESeq(dds, minReplicatesForReplace=3, parallel=TRUE, BPPARAM=MulticoreParam(4))
  #plotDispEsts(dds)
  #meanSdPlot(assay(normTransform(dds)))
  
  # Find differentially expressed genes
  if (shrink) { 
    res <- lfcShrink(dds, contrast=contrast, type="ashr") 
  } else { 
    res <- results(dds, contrast=contrast) 
  }
  genes <- dge$genes[match(rownames(res), dge$genes$gene_id_simplified),,drop=FALSE]
  res <- data.frame(genes, res)
  res$lfc <- res$log2FoldChange
  res$padj <- res$padj
  
  return(res)
}



#create design matrix
design <- model.matrix(as.formula(model), data=dge$samples)
colnames(design)[grepl("Intercept", colnames(design))] <- "Intercept"
write.table(design, file=file.path(outdir, paste0("DEA_design_matrix.txt")),  sep="\t", quote=F, row.names=T)


#specify contrasts of interest in design matrix
tmp <- paste0("group", unique(dge$samples$group)) # must use group variable for DEA
groups <- colnames(design)[which(colnames(design) %in% tmp)]
if (length(groups) >= 2) {
  contrasts <- as.data.frame(combn(groups, 2)) # create all pairwise combos
  contrasts <- apply(contrasts, 2, function(x) paste(x, collapse="-"))
  names(contrasts) <- gsub("group", "", gsub("-", "vs", contrasts)) # name all pairwise combos
} else {
  contrasts <- groups
  names(contrasts) <- gsub("group", "", gsub("-", "vs", contrasts))
}
design.contrasts <- makeContrasts(contrasts=contrasts, levels=design)
colnames(design.contrasts) <- names(contrasts)
write.table(as.data.frame(design.contrasts), file=file.path(outdir, paste0("DEA_contrast_matrix.txt")),  sep="\t", quote=F, row.names=T)


#generate standard output files
summary <- list()
for (i in 1:ncol(design.contrasts)) {
  contrast=design.contrasts[,i,drop=FALSE]
  name=colnames(design.contrasts)[i]
  if (grepl("DESeq2", method)) { res <- runDESeq2(dge, design, contrast, shrink=shrink) }
  res$contrast <- name

  #write out dea results
  write.table(res, file=file.path(outdir, paste0(name, "_differential_expression.txt")),  sep="\t", quote=F, row.names=F)
  summary[[i]] <- res
}


#write out dea summary statistics
summary <- do.call(rbind, summary)
summary <- summary %>%
  group_by(contrast) %>%
  summarize(upGenes=sum(lfc>=lfc_cutoff & padj<padj_cutoff, na.rm=T),
            downGenes=sum(lfc<=-lfc_cutoff & padj<padj_cutoff, na.rm=T))
write.table(summary, file=file.path(outdir, "DEA_summary.txt"),  sep="\t", quote=F, row.names=F)




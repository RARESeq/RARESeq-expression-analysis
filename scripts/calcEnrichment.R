#load required packages
required.packages <- c("tidyverse",
                       "data.table",
                       "MASS")
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
if (length(args) < 4) {
  stop("Expected arguments: 1) dge file 2) enrichment method 3) gene set file 4) out dir name")
}
dge_file <- args[1]
method_q <- args[2]
method <- strsplit(method_q, "_")[[1]][1]
q <- strsplit(method_q, "_")[[1]][2]
q <- as.numeric(ifelse(is.na(q), 0.9, q))
gene_set_file <- args[3]
dir_name <- args[4]
np <- 1000
a <- 0.05
set.seed(100)


#create output files
outdir <- file.path(dirname(dge_file), "enrichment", dir_name)
dir.create(outdir, showWarnings=F, recursive=T)


#read in zscores, and re-format
dgel <- readRDS(dge_file)


#read in gene file, and re-format
gene_set <- fread(gene_set_file, header=T, data.table=F)
if (!"gene_id_simplified" %in% colnames(gene_set)) {
  stop("Warning! Gene signatures file must contain gene_id_simplified column header.")
}

# adjust weight column
if ("weight" %in% colnames(gene_set)) {
} else if ("Weight" %in% colnames(gene_set)) {
  gene_set$weight <- gene_set$Weight
} else {
  gene_set$weight <- 1
  print("No weights found in signature file, assigning all gene weights to 1.")
}
#adjust gene set name column
if (!"gene_set_name" %in% colnames(gene_set)) {
  gene_set$gene_set_name <- dir_name # give gene set same name as dir
}



# calculates enrichment score per gene
calcEnrichment <- function(zscores, genes, weights, method="stouffer", q=0.9, np=1000) {
  set.seed(1234)
  
  # various methods for calculating enrichment
  if (method=="sum") {
    f <- function(z, w, q) { apply(z, 2, function(x) sum(x) / length(x)) }
  } else if (method=="stouffer") {
    f <- function(z, w, q) { apply(z, 2, function(x) sum(x*w) / sqrt(sum(w**2))) }
  } else if (method=="percentile") {
    f <- function(z, w, q) { apply(z, 2, function(x) quantile(x, q, na.rm=TRUE)) }
  }
  
  up_gns <- genes[weights>0]
  up_wts <- weights[weights>0]
  up_idx <- match(up_gns, rownames(zscores))
  
  dwn_gns <- rev(genes[weights<=0])
  dwn_wts <- -1*(rev(weights[weights<=0]))
  dwn_idx <- match(dwn_gns, rownames(zscores))
  
  if (sum(is.na(c(up_idx, dwn_idx))) > 0) {
    stop("Warning! Gene set ids do not match expression object")
  } else if (length(up_idx)==0) {
    stop("Warning! Gene set does not contain any genes with positive weights")
  } else {
    if (length(dwn_idx) > 0) {
      up_zscores.in <- zscores[up_idx,] # genes will be ordered from highest to lowest weight
      dwn_zscores.in <- zscores[dwn_idx,] # genes will be ordered from highest to lowest weight
      zscores.out <- zscores[-c(up_idx, dwn_idx),]
      zscores.out <- zscores.out[sample(1:nrow(zscores.out), replace=F),] # randomly order genes
      
      #calculate enrichment score
      test <- f(up_zscores.in, up_wts, q) + f(dwn_zscores.in, dwn_wts, q)
      
      #calculate null distribution
      null <- as.data.frame(matrix(nrow=np, ncol=ncol(zscores)))
      colnames(null) <- colnames(zscores)
      for (i in 1:nrow(null)) {
        up_scrambled <- sample_n(zscores.out, length(up_gns), replace=TRUE)
        dwn_scrambled <- sample_n(zscores.out, length(dwn_gns), replace=TRUE)
        null[i,] <- f(up_scrambled, up_wts, q) + f(dwn_scrambled, dwn_wts, q)
      }
    } else {
      up_zscores.in <- zscores[up_idx,] # genes will be ordered from highest to lowest weight
      zscores.out <- zscores[-up_idx,]
      zscores.out <- zscores.out[sample(1:nrow(zscores.out), replace=F),] # randomly order genes
      
      #calculate enrichment score
      test <- f(up_zscores.in, up_wts, q)
      
      #calculate null distribution
      null <- as.data.frame(matrix(nrow=np, ncol=ncol(zscores)))
      colnames(null) <- colnames(zscores)
      for (i in 1:nrow(null)) {
        up_scrambled <- sample_n(zscores.out, length(up_gns), replace=TRUE)
        null[i,] <- f(up_scrambled, up_wts, q)
      }
    }
  } 
  
  #calculate empirical pvalue based on null
  pvalues <- data.frame(ID=colnames(zscores), zscore=test, pvalue=NA)
  for (j in 1:ncol(null)) {
    p <- sum(null[,j] >= test[j])
    pvalues$pvalue[j] <- ifelse(p==0, 1/np, p/np)
  }
  
  #summarize as enrichment score
  pvalues$enrichment_score <- pvalues$zscore
  pvalues$enrichment_score[pvalues$pvalue > a] <- 0
  
  return(pvalues)
}







# calculate pvalues for each unique gene_set
combined <- list()
sets <- unique(gene_set$gene_set_name)
for (l in 1:length(dgel)) {
  dge <- dgel[[l]]
  run <- names(dgel)[l]
  zscores <- dge$z
  
  for (s in 1:length(sets)) {
    
    name <- sets[s]
    set <- gene_set[which(gene_set$gene_set_name==name & gene_set$gene_id_simplified %in% rownames(zscores)),]
    
    if (nrow(set)>0) {
      
      set <- set[order(set$weight, decreasing=T),]
      genes <- set$gene_id_simplified
      weights <- set$weight
      out <- calcEnrichment(zscores, genes, weights, method=method, q=q, np=np)
      combined[[length(combined)+1]] <- data.frame(gene_set_name=name, run=run, out)
      
    } 
  }
}
combined <- do.call(rbind, combined)



#generate standard output files
write.table(combined, file.path(outdir, "enrichment_scores.txt"), sep="\t", quote=F, row.names=F, col.names=T)




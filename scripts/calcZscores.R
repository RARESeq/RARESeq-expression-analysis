#load required packages
required.packages <- c("tidyverse",
                       "data.table")
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
if (length(args) < 1) {
  stop("Expected arguments: 1) dge object")
}
dge_file <- args[1]
set.seed(100)


#create output files
outdir <- dirname(dge_file)


#read in DGE object
dgel <- readRDS(dge_file)



setReference <- function(ref_log2nx) {
  ref <- data.frame(avg = apply(ref_log2nx, 1, mean), 
                    std = apply(ref_log2nx, 1, sd),
                    var = apply(ref_log2nx, 1, var))
  rownames(ref) <- rownames(ref_log2nx)

  # calculate avg and sd of lowest 5% of genes to avoid any genes have SD=0
  q <- quantile(ref$avg, 0.05, na.rm=T)
  if (q > min(ref$avg, na.rm=T)) {
    u <- which(ref$avg <= q)
  } else {
    q <- quantile(ref$avg, 0.2, na.rm=T)
    u <- which(ref$avg <= q)
  }
  ref$avg[u] <- mean(c(as.matrix(ref_log2nx[u,])), na.rm=T)
  ref$std[u] <- sd(c(as.matrix(ref_log2nx[u,])), na.rm=T)
  
  return(ref)
}





for (l in 1:length(dgel)) {
  dge <- dgel[[l]]
  
  #identify reference samples
  ref_i <- which(!is.na(dge$samples$Reference))
  if (length(ref_i)==0) { 
    print("Warning! No reference samples provided, so using all samples provided.")
    ref_i <- seq(1:length(dge$samples$ID))
  } 
  
  #calculate avg and std dev of reference expression
  ref <- setReference(dge$log2nx[,ref_i])
  
  #calculate zscores & output
  if (length(ref_i) < ncol(dge$log2nx)) {
    log2nx <- dge$log2nx[,-ref_i]
  } else {
    log2nx <- dge$log2nx
  }
  zscores <- as.data.frame((log2nx-ref$avg) / ref$std)
  rownames(zscores) <- rownames(log2nx)
  dge$z <- zscores
  dgel[[l]] <- dge
  }


saveRDS(dgel, gsub(".rds", "_withZscores.rds", dge_file))



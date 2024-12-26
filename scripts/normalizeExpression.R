#load required packages
required.packages <- c("tidyverse",
                       "data.table",
                       "ggpubr",
                       "edgeR",
                       "DESeq2",
                       "RUVSeq",
                       "mixtools")
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
#note that norm and correct methods can be "NONE
args=commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Expected arguments: 1) tximport object 2) gene info 3) sample info 5) normalization method 6) correction method 7) negative control genes for RUVg correction method (optional)")
}
txi_file <- args[1]
gene_file <- args[2]
sample_file <- args[3]
norm <- args[4]
correct <- args[5]
method <- paste(norm, correct, sep="_") 
if (length(args)>5) { neg_ctrl_file <- args[6] } 
rmGroup <- "Control"
set.seed(100)



print("normalizing expression...")


#create output directory
outdir <- file.path(dirname(txi_file), method)
dir.create(outdir, showWarnings=F, recursive=T)


#read in R data object
txi <- readRDS(file=txi_file)


#read in gene info
#this file should include only genes of interest (i.e. coding captured genes)
if (exists("gene_file")) {
  gene_info <- fread(gene_file, data.table=F, header=T)
  if (!"gene_id_simplified" %in% colnames(gene_info)) {
    stop("Warning! Gene info file must contain gene_id_simplified column header.")
  }
  
  if ("V1V2shared" %in% colnames(gene_info)) {
    shared_gene_info <- gene_info[which(!is.na(gene_info$V1V2shared)),]
  } 
  
} else {
  stop("Warning! Gene info file not found.")
}


#read in sample info
#this file should include only samples of interest
if (exists("sample_file")) {
  sample_info <- fread(sample_file, data.table=F, header=T)
  if (!"ID" %in% colnames(sample_info)) {
    stop("Warning! Sample info file must contain ID column header.")
  }
} else {
  stop("Warning! Gene info file not found.")
}
if (!"group" %in% colnames(sample_info)) { sample_info$group <- 1 }
if (!"DonorID" %in% colnames(sample_info)) { sample_info$DonorID <- unlist(lapply(strsplit(sample_info$ID, "[._-]"), function(x) x[1])) }


#read in negative control gene info, if provided
if (exists("neg_ctrl_file")) {
  neg_ctrl_info <- fread(neg_ctrl_file, data.table=F, header=T)
  if (!"gene_id_simplified" %in% colnames(gene_info)) {
    stop("Warning! Negative control info file must contain gene_id_simplified column header.")
  } else {
    neg_ctrl_ids <- neg_ctrl_info$gene_id_simplified
  }
} else {
  neg_ctrl_ids <- NULL
}







#normalize expression data
calcNX <- function(txi, gene_info, sample_info, norm) {
  set.seed(100)
  
  #extract counts, length, and TPM data
  genes <- intersect(rownames(txi$counts), gene_info$gene_id_simplified)
  samples <- intersect(colnames(txi$counts), sample_info$ID)
  if (length(samples) == 0) {
    stop("Warning! Provided sample IDs do not match")
  } else if (length(genes) == 0) {
    stop("Warning! Provided gene IDs do not match")
  } 
  cts <- txi$counts[genes, samples]
  len <- txi$length[genes, samples]
  tpm <- txi$abundance[genes, samples]
  sample_info <- sample_info[match(samples, sample_info$ID),,drop=FALSE]
  gene_info <- gene_info[match(genes, gene_info$gene_id_simplified),,drop=FALSE]
  ref_i <- which(!is.na(sample_info$Reference))

  if (norm == "TPM") {
    nx.in <- t(t(tpm)/colSums(tpm)) * 1e6 #rescales after filtering for genes of interest
    if (sum(is.na(nx.in))>0) { stop("Warning! Found NAs in input matrix")}
    lib.sizes <- colSums(nx.in)
    norm.factors <- calcNormFactors(nx.in, method="none")
    
    #nx <- nx.in
    log2nx <- log(nx.in + 1, 2)
    
  } else {
    if (norm == "TMM") {
      nx.in <- cts
      if (length(ref_i)==0) { 
        ref.in <- NULL 
      } else {
        ref.in <- as.numeric(apply(nx.in[,ref_i,drop=FALSE], 1, median, na.rm=TRUE))
      }
      if (sum(is.na(nx.in))>0) { stop("Warning! Found NAs in input matrix")}
      lib.sizes <- colSums(nx.in)
      norm.factors <- calcNormFactorsCustom(nx.in, method="TMMwsp", refColumn=ref.in)
    } else if (norm == "NONE") {
      nx.in <- cts
      if (sum(is.na(nx.in))>0) { stop("Warning! Found NAs in input matrix")}
      lib.sizes <- rep(1e6, ncol(nx.in))
      norm.factors <- rep(1, ncol(nx.in))
    } else {
      stop("Warning! Normalization method is invalid. Must be one of: TPM, TMM, NONE.")
    }
    #nx <- cpm(nx.in, lib.size=lib.sizes*norm.factors, log=FALSE)
    log2nx <- cpm(nx.in, lib.size=lib.sizes*norm.factors, log=TRUE)
  }

  #combine all info in new DGE obj
  dge <- DGEList(genes=as.data.frame(gene_info),
                 samples=as.data.frame(sample_info),
                 counts=as.matrix(cts))
  dge$samples$lib.size <- lib.sizes
  dge$samples$norm.factors <- norm.factors
  dge$samples$effective.lib.size <- lib.sizes*norm.factors
  dge$log2nx <- log2nx
  dge$length <- len
  return(dge)
}






#customize edgeR calcNormFactors function to use specified reference sample(s)
calcNormFactorsCustom <- function(object, ...)
  UseMethod("calcNormFactors")

calcNormFactors.DGEList <- function(object, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, 
                                    logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
  #	Scale normalization of RNA-Seq data, for DGEList objects
  #	Created 2 October 2014.  Last modified 2 June 2020.
{
  if(!is.null(object$offset)) warning("object contains offsets, which take precedence over library\nsizes and norm factors (and which will not be recomputed).")
  object$samples$norm.factors <- calcNormFactors(object=object$counts, lib.size=object$samples$lib.size, method=method, refColumn=refColumn, 
                                                 logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p)
  object
}

calcNormFactors.default <- function(object, lib.size=NULL, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, 
                                    logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
  #	Scale normalization of RNA-Seq data, for count matrices
  #	Mark Robinson, Gordon Smyth and edgeR team
  #	Created 22 October 2009. Last modified 2 June 2020.
{
  #	Check object
  x <- as.matrix(object)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)
  
  #	Check lib.size
  if(is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if(length(lib.size) != nsamples) {
      if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
      lib.size <- rep_len(lib.size,nsamples)
    }
  }
  
  #	Check method
  #	Backward compatability with previous name
  if(length(method)==1L && method=="TMMwzp") {
    method <- "TMMwsp"
    message("TMMwzp has been renamed to TMMwsp")
  }
  method <- match.arg(method)
  
  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]
  
  #	Degenerate cases
  if(nrow(x)==0 || nsamples==1) method="none"
  
  #	Calculate factors
  # select ref column using standard TMM strategy, i.e. which sample norm factor matches closest with the 75th percentile
  if (is.null(refColumn)) {
    #f50 <- suppressWarnings(.calcFactorQuantile(data=x, lib.size=lib.size, p=0.50))
    #f75 <- suppressWarnings(.calcFactorQuantile(data=x, lib.size=lib.size, p=0.75))
    #if (median(f75) < 1e-20) {
      refColumn <- which.max(colSums(sqrt(x)))
    #} else {
    #  refColumn <- which.min(abs(f75 - mean(f75)))
    #}
    print(paste0("Using reference sample: ", colnames(x)[refColumn]))
    ref <- x[,refColumn]
    libsize.ref <- lib.size[refColumn]
  } else {
    ref <- refColumn[!allzero]
    libsize.ref <- sum(refColumn)
  }
  f <- rep_len(NA_real_,nsamples)
  for(i in 1:nsamples) {
    f[i] <- .calcFactorTMMwsp(obs=x[,i], ref=ref, libsize.obs=lib.size[i], libsize.ref=libsize.ref, 
                              logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, name=colnames(x)[i])
  }

  #	Factors should multiple to one
  f <- f/exp(mean(log(f)))

  #	Output
  names(f) <- colnames(x)
  f
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	Created 16 Aug 2010. Last modified 12 Sep 2020.
{
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
  if(min(f)==0) warning("One or more quantiles are zero")
  f / lib.size
}

.calcFactorTMMwsp <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, name="")
  #	TMM with pairing of singleton positive counts between the obs and ref libraries
  #	Gordon Smyth
  #	Created 19 Sep 2018. Last modified 9 Jun 2020.
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  #	epsilon serves as floating-point zero
  eps <- 1e-14

  #	Identify zero counts
  pos.obs <- (obs > eps)
  pos.ref <- (ref > eps)
  npos <- 2L * pos.obs + pos.ref

  #	Remove double zeros and NAs
  i <- which(npos==0L | is.na(npos))
  if(length(i)) {
    obs <- obs[-i]
    ref <- ref[-i]
    npos <- npos[-i]
  }

  #	Check library sizes
  if(is.null(libsize.obs)) libsize.obs <- sum(obs)
  if(is.null(libsize.ref)) libsize.ref <- sum(ref)

  #	Pair up as many singleton positives as possible
  #	The unpaired singleton positives are discarded so that no zeros remain
  zero.obs <- (npos == 1L)
  zero.ref <- (npos == 2L)
  k <- (zero.obs | zero.ref)
  n.eligible.singles <- min( sum(zero.obs), sum(zero.ref))
  if(n.eligible.singles > 0L) {
    refk <- sort(ref[k],decreasing=TRUE)[1:n.eligible.singles]
    obsk <- sort(obs[k],decreasing=TRUE)[1:n.eligible.singles]
    obs <- c(obs[!k],obsk)
    ref <- c(ref[!k],refk)
  } else {
    obs <- obs[!k]
    ref <- ref[!k]
  }

  #	Any left?
  n <- length(obs)
  if(n==0L) return(1)

  #	Compute M and A values
  obs.p <- obs / libsize.obs
  ref.p <- ref / libsize.ref
  M <- log2( obs.p / ref.p )
  A <- 0.5 * log2( obs.p * ref.p )

  #	If M all zero, return 1
  if(max(abs(M)) < 1e-6) return(1)

  #	M order, breaking ties by shrunk M
  obs.p.shrunk <- (obs+0.5) / (libsize.obs+0.5)
  ref.p.shrunk <- (ref+0.5) / (libsize.ref+0.5)
  M.shrunk <- log2( obs.p.shrunk / ref.p.shrunk )
  o.M <- order(M, M.shrunk)

  #	A order
  o.A <- order(A)

  #	Trim
  loM <- as.integer(n * logratioTrim) + 1L
  hiM <- n + 1L - loM
  keep.M <- rep_len(FALSE,n)
  keep.M[o.M[loM:hiM]] <- TRUE
  loA <- as.integer(n * sumTrim) + 1L
  hiA <- n + 1L - loA
  keep.A <- rep_len(FALSE,n)
  keep.A[o.A[loA:hiA]] <- TRUE
  keep <- keep.M & keep.A
  M <- M[keep]

  #	Average the M values
  if(doWeighting) {
    obs.p <- obs.p[keep]
    ref.p <- ref.p[keep]
    v <- (1-obs.p)/obs.p/libsize.obs + (1-ref.p)/ref.p/libsize.ref
    w <- (1+1e-6) / (v+1e-6)
    TMM <- sum(w*M) / sum(w)
  } else {
    TMM <- mean(M)
  }

  2^TMM
}








#correct gene expression based on platelets
runRUV <- function(dge, correct, neg_ctrl_ids, k=1) {
  set.seed(100)
  
  #RUVg method used to normalize based on set of negative control genes
  #Requires pre-defined gene list(s)
  #K specifies how many factors of unwanted expression to estimate (1 recommended)
  if (correct=="RUVg") {
    
    if (length(neg_ctrl_ids) > 0) {
      
      neg_ctrl_ids <- intersect(neg_ctrl_ids, dge$genes$gene_id_simplified)
      #print(paste0(length(neg_ctrl_ids), " negative control genes used for RUVg normalization."))
      
      ruv.out <- RUVg(dge$log2nx, neg_ctrl_ids, k=k, isLog=TRUE, center=TRUE)
      w <- as.data.frame(ruv.out$W)
      dge$samples <- cbind(dge$samples, w)
      dge$log2nx <- ruv.out$normalizedCounts
      
    } else {
      stop("Warning! Negative control genes must be provided for RUVg correction.")
    }
    
  #RUVs method used to normalize based on set of negative control genes AND set of negative control samples
  #Requires pre-defined gene list(s)
  #Requires pre-defined reference samples
  #K specifies how many factors of unwanted expression to estimate (1 recommended)
  } else if (correct=="RUVs") {
    
    if (length(neg_ctrl_ids)==0) { neg_ctrl_ids <- dge$genes$gene_id_simplified }
    #print(paste0(length(neg_ctrl_ids), " negative control genes used for RUVs normalization."))
    
    if ("Reference" %in% colnames(dge$samples)) { 
      ref_ids <- t(data.matrix(which(!is.na(dge$samples$Reference)))) 
      #print(paste0(sum(ref_ids > 0), " negative control samples used for RUVs normalization."))
    } else {
      stop("Warning! Negative control samples must be provided for RUVs correction.")
    }
    
    ruv.out <- RUVs(dge$log2nx, neg_ctrl_ids, k=k, ref_ids, isLog=TRUE)
    w <- as.data.frame(ruv.out$W)
    dge$samples <- cbind(dge$samples, w)
    dge$log2nx <- ruv.out$normalizedCounts
    
  } else if (correct=="NONE") {
    # no correction done
  } else {
    stop("Warning! Correction method is invalid. Must be one of: NONE, RUVg, RUVs.")
  }

  return(dge)
}







#remove biological replicates
removeReplicates <- function(dge, rm="none") {
  
  if (!"DonorID" %in% colnames(dge$samples)) {
    stop("Warning! Sample info file must contain DonorID column header in order to remove replicates.")
  }
  
  #calculate sample correlations within group
  withinGroupR <- list()
  for (i in 1:nrow(dge$samples)) {
    id <- dge$samples$ID[i]
    grp <- dge$samples$group[i]
    grp.ids <- which(dge$samples$ID!=id & dge$samples$group==grp)
    R <- cor(dge$log2nx[,id], apply(dge$log2nx[,grp.ids,drop=FALSE], 1, mean))
    withinGroupR[[i]] <- data.frame(ID=id, withinGroupR=R)
  }
  withinGroupR <- do.call(rbind, withinGroupR) %>%
    inner_join(dge$samples) %>%
    group_by(DonorID, group) %>%
    mutate(bestReplicate=withinGroupR==max(withinGroupR))
  dge$samples$bestReplicate <- withinGroupR$bestReplicate[match(dge$samples$ID, withinGroupR$ID)]
  
  if (!(rm=="none"|rm=="n")) {
    if (rm=="all"|rm=="a") {
      rm_ids <- which(!dge$samples$bestReplicate)
    } else {
      rm_ids <- which(dge$samples$group %in% rm & !dge$samples$bestReplicate)
    }
    if (length(rm_ids) > 0) {
      dge$samples <- dge$samples[-rm_ids,]
      dge$counts <- dge$counts[,-rm_ids]
      dge$log2nx <- dge$log2nx[,-rm_ids]
      dge$length <- dge$length[,-rm_ids]
    }
  }
  return(dge) 
}











#if V1 and V2 panels used
#need to normalize separately
if ("Selector" %in% colnames(sample_info)) {
  
  dge <- calcNX(txi, shared_gene_info, sample_info, norm)
  dge <- runRUV(dge, correct, neg_ctrl_ids) 
  dge <- removeReplicates(dge, rm=rmGroup)
  dgel <- list(dge)
  names(dgel)[1] <- "All"
  
  if (any(grepl("V1", sample_info$Selector))) {
    
    v1_dge <- calcNX(txi, gene_info, sample_info[which(sample_info$Selector=="V1"),], norm)
    v1_dge <- runRUV(v1_dge, correct, neg_ctrl_ids) 
    v1_dge <- removeReplicates(v1_dge, rm=rmGroup)
    dgel[[2]] <- v1_dge
    names(dgel)[2] <- "V1"
  } 
  
  if (any(grepl("V2", sample_info$Selector))) {
  
    v2_dge <- calcNX(txi, shared_gene_info, sample_info[which(sample_info$Selector=="V2"),], norm)
    v2_dge <- runRUV(v2_dge, correct, neg_ctrl_ids) 
    v2_dge <- removeReplicates(v2_dge, rm=rmGroup)
    dgel[[3]] <- v2_dge
    names(dgel)[3] <- "V2"
  }
  
  
} else {
  
  dge <- calcNX(txi, gene_info, sample_info, norm)
  dge <- runRUV(dge, correct, neg_ctrl_ids) 
  dge <- removeReplicates(dge, rm=rmGroup)
  dgel <- list(dge)
  names(dgel)[1] <- "All"
  
}






#generate output text files
saveRDS(dgel, file=file.path(outdir, paste0(method, ".rds")))





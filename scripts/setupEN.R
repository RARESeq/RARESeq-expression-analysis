required.packages <- c("tidyverse",
                       "data.table",
                       "glmnet",
                       "caret",
                       "limma")
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
if (length(args) == 0) {
  stop("Expected arguments: 1) dge object 2) mutation calls file (optional)") 
}
dge_file <- args[1]
if (length(args)>1) { calls_file <- args[2] }
case_grps=c("LUAD")
ctrl_grps=c("Control")
nfolds <- 10
k <- 2
q <- 0.6
balanced <- TRUE
set.seed(15000)


####required####
dir <- "../gene_lists"
tcga_luad_df=fread(file.path(dir, "TCGA_LUAD_feature_weights.txt"), data.table=F)
tcga_luad_df$weight=((tcga_luad_df$Cancer-min(tcga_luad_df$Cancer)) + abs(tcga_luad_df$CancervscfRNA_LFC)) * sign(tcga_luad_df$CancervscfRNA_LFC)
ss_df=fread(file.path(dir, "sex_specific_genes.txt"), data.table=F)
if (length(args)>1) { varid_df=fread(file.path(dir, "NCCN_VARID.txt"), data.table=F) }



#create output dir
outdir <- file.path(dirname(dge_file),"caret")
dir.create(outdir, showWarnings=F, recursive=T)



#read in data
dgel <- readRDS(dge_file)
dge <- dgel[["All"]]
info <- dge$samples
genes <- dge$genes
gene_ids <- genes$gene_id_simplified



#select samples
selected_samples <- info[which(info$group %in% c(case_grps, ctrl_grps)),,drop=FALSE]
if ("Reference" %in% colnames(info)) { selected_samples <- selected_samples[is.na(selected_samples$Reference),,drop=FALSE] }
if ("GenotyperReference" %in% colnames(info)) { selected_samples <- selected_samples[is.na(selected_samples$GenotyperReference),,drop=FALSE] }
selected_samples$true_class <- 1*(selected_samples$group %in% case_grps)
ids <- selected_samples$ID
#print(ids)



#select features
selected_genes <- tcga_luad_df %>% 
  filter(gene_id_simplified %in% gene_ids) %>% 
  filter(Cancer_K==k | cfRNA_K==k) %>%
  filter(CancervscfRNA_LFC <= quantile(tcga_luad_df$CancervscfRNA_LFC, 1-q) |
         CancervscfRNA_LFC >= quantile(tcga_luad_df$CancervscfRNA_LFC, q)) %>%
  filter(!gene_id_simplified %in% ss_df$gene_id_simplified) %>% 
  select(gene_id_simplified, gene_symbol, weight)
colnames(selected_genes)[1]="feature_id"
selected_genes$weight=1 / (selected_genes$weight + 20)

gene_ids <- selected_genes$feature_id
#print(gene_ids)


#remove features with zero variance
#and correct for selector differences using limma removeBatchEffect()
exp <- dge$log2nx[match(gene_ids, rownames(dge$log2nx)),
                  match(ids, colnames(dge$log2nx))]

#zv <- apply(exp, 1, var) == 0 
nzv <- apply(exp, 1, var) < 1e-3
exp <- exp[!nzv,]
gene_ids <- gene_ids[!nzv]
selected_genes <- selected_genes[!nzv,]

batch=selected_samples$Selector
mod=model.matrix(~as.factor(true_class), data=selected_samples)
exp=removeBatchEffect(exp,batch=batch,design=mod)








#read in mut data
if (exists("calls_file")) {
  calls=fread(calls_file,data.table=FALSE)
  
  if (!"ID" %in% colnames(calls)) {
    stop("Warning! Genotyping file must contain ID column header.")
  } else if (length(setdiff(ids, calls$ID)) > 0) {
    print(setdiff(ids, calls$ID))
    stop("Warning! Genotyping file missing sample IDs.")
  } else {
    
    if (!"VarID" %in% colnames(calls)) {
      varidsbypos=varid_df[!is.na(varid_df$START) & !is.na(varid_df$END),]
      varidsclean=varid_df[is.na(varid_df$START) & is.na(varid_df$END), c("CHR","POS","REF","VAR","GENE","VarID")]
      
      calls=calls %>% full_join(varidsclean)
      for (i in 1:nrow(varidsbypos)) {
        id=varidsbypos$VarID[i]
        gene=varidsbypos$GENE[i]
        start=varidsbypos$START[i]
        end=varidsbypos$END[i]
        calls$VarID[calls$GENE==gene & ((calls$POS>start & calls$POS<end) | is.na(calls$POS))]=id
      }
    }
    
    fusions <- calls$VarID[which(calls$CallType=="Fusion")]
    indels <- calls$VarID[which(calls$CallType=="Indel")]
    snvs <- calls$VarID[which(calls$CallType=="SNV")]
    
    calls$AF[!calls$Called]=0
    calls$AF[which(calls$AF>1)]=1
    calls=calls %>%
      filter(!is.na(ID)) %>%
      filter(!is.na(VarID)) %>%
      group_by(ID, VarID) %>%
      summarize(maxAF=max(AF, na.rm=TRUE)) %>%
      spread(VarID, maxAF) %>%
      as.data.frame()
    calls[is.na(calls)]=0
    rownames(calls)=calls[,1,drop=TRUE]
    calls=calls[,-1,drop=FALSE]
    
    calls$ANYEGFR=apply(calls[,grepl("EGFR",colnames(calls)),drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANYKRAS=apply(calls[,grepl("KRAS",colnames(calls)),drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANYMET=apply(calls[,grepl("MET",colnames(calls)),drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANYFUSION=apply(calls[,colnames(calls) %in% fusions,drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANYINDEL=apply(calls[,colnames(calls) %in% indels,drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANYSNV=apply(calls[,colnames(calls) %in% snvs,drop=FALSE], 1, function(x) ifelse(sum(x)>0,1,0))
    calls$ANY=apply(calls, 1, function(x) ifelse(sum(x)>0,1,0))
    
  }
  
} else {
  calls <- matrix(nrow=length(ids),ncol=0)
  rownames(calls) <- ids
  
}




#set xall
#combine expression & muts features
xall=cbind(t(exp[,match(ids, colnames(exp)),drop=FALSE]),
           calls[match(ids, rownames(calls)),,drop=FALSE])

features=as.data.frame(rbindlist(list(selected_genes, 
                                      data.frame(feature_id=colnames(calls), 
                                                 gene_symbol=colnames(calls), 
                                                 weight=rep(min(selected_genes$weight),ncol(calls)))),
                                 fill=TRUE))
rownames(features) <- features$feature_id
  

#set yall
yall=selected_samples$true_class


#split samples into train and test sets
train <- as.logical(selected_samples$EN=="train")
xtrain=xall[train,]
ytrain=yall[train]
xval=xall[!train,]
yval=yall[!train]



#split samples into random folds
balanced.folds <- function(y, nfolds=min(min(table(y)), 10)) {
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)     
  nfolds <- max(nfolds, 2) # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)  # nice we to get the ids in a list, split by class
  
  # make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds)*nfolds, length(totals))
  for (i in seq(totals)) {
    if (length(yids[[i]])>1) { bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
    if (length(yids[[i]])==1) { bigmat[seq(totals[i]), i] <- yids[[i]] }
  }
  # reshape the matrix
  smallmat <- matrix(bigmat, nrow=nfolds)
  
  #now do a clever sort to mix up the NAs
  #smallmat <- permute.rows(t(smallmat))
  smallmat <- t(smallmat)
  smallmat <- smallmat[sample(nrow(smallmat)),]
  
  #now a clever unlisting
  #note the "clever" unlist doesn't work when there are no NAs
  apply(smallmat, 2, function(x) x[!is.na(x)])
  res <- vector("list", nfolds)
  for(j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  return(res)
}

if (balanced) {
  folds=balanced.folds(ytrain, nfolds)
  foldss=rep(NA, length(ytrain))
  for (i in 1:nfolds) { foldss[folds[[i]]]=i }
  folds=foldss
} else {
  folds=sample(1:nfolds, size=length(ytrain), replace=TRUE)
}



#center and scale all features
#zv <- apply(xtrain, 2, var) == 0 
nzv <- apply(xtrain, 2, var) < 1e-3
xtrain <- xtrain[,!nzv,drop=FALSE]
xtrainpp <- preProcess(xtrain, method=c("center", "scale"))
xtrain <- predict(xtrainpp, newdata=xtrain)

xval <- xval[,!nzv,drop=FALSE]
features <- features[!nzv,,drop=FALSE]



#write out train and test data
dge <- list(genes=selected_genes,
            samples=selected_samples,
            xtrain=data.matrix(xtrain),
            ytrain=ytrain,
            xval=data.matrix(xval),
            yval=yval,
            folds=folds,
            features=features,
            preProcess=xtrainpp)
dgel <- list()
dgel[[1]] <- dge
names(dgel)[1] <- "All"
saveRDS(dgel, file.path(outdir, gsub(".rds", "_forCaret.rds", basename(dge_file))))





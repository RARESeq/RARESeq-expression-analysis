required.packages <- c("tidyverse",
                       "data.table",
                       "glmnet",
                       "caret",
                       "limma",
                       "pROC")
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
  stop("Expected arguments: 1) dge object 2) caret model output")
}
dge_file <- args[1]
mdl_file <- args[2]
set.seed(15000)


#create output dir
outdir <- dirname(dge_file)


#read in data
dgel <- readRDS(dge_file)
dge <- dgel[["All"]]
info <- dge$samples
xval <- dge$xval
yval <- dge$yval
pp <- dge$preProcess


# pre-process val data
xval <- predict(pp, newdata=xval)


#read in final model
model <- readRDS(mdl_file)


# make predictions
sample_pred <- data.frame(ID=rownames(xval), 
                          true_class=yval,
                          response=predict(model, newdata=xval, type="prob", s="lambda.min")[,2])
write.table(sample_pred, file.path(outdir, "model_sample_predictions_validation.txt"), row.names=FALSE, quote=FALSE, sep="\t")


val_roc <- roc(sample_pred$true_class, sample_pred$response)
aucs <- data.frame(train="Validation", 
                   auc=val_roc$auc)
write.table(aucs, file.path(outdir, "model_metrics_validation.txt"), row.names=FALSE, quote=FALSE, sep="\t")




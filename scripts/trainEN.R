required.packages <- c("tidyverse",
                       "data.table",
                       "glmnet",
                       "caret",
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
  stop("Expected arguments: 1) dge object")
}
dge_file <- args[1]
stand <- FALSE
set.seed(15000)


#create output dir
outdir <- dirname(dge_file)



#read in input file
dgel <- readRDS(dge_file)
dge <- dgel[["All"]]
xtrain <- dge$xtrain
ytrain <- dge$ytrain
features <- dge$features
featurewts <- features$weight
folds <- dge$folds
nfolds <- length(unique(folds))



#establish hyperparameter tuning grid
#alphas <- 1
alphas <- seq(0, 1, length=11)
len <- 10
init <- glmnet(xtrain, ytrain, family="binomial", alpha=0.5, nlambda=len)
lambdas <- unique(init$lambda)
lambdas <- lambdas[-c(1, length(lambdas))]
lambdas <- lambdas[1:min(length(lambdas), len)]
fitGrid <- expand.grid(alpha=alphas, lambda=lambdas) # creates grid of desired alphas and lambdas
fitControl <- trainControl("cv", number=5, classProbs=TRUE, summaryFunction=twoClassSummary) # runs 5-fold CV for every combo of alpha/lambda in tuneGrid



#train model 
sample_pred <- list()
feature_sel <- list()
metrics <- list()
rocobj <- list()
if (nfolds > 1) {
  for (i in 1:nfolds) {
    print(paste0("Fold ", i))
    
    xout <- xtrain[folds!=i,,drop=FALSE]
    yout <- ytrain[folds!=i]
    xin <- xtrain[folds==i,,drop=FALSE]
    yin <- ytrain[folds==i]
    
    model <- train(xout, make.names(yout), method="glmnet", family="binomial", metric="ROC", 
                   penalty.factor=featurewts, lower.limits=0, standardize=stand, 
                   trControl=fitControl, tuneGrid=fitGrid)

    # select best model (i.e. highest accuracy)
    alpha <- model$bestTune$alpha
    lambda <- model$bestTune$lambda
    
    # extract selected features and associated coefficients 
    feature_coef <- as.numeric(coef(model$finalModel, lambda)[-1])
    feature_i <- feature_coef!=0
    if(sum(feature_i)>0) {
      feature_sel[[i]] <- data.frame(fold=i,
                                     features[feature_i,,drop=FALSE], 
                                     coef=feature_coef[feature_i])
      
      # extract info for samples & assess performance
      sample_pred[[i]] <- data.frame(fold=i,
                                     ID=rownames(xin),
                                     true_class=yin,
                                     response=predict(model, xin, type="prob", s=lambda)[,2]) 

    }
  }
}



# generate final model using FULL training set
# above used to estimate the performance of this model
print("All folds")

model <- train(xtrain, make.names(ytrain), method="glmnet", family="binomial", metric="ROC", 
               penalty.factor=featurewts, lower.limits=0, standardize=stand, 
               trControl=fitControl, tuneGrid=fitGrid)

# select best model (i.e. highest accuracy)
alpha <- model$bestTune$alpha
lambda <- model$bestTune$lambda

# extract selected features and associated coefficients 
feature_coef <- as.numeric(coef(model$finalModel, lambda)[-1])
feature_i <- feature_coef!=0
feature_sel[[length(feature_sel)+1]] <- data.frame(fold="full",
                                                   features[feature_i,,drop=FALSE], 
                                                   coef=feature_coef[feature_i])

# extract info for samples & assess performance
sample_pred[[length(sample_pred)+1]] <- data.frame(fold="full",
                                                   ID=rownames(xtrain),
                                                   true_class=ytrain,
                                                   response=predict(model, xtrain, type="prob", s=lambda)[,2]) 






# write out model performance metrics
feature_sel <- do.call(rbind, feature_sel)
write.table(feature_sel %>% filter(fold=="full"), file=file.path(outdir, "model_selected_features.txt"), quote=F, row.names=F, sep="\t")


sample_pred <- do.call(rbind, sample_pred)
write.table(sample_pred %>% filter(fold=="full"), file=file.path(outdir, "model_sample_predictions.txt"), quote=F, row.names=F, sep="\t")

trainroc <- roc(sample_pred$true_class[sample_pred$fold=="full"], sample_pred$response[sample_pred$fold=="full"])
cvroc <- roc(sample_pred$true_class[sample_pred$fold!="full"], sample_pred$response[sample_pred$fold!="full"])

metrics <- data.frame(train=c("Train", "CV"), auc=c(trainroc$auc, cvroc$auc))
write.table(metrics, file=file.path(outdir, "model_metrics.txt"), quote=F, row.names=F, sep="\t")




# write out final model
saveRDS(model, file=file.path(outdir, "model.rds"))



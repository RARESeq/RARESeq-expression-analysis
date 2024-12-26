#load required packages
required.packages <- c("tidyverse",
                       "data.table",
                       "tximport")
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
  stop("Expected arguments: 1) directory with RSEM expression files 2) sample info file ")
}
indir <- args[1]
sample_file <- args[2]
set.seed(100)


#create output dir
outdir <- file.path("output")
dir.create(outdir, showWarnings=F, recursive=T)


#read in input file
#requires absolute paths
samples <- fread(sample_file, data.table=F)


#check input file contains correct cols
if (is.null(samples$ID)) {
  stop("Your input file must contain a column named ID (case-sensitive)")
} else {
  files <- paste0(indir, "/Sample_", samples$ID, "_cfrna.genes.results")
  names(files) <- samples$ID

  print(paste0("Found ", sum(file.exists(files)), " samples of ", length(files), " total samples"))
  print(paste0("Missing samples: ", paste(samples$ID[!file.exists(files)], collapse=",")))
}



#read in expression data
txi <- tximport(files[file.exists(files)], type="rsem", txIn=FALSE, txOut=FALSE)
txi$length[which(txi$length==0)] <- 1.0 # required for deseq to work



#remove gene id version numbers
ids=unlist(lapply(strsplit(rownames(txi$counts), "\\."), function(x) x[1]))
rownames(txi$counts)=ids
rownames(txi$abundance)=ids
rownames(txi$length)=ids


#write out txi object
saveRDS(txi, file=file.path(outdir, "tximport.rds"))

#write out counts matrix
write.table(data.frame(gene_id_simplified=ids, txi$counts), file.path(outdir, "gene_expected_count.txt"), sep="\t", quote=F, row.names=F, col.names=T)



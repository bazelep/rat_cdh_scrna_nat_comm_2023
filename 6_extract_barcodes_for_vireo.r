#Summary: Extract barcodes for creating per-cell sam files
#Conda: cdh


library(Seurat)
library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)

source("utilities.r")



#need to identify and remove doublets

#create folder for donor analysis
donors.dir = "../OutputData/Donors"
if(!dir.exists(donors.dir)){dir.create(donors.dir)}



##Filter settings
f.genes = 500
f.umi = 500
f.mt = 0.2
f.complexity = 0.8
f.rb = 0.2

filter.string = paste0(
  "Filtered.G",f.genes,
  ".U",f.umi,
  ".M",f.mt,
  ".R",f.rb,
  ".C",f.complexity
)

filtered_seurat = readRDS(paste0("../RData/",filter.string,".rds"))
samples = unique(filtered_seurat$sample)

plan("multisession", workers = 40)

barcodes.res = foreach(s = samples) %do% {

  barcodes = subset(filtered_seurat@meta.data, sample == s)$barcode

  sample.res = foreach(barcode = barcodes) %dopar% {
    barcode.parts = strsplit(barcode, "_")[[1]]
    return(barcode.parts[length(barcode.parts)])
  }
  return(unlist(sample.res))
}
names(barcodes.res) = samples



sapply(
  names(barcodes.res), 
  function(s){
    write.table(
      barcodes.res[[s]], 
      file = paste0(donors.dir,"/",s,".barcodes.csv"), 
      sep=",", 
      row.names=F,
      quote=F,
      col.names=F
    )
  }
)


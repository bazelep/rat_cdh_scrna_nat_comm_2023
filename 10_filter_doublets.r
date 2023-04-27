#Summary: Load Vireo results and identify and remove doublets
#Conda: cdh

library(Seurat)
library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)

source("utilities.r")


donors.dir = "../OutputData/Donors/"
pc.dims = 30

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


#split object by sample
samples.list = SplitObject(filtered_seurat, split.by = "sample")

#function to strip off original barcode from current one with sample name affixed
get.original.barcodes = function(barcode){
    parts = strsplit(barcode,"_")[[1]]
    return(parts[length(parts)])
}


#assign Vireo doublets and donors to split samples
sample.res = foreach(sample = names(samples.list)) %do% {
    print(sample)

    #barcode needs updated sample name
    #grab original barcode
    samples.list[[sample]]$original.barcode = sapply(samples.list[[sample]]$barcode,get.original.barcodes)

    #combine newer sample name with original barcode
    samples.list[[sample]]$sample.barcode = paste0(
        samples.list[[sample]]$sample,
        "_",
        samples.list[[sample]]$original.barcode
    )

    #load vireo donor assignments    
    sample.donors.file = paste0(donors.dir,"/",sample,"/vireo.out/donor_ids.tsv")
    print(sample.donors.file)

    sample.donors = read.table(sample.donors.file, header = T, sep = "\t")
    sample.donors$barcode = gsub("\\.","_",sample.donors$cell)
    print(table(sample.donors$best_singlet))
    donors = unique(sample.donors$best_singlet)

    sample.doublets = subset(sample.donors, donor_id == "doublet")$barcode

    sample.barcodes = as.vector(samples.list[[sample]]$sample.barcode)

    donor.res = foreach(donor = donors) %do% {
        print(donor)

        donor.rows = subset(sample.donors, best_singlet == donor)

        barcode.res = foreach(sample.bc = sample.barcodes) %do% {

            sample.bc.row = rownames(subset(samples.list[[sample]]@meta.data, sample.barcode == sample.bc))

            if(sample.bc %in% donor.rows$barcode){
                samples.list[[sample]]@meta.data[sample.bc.row,"Donor"] = donor
            }

            if(sample.bc %in% sample.doublets){
                samples.list[[sample]]@meta.data[sample.bc.row,"Doublet"] = 1
            } else {
                samples.list[[sample]]@meta.data[sample.bc.row,"Doublet"] = 0
            }
        }
    }
}


#combine updated seurat list
filtered_seurat_withdoublets = merge(samples.list[[1]],samples.list[-1])


xtabs(~ sample + Donor, data = filtered_seurat_withdoublets@meta.data, addNA = T)
xtabs(~ sample + Doublet, data = filtered_seurat_withdoublets@meta.data, addNA = T)


#create cell count table, updated with Donors,Doublets
cell.counts = store.counts(so = filtered_seurat_withdoublets, stage = paste0(filter.string,".WithDoublets"), ID.col = "sample")
write.csv(cell.counts, file = paste0("../OutputData/QC/",filter.string,".WithDoublets.csv"), row.names = F)






#create UMAPs of doublet locations

#grab doublet cells for highlighting in UMAP
doublets = rownames(subset(filtered_seurat_withdoublets, Doublet == 1)@meta.data)

#need normalized data - will do this in next script as well
#normalize, regress out MT variation, find variable features, scale, run pca
# no established cell-cycle genes for rat, so skipping correction
filtered_seurat_withdoublets = SCTransform(
    filtered_seurat_withdoublets, 
    vars.to.regress = c("mitoRatio"), 
    variable.features.n = 3000,
    ncells = 5000,
    method = "qpoisson",
    batch_var = "sample"
)
filtered_seurat_withdoublets = RunPCA(filtered_seurat_withdoublets, assay = "SCT")
filtered_seurat_withdoublets = RunUMAP(filtered_seurat_withdoublets, dims = 1:pc.dims,reduction = "pca", return.model = T)


#highlight doublet cells in UMAP
pdf(file=paste0("../Graphs/QC/",filter.string,".Doublets.pdf"))
p = DimPlot(
    filtered_seurat_withdoublets, 
    cells.highlight = doublets,
    split.by = "sample"
)
print(p)
dev.off()

#remove cells with NAs for Donor assignment - shouldn't be any
dim(filtered_seurat_withdoublets)
filtered_seurat_withdoublets = filtered_seurat_withdoublets[,!is.na(filtered_seurat_withdoublets@meta.data$Donor)]
dim(filtered_seurat_withdoublets)

#remove doublets
filtered_seurat_nodoublets = subset(filtered_seurat_withdoublets, Doublet == 0)
dim(filtered_seurat_nodoublets)

table(filtered_seurat_nodoublets$sample)

#store filtered counts
cell.counts = store.counts(so = filtered_seurat_nodoublets, stage = paste0(filter.string,".NoDoublets"), ID.col = "sample")
write.csv(cell.counts, file = paste0("../OutputData/QC/",filter.string,".NoDoublets.csv"), row.names = F)

#save filtered object
saveRDS(filtered_seurat_nodoublets, file = paste0("../RData/",filter.string,".donors.nodoublets.rds"))



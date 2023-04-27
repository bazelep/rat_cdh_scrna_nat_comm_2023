#Summary: Normalize counts and correct for nuisance variation - re-normalize after removing doublets
#Conda: cdh

library(Seurat)


set.seed(100)
options(width=150)


#set filter string for filtered seurat object
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


#load data
filtered_seurat = readRDS(file = paste0("../RData/",filter.string,".donors.nodoublets.rds"))

so.list = SplitObject(filtered_seurat, split.by = "sample")

#normalize, regress out MT variation, find variable features, scale
# no established cell-cycle genes for rat, so skipping correction

for(samp in names(so.list)){
    print(samp)
    so.list[[samp]] = SCTransform(
        so.list[[samp]], 
        vars.to.regress = c("mitoRatio"), 
        variable.features.n = 3000,
        ncells = 5000,
        method = "qpoisson"
    )
}

#save normalized object
saveRDS(
    so.list, 
    file = paste0("../RData/",filter.string,".donors.nodoublets.SCTnormalized_UMIcount,MTratio.rds")
)







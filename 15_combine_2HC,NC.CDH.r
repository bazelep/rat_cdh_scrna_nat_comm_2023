#Summary: Combine integrated healthy controls with nitrofen control and nitrofen cdh
#Conda: cdh


library(Seurat)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=100)



healthy.samples.file = "../RData/2HealthyControls.integrated.clustered.rds"
nitrofen.samples.file = "../RData/2NitrofenSamples.clustered.rds"



healthy.samples = readRDS(healthy.samples.file)
nitrofen.samples = readRDS(nitrofen.samples.file)

#assign experimental group status
healthy.samples$Status = "2HC"
nitrofen.samples[["E215_NitrofenControl"]]$Status = "NC"
nitrofen.samples[["E215_NitrofenCDH"]]$Status = "CDH"


so.list = c(healthy.samples, nitrofen.samples)
names(so.list)[1] = "2HC"


#use status to name so.list
names(so.list) = sapply(
    names(so.list),
        function(x){
            unique(so.list[[x]]$Status)[1]
        }
    )



saveRDS(so.list, file = "../RData/combined.samples.clustered.rds")


load("RData/DiffBind2/li.diff.RData.minOverlap2")

dba.report(li.diff[["female"]],method=DBA_DESEQ2) %>% nrow
dba.report(li.diff[["female"]],method=DBA_EDGER) %>% nrow

dba.report(li.diff[["male"]],method=DBA_DESEQ2) %>% nrow
dba.report(li.diff[["male"]],method=DBA_EDGER) %>% nrow

li.diff[["male"]] %>% names
li.diff[["male"]][["filter"]]
li.diff[["male"]][["attributes"]]
li.diff[["male"]][["merged"]] %>% nrow
li.diff[["male"]][["merged"]][1:10,]
li.diff[["male"]][["peaks"]][[1]][1:10,]

li.diff[["male"]][["contrasts"]][[1]] %>% names

li.diff[["male"]][["contrasts"]][[1]][["group1"]]
li.diff[["male"]][["contrasts"]][[1]][["group2"]]
li.diff[["male"]][["contrasts"]][[1]][["name1"]]
li.diff[["male"]][["contrasts"]][[1]][["name2"]]

li.diff[["male"]][["contrasts"]][[1]][["edgeR"]] %>% names
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$design
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$common.dispersion
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$pseudo.counts[1:10,]
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$AveLogCPM[1:10] 
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$pseudo.lib.size
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$bFullLibrarySize
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$bSubControl
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$span
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$prior.df
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$GLM %>% names
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$LRT %>% names
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$GLM
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$LRT
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$GLM$samples
li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$LRT$samples

dt.table<-li.diff[["male"]][["contrasts"]][[1]][["edgeR"]]$LRT$table %>% data.table
dt.table[PValue<0.05]
dt.table[,Adj.PValue:=p.adjust(PValue,method="BH")]
dt.table[Adj.PValue<0.05]


##
## DESeq2
##
li.diff[["male"]][["contrasts"]][[1]][["DESeq2"]] %>% names
li.diff[["male"]][["contrasts"]][[1]][["DESeq2"]]$bFullLibrarySize
li.diff[["male"]][["contrasts"]][[1]][["DESeq2"]]$bSubControl
li.diff[["male"]][["contrasts"]][[1]][["DESeq2"]]$facs
dt.table<-li.diff[["male"]][["contrasts"]][[1]][["DESeq2"]]$de %>% data.table
dt.table[padj<0.05]
dt.table[pval<0.05]


foo.edgeR<-
    li.diff[["male"]][["contrasts"]][[1]]$edgeR
class(foo.edgeR)

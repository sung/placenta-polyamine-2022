library(DiffBind)
load("RData/DiffBind2/dbo.cnt.RData")
dba.show(dbo.cnt)
dba.normalize(dbo.cnt, bRetrieve=T,method=DBA_ALL_METHODS)

li.diff<-list() # list of DBA

########################
## pre-version 3 mode ##
########################
dbo.cnt.old<-dbo.cnt
dbo.cnt.old$config$design<-FALSE
li.diff[['pre3']][["female"]]<-dba(dbo.cnt.old,dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3']][["male"]]<- dba(dbo.cnt.old,dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3']]<-lapply(li.diff[['pre3']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

###############################################################################
## pre-version3 with 'RiP' normalization with native method (i.e. TMM & RLE) ##
###############################################################################
dbo.cnt.old<-dbo.cnt
li.diff[['pre3-rip']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-rip']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-rip']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-rip']]<-lapply(li.diff[['pre3-rip']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_PEAKREADS)
li.diff[['pre3-rip']]<-lapply(li.diff[['pre3-rip']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

li.diff[['pre3-rip-pair']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-rip-pair']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-rip-pair']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-rip-pair']]<-lapply(li.diff[['pre3-rip-pair']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_PEAKREADS)
li.diff[['pre3-rip-pair']]<-lapply(li.diff[['pre3-rip-pair']], dba.contrast, design="~Replicate+Factor" )
li.diff[['pre3-rip-pair']]<-lapply(li.diff[['pre3-rip-pair']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

############################################################################
## pre-version3 with 'default' normalization ('full' size & 'lib' method) ##
############################################################################
dbo.cnt.old<-dbo.cnt
li.diff[['pre3-default']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-default']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-default']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-default']]<-lapply(li.diff[['pre3-default']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

li.diff[['pre3-default-pair']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-default-pair']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-default-pair']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-default-pair']]<-lapply(li.diff[['pre3-default-pair']], dba.contrast, design="~Replicate+Factor" )
li.diff[['pre3-default-pair']]<-lapply(li.diff[['pre3-default-pair']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

######################################################################################
## pre-version3 with 'background' normalization with native method (i.e. TMM & RLE) ##
######################################################################################
dbo.cnt.old<-dbo.cnt
li.diff[['pre3-bg']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-bg']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-bg']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-bg']]<-lapply(li.diff[['pre3-bg']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['pre3-bg']]<-lapply(li.diff[['pre3-bg']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

li.diff[['pre3-bg-pair']][["female"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$F) # isa DBA
li.diff[['pre3-bg-pair']][["male"]]<- dba(dbo.cnt.old, dbo.cnt.old$masks$M) # isa DBA
li.diff[['pre3-bg-pair']][["male-n10"]]<-dba(dbo.cnt.old, dbo.cnt.old$masks$M & !dbo.cnt.old$mask$Replicate.2) # male-n10
li.diff[['pre3-bg-pair']]<-lapply(li.diff[['pre3-bg-pair']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['pre3-bg-pair']]<-lapply(li.diff[['pre3-bg-pair']], dba.contrast, design="~Replicate+Factor")
li.diff[['pre3-bg-pair']]<-lapply(li.diff[['pre3-bg-pair']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

#####################################
## Ver3 summit count with default ###
#####################################
load("RData/DiffBind3.0.13/dbo.RData")
load("RData/DiffBind3.0.13/dbo.con.RData")
consensus_peaks <- dba.peakset(dbo.con, bRetrieve=TRUE) # isa GRanges if occupancy matrix (log10p0val based)
dbo.cnt.summit<-dba.count(dbo, peaks=consensus_peaks, bSubControl=TRUE) 
dba.show(dbo.cnt.summit)

load("RData/DiffBind3.0.13/dbo.con.RData")

li.diff[['ver3-summit-default']][["female"]]<-dba(dbo.cnt.summit, dbo.cnt.summit$masks$F) # isa DBA
li.diff[['ver3-summit-default']][["male"]]<- dba(dbo.cnt.summit, dbo.cnt.summit$masks$M) # isa DBA
li.diff[['ver3-summit-default']]<-lapply(li.diff[['ver3-summit-default']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

########################################################
## Ver3 summit count with 'background' normalization  ##
## with native method (i.e. TMM & RLE)                ##
########################################################
li.diff[['ver3-summit-bg']][["female"]]<-dba(dbo.cnt.summit, dbo.cnt.summit$masks$F) # isa DBA
li.diff[['ver3-summit-bg']][["male"]]<- dba(dbo.cnt.summit, dbo.cnt.summit$masks$M) # isa DBA
li.diff[['ver3-summit-bg']]<-lapply(li.diff[['ver3-summit-bg']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['ver3-summit-bg']]<-lapply(li.diff[['ver3-summit-bg']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

li.diff[['ver3-summit-bg-pair']][["female"]]<-dba(dbo.cnt.summit, dbo.cnt.summit$masks$F) # isa DBA
li.diff[['ver3-summit-bg-pair']][["male"]]<- dba(dbo.cnt.summit, dbo.cnt.summit$masks$M) # isa DBA
li.diff[['ver3-summit-bg-pair']]<-lapply(li.diff[['ver3-summit-bg-pair']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['ver3-summit-bg-pair']]<-lapply(li.diff[['ver3-summit-bg-pair']], dba.contrast, design="~Replicate+Factor")
li.diff[['ver3-summit-bg-pair']]<-lapply(li.diff[['ver3-summit-bg-pair']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

###############################
## Ver3 peaks with 'default' ##
###############################
dbo.cnt.peak<-dba.count(dbo, peaks=consensus_peaks, bSubControl=TRUE, summits=FALSE,bUseSummarizeOverlaps=FALSE,filter=0,minCount=1,mqpQCth=10) # this applies all samples using the consensus peaks (see above)
dba.show(dbo.cnt.peak)

li.diff[['ver3-pk-default']][["female"]]<-dba(dbo.cnt.peak, dbo.cnt.peak$masks$F) # isa DBA
li.diff[['ver3-pk-default']][["male"]]<- dba(dbo.cnt.peak, dbo.cnt.peak$masks$M) # isa DBA
#li.diff[['ver3-pk-default']]<-lapply(li.diff[['ver3-pk-default']], dba.contrast, design=FALSE) # pre3 mode
li.diff[['ver3-pk-default']]<-lapply(li.diff[['ver3-pk-default']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

#################################
## Ver3-peak with 'background' ##
#################################
li.diff[['ver3-pk-bg']][["female"]]<-dba(dbo.cnt.peak, dbo.cnt.peak$masks$F) # isa DBA
li.diff[['ver3-pk-bg']][["male"]]<- dba(dbo.cnt.peak, dbo.cnt.peak$masks$M) # isa DBA
li.diff[['ver3-pk-bg']]<-lapply(li.diff[['ver3-pk-bg']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['ver3-pk-bg']]<-lapply(li.diff[['ver3-pk-bg']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

li.diff[['ver3-pk-bg-pair']][["female"]]<-dba(dbo.cnt.peak, dbo.cnt.peak$masks$F) # isa DBA
li.diff[['ver3-pk-bg-pair']][["male"]]<- dba(dbo.cnt.peak, dbo.cnt.peak$masks$M) # isa DBA
li.diff[['ver3-pk-bg-pair']]<-lapply(li.diff[['ver3-pk-bg-pair']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE)
li.diff[['ver3-pk-bg-pair']]<-lapply(li.diff[['ver3-pk-bg-pair']], dba.contrast, design="~Replicate+Factor")
li.diff[['ver3-pk-bg-pair']]<-lapply(li.diff[['ver3-pk-bg-pair']], dba.analyze, method=DBA_ALL_METHODS,bBlacklist=F,bGreylist=F)

##
## NO. of DBR
##
lapply(li.diff[['pre3']], dba.show, bContrasts=TRUE)

lapply(li.diff[['pre3-rip']], dba.show, bContrasts=TRUE)
lapply(li.diff[['pre3-rip-pair']], dba.show, bContrasts=TRUE)

lapply(li.diff[['pre3-default']], dba.show, bContrasts=TRUE)
lapply(li.diff[['pre3-default-pair']], dba.show, bContrasts=TRUE)

lapply(li.diff[['pre3-bg']], dba.show, bContrasts=TRUE)
lapply(li.diff[['pre3-bg-pair']], dba.show, bContrasts=TRUE)

lapply(li.diff[['pre3-rip']][c('male','male-n10')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
lapply(li.diff[['pre3-rip-pair']][c('male','male-n10')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)

lapply(li.diff[['pre3-default']][c('male','male-n10')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
lapply(li.diff[['pre3-default-pair']][c('male','male-n10')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)

lapply(li.diff[['pre3-bg']][c('male','male-n10')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
lapply(li.diff[['pre3-bg-pair']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)

#
lapply(li.diff[['ver3-pk-default']], dba.show, bContrasts=TRUE)
lapply(li.diff[['ver3-pk-bg']], dba.show, bContrasts=TRUE)
lapply(li.diff[['ver3-pk-bg-pair']], dba.show, bContrasts=TRUE)

lapply(li.diff[['ver3-pk-bg']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
lapply(li.diff[['ver3-pk-bg-pair']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)

#
lapply(li.diff[['ver3-summit-default']], dba.show, bContrasts=TRUE)
lapply(li.diff[['ver3-summit-bg']], dba.show, bContrasts=TRUE)
lapply(li.diff[['ver3-summit-bg-pair']], dba.show, bContrasts=TRUE)

lapply(li.diff[['ver3-summit-bg']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
lapply(li.diff[['ver3-summit-bg-pair']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)

##
## design
##
lapply(li.diff[['pre3-rip']], dba.show, bDesign=T)
lapply(li.diff[['pre3-rip-pair']], dba.show, bDesign=T)


##
## Normalization parameter
##
lapply(li.diff[['pre3']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']]

lapply(li.diff[['pre3-rip']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]

lapply(li.diff[['pre3-default']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]
lapply(li.diff[['pre3-default']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['DESeq2']]
lapply(li.diff[['pre3-default-pair']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]

lapply(li.diff[['pre3-bg']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]
lapply(li.diff[['pre3-bg']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['DESeq2']]
lapply(li.diff[['pre3-bg-pair']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['DESeq2']]

lapply(li.diff[['ver3-pk-default']], dba.normalize,bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]
lapply(li.diff[['ver3-pk-bg']], dba.normalize,bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['DESeq2']]
lapply(li.diff[['ver3-pk-bg-pair']], dba.normalize,bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['DESeq2']]

lapply(li.diff[['ver3-summit-default']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]
lapply(li.diff[['ver3-summit-bg']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[['male']][['edgeR']]


##
## DESeq2 details
##
li.diff[['pre3-bg']][['male']]$config %>% names
li.diff[['ver3-pk-bg']][['male']]$config %>% names

dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_DESEQ2) # isa DESeqDataSet
dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_DESEQ2) %>% design
dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_DESEQ2) %>% rowData
dba.analyze(li.diff[['ver3-pk-bg']][['male']], bRetrieveAnalysis=DBA_DESEQ2) %>% rowData

##
## edgeR details
##
dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_EDGER)
dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_EDGER) %>% names
dba.analyze(li.diff[['pre3-bg']][['male']], bRetrieveAnalysis=DBA_EDGER)$method

dba.analyze(li.diff[['ver3-pk-bg']][['male']], bRetrieveAnalysis=DBA_EDGER) %>% class

###
### Plot
###

dba.plotMA(li.diff[['pre3']][['male']],bNormalized=F,dotSize=2,method=DBA_DESEQ2)
dba.plotMA(li.diff[['pre3']][['male']],bNormalized=T,dotSize=2,method=DBA_DESEQ2)

dba.plotMA(li.diff[['pre3']][['female']],bNormalized=F,dotSize=1)

###
dba.peakset( dba.count(dbo.cnt.old, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE)[1:10,1:6]
dba.peakset( dba.count(dbo.cnt.peak, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE)[1:10,1:6]
dba.peakset( dba.count(dbo.cnt.summit, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE)[1:10,1:6]


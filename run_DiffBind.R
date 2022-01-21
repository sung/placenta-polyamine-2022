##
## Samples
##
dt.samples<-fread("~/results/SLX-16110-16111.Homo_sapiens.v1/Meta/samples.diffbind.csv")
dt.samples.narrow<-fread("~/results/SLX-16110-16111.Homo_sapiens.v1/Meta/samples.diffbind.narrowPeak.csv")
names(dt.samples)

#############
## SET ENV ##
#############
library(BiocParallel)
chip.env<-list()
chip.env['FR.size']=148 # default: 125; `d`? in macs2
chip.env['ref']='hg38' 
chip.env['MAPQ']=10
chip.env['blacklist']=file.path("~/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/GRCh38.blacklist.bed")
register(MulticoreParam(nrow(dt.samples)))

##
library(ChIPQC)
##
cqo<- ChIPQC(dt.samples, 
             annotation=chip.env$ref, 
             chromosomes=c("21"),
             blacklist=chip.env$blacklist,
             #consensus=T,
             #bCount=T,
             mapQCth=chip.env$MAPQ
)

print(cqo)
QCdba(cqo)
QCmetrics(cqo)
QCsample(cqo)[[1]]
peaks(cqo) # isa GRangeList for a chosen chr
#
ChIPQCreport(cqo, reportName="ChIP-Seq QC report: SLX-16110-16111", reportFolder="~/results/SLX-16110-16111.Homo_sapiens.v1/ChIPQC",facet=F, facetBy=c("Condition"))


####################
library(DiffBind) ##
####################
if(FALSE){
dbo.qc<-ChIPQC::QCdba(cqo) # directly from the ChIPQC object (using default parameters of dba)
dbo.np<-DiffBind::dba(sampleSheet=dt.samples.narrow,
         minOverlap=round(nrow(dt.samples.narrow)*.1), # 10% of samples
         #scoreCol=7, # fold_enrichment
         #scoreCol=8, # -log10(pval) - default for 'narrow'
         #scoreCol=9, # -log10(qval)
         config=data.frame(fragmentSize=chip.env$FR.size,
                           DataType="DBA_DATA_GRANGES",
                           minQCth=chip.env$MAPQ)) # isa 'DBA' - this reads bed (or xls) macs files
}
dbo<-DiffBind::dba(sampleSheet=dt.samples,
         minOverlap=round(nrow(dt.samples)*.1), # 10% of samples
         #scoreCol=7, # -log10(pval) - default for 'macs'
         config=data.frame(fragmentSize=chip.env$FR.size,
                           DataType="DBA_DATA_GRANGES",
                           minQCth=chip.env$MAPQ)) # isa 'DBA' - this reads bed (or xls) macs files
print(dbo)
dba.show(dbo)
dba.show(dbo, attributes=c(DBA_TISSUE,DBA_CONDITION,DBA_REPLICATE,DBA_CALLER))
dba.plotHeatmap(dbo) # same with `plot(dbo)`

sort(names(dbo))
dbo$config
dbo$class[,1:2]
dbo$chrmap
sapply(dbo$peaks, nrow) # NO. of peaks per sample
sapply(dbo$peaks, nrow) == sapply(dbo.np$peaks, nrow) # NO. of peaks per sample
dim(dbo$called) # a matrix of peak called (1) or not (0)
dbo$called[1:3,]
table(apply(dbo$called,1,sum)>=dbo$minOverlap)
table(apply(dbo$called,1,sum)>=dbo$minOverlap)[["TRUE"]]==nrow(dbo$binding)

# binding occupancy score (confidence score) - called at least by minOverlap
# NB, CHR is not real chr names - they are integer assgined for dbo$chrmap
# where the score comes from? - it shows normalised scores with 1 as max (i.e. mat(score)/max(score))
dim(dbo$binding)
dbo$binding[1:3,] 
dbo$peaks[[8]][1:3,] 
data.table(dbo$binding)[,.N,CHR]
# TODO:  <14-10-20, yourname> # check below - they differ
dim(dbo.np$binding)
dbo.np$binding[1:3,]  

## peaks (from macs result files - either narrowPeak or .xls)
data.table(dbo$peaks[[1]])[,.(min(X.log10.pvalue.), max(X.log10.pvalue.))]
data.table(dbo.np$peaks[[1]])[,.(min(V7), max(V7))]

dbo$attributes

dbo$merged[1:3,]
dim(dbo$called)
dbo$merged[1:3,]
dbo$called[1:3,]
nrow(dbo$merged) == nrow(dbo$called)
data.table(dbo$merged)[,.(.N,min(END-START),max(END-START),median(END-START),mean(END-START),sd(END-START))]
hist(dbo$merged[,3]-dbo$merged[,2])
dbo$totalMerged

dba.plotPCA(dbo, attributes=DBA_FACTOR, label=DBA_ID) # based on binding occupancy matrix

##
# mask
##
names(dbo$mask)
dba.show(dbo,dbo$mask$F)
dba.show(dbo,dbo$mask$FV)
split(dt.samples[,-c("bamReads","bamControl","Peaks")], dt.samples$Factor)
split(dt.samples[,-c("bamReads","bamControl","Peaks")], dt.samples$Tissue)

##
## Consensus peaksets
##
print(olap.rate<-dba.overlap(dbo,mode=DBA_OLAP_RATE)) #  No. of peaks by minOverlap (1..)
olap.rate[1] == dbo$totalMerged  
olap.rate[2] == nrow(dbo$binding) 

plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')


dba.overlap(dbo,dbo$masks$FD,mode=DBA_OLAP_RATE) # for FD samples (n=5) only
dba.plotVenn(dbo,dbo$masks$FD)

dbo.consensus<-dba.peakset(dbo, consensus=c(DBA_TISSUE,DBA_CONDITION),minOverlap=2)
names(dbo.consensus$mask)
dba.show(dbo.consensus)
dba.show(dbo.consensus,dbo.consensus$mask$Consensus) # only for the consensus dataset (i.e. peaksets supported by at least 2 samples within factor)
dbo.consensus<-dba(dbo.consensus, mask=dbo.consensus$masks$Consensus, minOverlap=1) # re-construct DBA
dba.show(dbo.consensus)
consensus_peaks <- dba.peakset(dbo.consensus, bRetrieve=TRUE) # isa GRanges
consensus_peaks

####################
## Affiniy scores ##
####################
dbo.cnt<-dba.count(dbo)
dbo.cnt1<-dba.count(dbo, peaks=consensus_peaks) # this applies all samples using the consensus peaks (see above)
dbo.cnt2<-dba.count(dbo, bUseSummarizeOverlaps=TRUE,bRemoveDuplicates=TRUE) #fragmentSize ignored when bUseSummarizeOverlaps is TRUE in dba.count

sort(names(dbo.cnt1))
data.table(dbo.cnt1$merged)[,.(.N,min(END-START),max(END-START),median(END-START),mean(END-START),sd(END-START))]
dbo.cnt1$totalMerged
dbo.cnt2$totalMerged

dim(dbo.cnt1$binding)
dbo.cnt1$binding[1:3,] # binding affinity score - normalised?
dbo.cnt2$binding[1:3,] # binding affinity score

dbo.cnt1$score # ??
dbo.cnt1$SN # FRiP

dba.plotHeatmap(dbo.cnt)
#
dba.plotPCA(dbo.cnt, attributes=DBA_FACTOR, label=DBA_ID) # based on binding matrix
dba.plotPCA(dbo.cnt, attributes=DBA_TISSUE, label=DBA_ID)
dba.plotPCA(dbo.cnt, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotPCA(dbo.cnt, mask=dbo.cnt$mask$M, attributes=DBA_CONDITION, label=DBA_ID)

###########################
## Differential Binding  ##
###########################
dbo.cnt<-dba.contrast(dbo.cnt, categories=DBA_CONDITION)
#dbo.cnt<-dba.contrast(dbo.cnt, categories=DBA_FACTOR)
dbo.cnt<-dba.analyze(dbo.cnt)

dba.normalize(dbo)

# one-go
li.diff<-list() # list of DBA
li.diff[["female"]]<-dba(dbo.cnt,dbo.cnt$masks$F)
li.diff[["male"]]<- dba(dbo.cnt,dbo.cnt$masks$M)
lapply(li.diff, dba.show)
lapply(li.diff, function(i) i$binding[1:3,])
#li.diff<-lapply(li.diff, dba.contrast, categories=DBA_CONDITION)
li.diff<-lapply(li.diff, dba.contrast)
li.diff<-lapply(li.diff, dba.analyze,method=DBA_ALL_METHODS)

lapply(li.diff, print)
lapply(li.diff, dba.show,bContrasts=T)
lapply(li.diff, dba.report,method=DBA_ALL_METHODS,bDB=TRUE,bAll=TRUE)
lapply(li.diff, dba.report,method=DBA_ALL_METHODS, DataType=DBA_DATA_GRANGES)
lapply(li.diff, dba.report,method=DBA_EDGER, DataType=DBA_DATA_GRANGES,th=0.01)
lapply(li.diff, dba.report,method=DBA_DESEQ2, DataType=DBA_DATA_GRANGES,th=0.01)

dba.plotHeatmap(li.diff[["female"]]) # same with `plot(dbo)`
dba.plotHeatmap(li.diff[["male"]]) # same with `plot(dbo)`

dba.plotBox(li.diff[['male']])


##
## Annotation via ChipSeeker
##
libary(ChIPseeker)

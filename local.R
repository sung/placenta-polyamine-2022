library(yaml)
library(BiocParallel)
register(MulticoreParam(32))
library(DiffBind) ##
library(kableExtra)
library(ggplot2)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)

config<-read_yaml('_config.yml') # isa 'list'
project.dir<-file.path('~/results',paste(config$SLX,config$species,config$ver,sep="."))
stopifnot(dir.exists(project.dir))
sample.file<-file.path(project.dir,'Meta','samples.diffbind.csv')
sample.file.np<-file.path(project.dir,'Meta','samples.diffbind.narrowPeak.csv')
group.file<-file.path(project.dir,'Meta','group.diffbind.csv')
stopifnot(file.exists(sample.file))
stopifnot(file.exists(sample.file.np))
stopifnot(file.exists(group.file))

##
## DiffBind Version
##
#RData_dir=ifelse(packageVersion('DiffBind') > 3,file.path('RData','DiffBind3'),file.path('RData','DiffBind3'))
RData_dir=file.path('RData',paste0("DiffBind",packageVersion('DiffBind'))) # 3.0.13 2021-02-16

##
## Samples
##
dt.samples<-fread(sample.file) # to be processed by DiffBind (xls formats by MACS2 for peaks)
dt.samples.narrow<-fread(sample.file.np) # to be processed by DiffBind (narrowPeak formarmat by MACS2 for peaks)
dt.groups<-fread(group.file) # peaks called by pulling (or merging) samples of the same group (FV, FD, MV, MD)

##
## ChIPQC
##
init_ChIPQC<-function(my.samples=dt.samples){
    library(ChIPQC)
    cqo<- ChIPQC(my.samples, 
                annotation=config$ref, 
                chromosomes=c("21"), # this chr for QC
                blacklist=config$blacklist,
                #consensus=T,
                #bCount=T,
                mapQCth=config$MAPQ
    )
    return(cqo)
}

#
# initialize 'dbo' object
init_DiffBind<-function(my.samples=dt.samples, my.type="sample"){
    if(my.type=="sample"){
        my.RData<-file.path(RData_dir,'dbo.RData')
        minOverlap=round(nrow(my.samples)*.1) # 10% of samples (2 in this study)
    }else{
        my.RData<-file.path(RData_dir,'dbo.group.RData')
        minOverlap=1
    }
    if(file.exists(my.RData)){
        message(paste("loading",my.RData))
        attach(my.RData)
    }else{
        #set diffbinf object
        my.dbo<-DiffBind::dba(sampleSheet=my.samples, 
                              minOverlap=minOverlap,
                #scoreCol=7, # -log10(pval) - default for 'macs'
                config= data.frame(fragmentSize=config$FR.size,
                                   AnalysisMethod=config$AnalysisMethod,
                                minQCth=config$MAPQ,
                                doBlacklist=TRUE,
                                doGreylist=FALSE)
        ) # isa 'DBA' - this reads bed (or xls) macs files
        message(paste("saving",my.RData))
        if(my.type=="sample"){
            dbo<-my.dbo
            dbo<-dba.blacklist(dbo) #DiffBind3
            save(dbo,file=my.RData)
        }else{
            dbo.group<-my.dbo
            dbo.group<-dba.blacklist(dbo.group) #DiffBind3
            save(dbo.group,file=my.RData)
        }
    }
    attach(my.RData)
    message("dbo attached")
}


# consens peaks within the factors
init_DiffBind_con<-function(){
    if(!exists('dbo')){init_DiffBind()} # load 'dbo'
    stopifnot(exists('dbo'))

    my.RData<-file.path(RData_dir,'dbo.con.RData')
    if(file.exists(my.RData)){
        message(paste("loading",my.RData))
        attach(my.RData)
    }else{
        # find consensus peak
        dbo.con<-dba.peakset(dbo, consensus=c(DBA_TISSUE,DBA_CONDITION),minOverlap=config$minOverlap) # config$minOverlap=2
        #2020-12-02: more conservative approach
        #dbo.con<-dba.peakset(dbo, consensus=c(DBA_TISSUE,DBA_CONDITION),minOverlap=3) # config$minOverlap=2
        dba.show(dbo.con,dbo.con$mask$Consensus) # only for the consensus dataset (i.e. peaksets supported by at least 2 samples within factor)
        dbo.con<-dba(dbo.con, mask=dbo.con$masks$Consensus, minOverlap=1) # re-construct DBA
        dba.show(dbo.con)

        message(paste("saving",my.RData))
        save(dbo.con,file=my.RData)
        attach(my.RData)
    }
    message("dbo.con attached")
}

#
init_DiffBind_cnt<-function(){
    if(!exists('dbo.con')){init_DiffBind_con()} # load 'dbo.con'
    stopifnot(exists('dbo.con'))

    my.RData<-file.path(RData_dir,'dbo.cnt.RData')
    if(file.exists(my.RData)){
        message(paste("loading",my.RData))
        attach(my.RData)
    }else{
        ## 
        dba.peakset(dbo.con,bRetrieve=T)
        dba.peakset(dbo.con,1:2,bRetrieve=T)
        #
        consensus_peaks <- dba.peakset(dbo.con, bRetrieve=TRUE) # isa GRanges if occupancy matrix (log10p0val based)
        ## Affiniy scores ##
        dbo.cnt.summit<-dba.count(dbo, peaks=consensus_peaks,filter=0,minCount=0,mapQCth=10) # this applies all samples using the consensus peaks (see above). counts based on summit
        dbo.cnt.peak<-dba.count(dbo, peaks=consensus_peaks,summits=FALSE,filter=0,minCount=1,mapQCth=10) # this applies all samples using the consensus peaks (see above)
        dbo.cnt<-dbo.cnt.peak

        dba.peakset(dbo.cnt, bRetrieve=T) # score=DBA_SCORE_TMM_MINUS_FULL (DiffBind2); DBA_SCORE_NORMALIZED (DiffBind3)
        dba.peakset( dba.count(dbo.cnt, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE)
        dba.peakset( dba.count(dbo.cnt, peaks=NULL, score=DBA_SCORE_READS_FOLD, bLog=T), bRetrieve=TRUE)

        foo<-dba.peakset( dba.count(dbo.cnt, peaks=NULL, score=DBA_SCORE_READS_MINUS), bRetrieve=TRUE)
        foo[seqnames(foo)=="8" & start(foo)==95099279 & end(foo)==95099683,grepl("F*D",colnames(mcols(foo)))]
        foo[seqnames(foo)=="8" & start(foo)==95099279 & end(foo)==95099683,grepl("F*V",colnames(mcols(foo)))]

        message(paste("saving",my.RData))
        save(dbo.cnt,file=my.RData)
        attach(my.RData)
    }
    message("dbo.cnt attached")
}

init_DiffBind_diff<-function(){
    message('making li.diff object...')
    if(!exists('dbo.cnt')){init_DiffBind_cnt()} # load 'dbo.cnt'
    stopifnot(exists('dbo.cnt'))

    my.RData<-file.path(RData_dir,'li.diff.RData')
    if(file.exists(my.RData)){
        message(paste("loading",my.RData))
        attach(my.RData)
    }else{
        li.diff<-list() # list of DBA

        li.diff[['bg']][["female"]]<-dba(dbo.cnt,dbo.cnt$masks$F) # isa DBA
        li.diff[['bg']][["male"]]<- dba(dbo.cnt,dbo.cnt$masks$M) # isa DBA

        ##li.diff<-lapply(li.diff, dba.contrast,design="~Replicate+Factor")
        #li.diff<-lapply(li.diff, dba.contrast)
        #li.diff<-lapply(li.diff, dba.analyze,method=DBA_ALL_METHODS)

        li.diff[['bg']]<-lapply(li.diff[['bg']], dba.normalize, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,library=DBA_LIBSIZE_FULL,background=TRUE) # 'background' normalisation via dba.normalize from DiffBind3
        li.diff[['bg']]<-lapply(li.diff[['bg']], dba.contrast, design="~Replicate+Factor")
        li.diff[['bg']]<-lapply(li.diff[['bg']], dba.analyze, method=DBA_ALL_METHODS)

        lapply(li.diff[['bg']], dba.show, bContrasts=TRUE)
        lapply(li.diff[['bg']][c('male')], dba.report, method=DBA_ALL_METHODS,bDB=T,bGain=T,bLoss=T)
        lapply(li.diff[['bg']], dba.normalize, bRetrieve=T,method=DBA_ALL_METHODS)[["male"]]

        message(paste("saving",my.RData))
        save(li.diff,file=my.RData)
        attach(my.RData)
    }
    message("li.diff attached")
}

## to be used by ChIPseeker
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
GenomeInfoDb::seqlevelsStyle(txdb)="NCBI" # based on 'X' not 'chrX'

# GRange summit of the consensus_peaks
# to be used by deepTools to plot bas-base depth coverage of peak by FD, FV, MD, MV
if(FALSE){
    consensus_peaks
    foo<-GRanges(seqnames=seqnames(consensus_peaks),
            IRanges(
            start=start(consensus_peaks)+round(width(consensus_peaks)/2),
            end=start(consensus_peaks)+round(width(consensus_peaks)/2)+1,
            ))
    export.bed(foo,file.path(project.dir,'DiffBind/consensus_peaks.bed'))
}

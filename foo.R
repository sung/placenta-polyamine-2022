
gl.dbr<-lapply(li.diff, dba.report, method=DBA_DESEQ2, DataType=DBA_DATA_GRANGES)
gl.dbr<-lapply(li.diff, dba.report, method=DBA_EDGER, DataType=DBA_DATA_GRANGES)

li.peakAnno<-lapply(gl.dbr['male'], ChIPseeker::annotatePeak, tssRegion=c(-2e3, 2e3), TxDb=txdb, annoDb="org.Hs.eg.db", verbose=F) # list of 'csAnno'

dt.foo<-data.table(data.frame(li.peakAnno[['female']]@anno))
dt.foo[Fold>0] # FV > FD
dt.foo[Fold<0] # FV < FD

dt.bar<-data.table(data.frame(li.peakAnno[['male']]@anno))
dt.bar[FDR<1e-6,.N,Fold>0] # MV > MD
dt.bar[FDR<1e-6]

dt.bar[order(-Fold)]


##
## Count
##
names(dbo.cnt)
dbo.cnt$peaks[[1]][1:3,]
sapply(dbo.cnt$peaks, function(i) sum(i$Reads))
sapply(dbo.cnt$peaks, function(i) sum(i$cReads))

dbo.cnt$binding[1:3,] # affinity matrix (use below)

dba.peakset(dbo.cnt,  bRetrieve=TRUE) # # default: DBA_SCORE_NORMALIZED

dba.peakset(li.diff[['female']], bRetrieve=T) # same above (but for female)
li.diff[['female']]$binding[1:3,]
li.diff[['female']]$peaks[[1]][1:3,]
names(li.diff[['female']]$contrast[[1]])
li.diff[['female']]$contrast[[1]]$group1
names(li.diff[['female']])

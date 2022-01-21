library(ChIPseqSpikeInFree) # tested on R 3.6.1

project.dir<-file.path('~/results/SLX-16110-16111.Homo_sapiens.v1')
dt.samples<-fread(file.path(project.dir,'Meta/samples.diffbind.csv'))
dt.samples[1:3]
male.meta<-dt.samples[Tissue=="M",.(ID=bamReads,GROUP=Factor,ANTIBODY="H3K27ac")]
female.meta<-dt.samples[Tissue=="F",.(ID=bamReads,GROUP=Factor,ANTIBODY="H3K27ac")]
male.bams<-male.meta$ID
female.bams<-female.meta$ID

male.meta$ID<-basename(male.meta$ID)
female.meta$ID<-basename(female.meta$ID)

male.meta.file=file.path(project.dir,'ChIPseqSpikeInFree/Meta/male.meta.txt')
fwrite(male.meta,file=male.meta.file,sep="\t")
female.meta.file=file.path(project.dir,'ChIPseqSpikeInFree/Meta/female.meta.txt')
fwrite(female.meta,female.meta.file,sep="\t")

ChIPseqSpikeInFree(bamFiles = male.bams, chromFile = "hg38", metaFile = male.meta.file, prefix = "Male", ncores=16)
ChIPseqSpikeInFree(bamFiles = female.bams, chromFile = "hg38", metaFile = female.meta.file, prefix = "Female", ncores=16)

male.bams<-male.meta$ID
bamFiles<-male.bams
bamlist <- BamFileList(bamFiles, yieldSize = 5e4)
names(bamlist)

male.bams<-male.meta$ID 
female.bams<-female.meta$ID 



data="Male_rawCounts.txt"
metaFile=male.meta.file
by=0.05
binSize=1000
ncores=16
parsedDF <- ParseReadCounts(data="Male_rawCounts.txt",metaFile=male.meta.file, prefix="Male",ncores=16)

intersect(c("A","B","C"),c("A","B"))

    metaFile <- ReadMeta(metaFile)
metaFile$ID
colnames(data)
colnames(data)
$ID
, colnames(data))
  kept <- intersect(metaFile$ID, colnames(data))

rownames(metaFile[1:2,])

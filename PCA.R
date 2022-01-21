source("libs/theme_publish.R")
# based on 02-calling.Rmd

dbo.cnt %>% dba.show
dba.peakset(dbo.cnt,bRetrieve=T)

# a legacy code from DiffBind
p<-dba.plotPCA(dbo.cnt, attributes=DBA_FACTOR, label=DBA_ID) # based on binding affinity matrix # isa 'trellis'
p %>% class # isa 'trellis'

# in-house code to reproduce the PCA plot via ggplot 
mat.cnt<-dbo.cnt$binding[,dba.show(dbo.cnt)$ID]
mat.cnt %>% dim
mat.cnt[1:10,]
bad.peak <- rowSums(mat.cnt) ==0
mat.cnt<-mat.cnt[!bad.peak,]

# set min value
mat.cnt[mat.cnt<=1] <- 1 # min=1
mat.cnt[1:10,]

# log2 transform
mat.cnt <- log2(mat.cnt)
mat.cnt[1:10,]

## PCA via prcomp
p.pc<-(mat.cnt %>% t %>% prcomp)
p.pc$x[,c("PC1","PC2")]
#plot(p.pc$x[,c("PC1","PC2")])
dt.foo<-merge(
    data.table(ID=rownames(p.pc$x),p.pc$x[,c("PC1","PC2")]),
    dba.show(dbo.cnt)
)

my.sex.col<-ggsci::pal_aaas("default")(2)
names(my.sex.col)<-c("M","F")
p1<-ggplot(dt.foo, aes(PC1,PC2,col=Tissue,shape=Condition)) +
    geom_point(size=5,stroke=3) +
    scale_shape_manual(values = c(`D`=21,`V`=19)) +
    scale_color_manual(values=my.sex.col) +
    theme_Publication() + theme(legend.position="")

file.name=file.path("docs5/PCA.norm.cnt")
pdf(file=paste(file.name,'pdf', sep ='.'), width=9, height=9,title="PCA plot based on the normalised peak counts")
p1
dev.off()


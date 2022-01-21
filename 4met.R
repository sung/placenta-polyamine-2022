library(data.table)

dt.foo<-fread("~/Downloads/4met_cor.csv")
dt.cast.foo<-dcast.data.table(dt.foo, transcript+hgnc_symbol+description~metabolite,value.var=c("rho","p","CI_low","CI_high","S"))

fwrite(dt.cast.foo, file="~/Downloads/4met_cor_cast.csv")

dt.foo[metabolite %in% dt.foo[p<0.05,.N,metabolite][N>30]$metabolite & p<0.05,.N,transcript]

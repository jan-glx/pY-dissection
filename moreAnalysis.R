library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(data.table)


df = fread("../Data/valid_interactions_grossmann.tsv")
df[,kinase_dep:=strsplit(kinase_dep,", ")]
kin_hist=df[,.(kinase=do.call(c,kinase_dep)),by=Interaction.identifier.s.][,.N,keyby=kinase]


df[,kinase_dep_factor:=lapply(kinase_dep,factor,levels=kin_hist[,kinase])]





dt=df[,.(kinase=do.call(c,kinase_dep)),by=Interaction.identifier.s.]
dt[,kinase:=as.factor(kinase)]
all_kinases=levels(dt[,kinase])
n_kinases=length(all_kinases)





a=dt[,.(ii=rep(kinase,each=.N),jj=rep(kinase,.N)),by=Interaction.identifier.s.][,.N,keyby=.(ii,jj)]
jac=matrix(0,ncol=n_kinases,nrow = n_kinases,dimnames = list(all_kinases,all_kinases))
jac[cbind(a[,ii],a[,jj])]=a[,N]
jac=jac+1/n_kinases^2
jac=jac/(matrix(kin_hist[,N]+1/n_kinases,n_kinases,n_kinases)-jac+ matrix(kin_hist[,N]+1/n_kinases,n_kinases,n_kinases,byrow=T))
format(round(jac,digits = 3),scientific=F)
fit <- hclust(as.dist(-log(jac)))
plot(fit)


retested=setkey(df,Interaction.identifier.s.)[dt[!(kinase %in% c(the_nine,"")),
                                                 .N,
                                                 by=Interaction.identifier.s.
                                                 ][order(N)]]

retested= retested[,.(kinase_dep,N,canonic_ID_bait,canonic_ID_prey)]
setkey(uniprot_data,ID)
setkey(retested,canonic_ID_bait)
retested=setkey(setkey(retested[uniprot_data[,.(ID,name_bait = Gene.names...primary..)],
         nomatch=0
         ],canonic_ID_prey)[uniprot_data[,.(ID,name_prey = Gene.names...primary..)],
           nomatch=0
           ],name_bait,name_prey)
retested
#  [,.(name_bait,name_prey,kinase_dep)]  # [28,kinase_dep]

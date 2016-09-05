source("misc.R")
#loadings all sites
all_sites <- freadlist("../Data/psp_all_sites_all_info.tsv")[!str_detect(SUB_ACC_ID, "-")]
setnames(all_sites, "SUB_ACC_ID", "ac_uniprot")
setkey(all_sites, "ac_uniprot", "res")


#jointing allsites into surface accesiblity
all_sites <- nacc[all_sites, on=c("ac_uniprot","res")]

ggplot(all_sites[!is.na(acc_s)],aes(x=acc_s,color=myPSP))+stat_ecdf()
wilcox.test(all_sites[!is.na(acc_s)&is.na(myPSP),acc_s],all_sites[!is.na(acc_s)&myPSP,acc_s])


# loading pfam
pfam <- fread("zcat ../Data/pfam_instances_uniprot.tsv.gz")

dt  <- all_sites[,.(ac_uniprot, res,myPSP)][pfam,,allow.cartesian=TRUE,on="ac_uniprot",nomatch=0][res<=end_seq&start_seq<=res][, c("start_seq", "end_seq", "e_value", "name") := NULL]

dt4<- dt[all_sites[,.(ac_uniprot, res, myPSP)], on= c("ac_uniprot", "res", "myPSP")][is.na(ac_pfam),ac_pfam:="none"]



# analyisis of pfam domain information
res <- dcast(dt4[,.(k=.N), keyby=.(ac_pfam,myPSP)
                 ][all_sites[,.(n=.N),keyby=myPSP],on="myPSP"
                   ], ac_pfam ~ myPSP,value.var=c("k","n"), fill=0)

res[,c("p.value","lower","upper","or.est"
):=(
  with(fisher.test(matrix(c(k_TRUE, k_NA, n_TRUE, n_NA),ncol=2)),
       list(p.value=p.value,
            lower=conf.int[1],
            upper=conf.int[2],
            or.est=estimate))
), by=ac_pfam]
setorder(res,p.value)
res[,p.value.bonferoni:=p.value*nrow(res)][,p.value.benjamini.hochberg:=p.value*nrow(res)/.I]

ress <- res[rev(0<cumsum(rev(p.value.benjamini.hochberg<0.05)))]

setorder(ress,-lower)
head(ress, 10)
cat(head(ress, 10)[,ac_pfam])
# Immunoreceptor tyrosine-based activation motif, Neuraxin and MAP1B repeat, Repeat in HS1/Cortactin
# Phosphoprotein associated with glycosphingolipid-enriched, Ephrin type-A receptor 2 transmembrane domain,
# Histone-like transcription factor (CBF/NF-Y) and archaeal histone, Tubulin C-terminal domain
# GAGE protein, Tubulin/FtsZ family, GTPase domain
setorder(ress,upper)
head(ress, 10)
cat(head(ress, 10)[,ac_pfam])
# Rhodopsin-like receptors, Rhodopsin-like receptors, Immunoglobulin domain, 
# Ubiquitin carboxyl-terminal hydrolase, Immunoglobulin I-set domain, Ankyrin repeat 





#library(ggplot2)
#all_sites[,myPSP2:=if (isTRUE(any(myPSP))) sapply(myPSP, isTRUE) else myPSP, by=ac_uniprot]
#ggplot(all_sites,aes(x=IUPred_prob,color=myPSP2))+stat_ecdf()+geom_vline(x=0.5,linetype="dotted")
cat("IUPRED and phosphorylation")
fisher.test(matrix(all_sites[,.N,keyby=.(IUPred_prob>0.5,myPSP)][,N],ncol=2))


cat("IUPRED and phosphorylation")
fisher.test(matrix(all_sites[!is.na(acc_s),.N,keyby=.(acc_s>0.5,myPSP)][,N],ncol=2))












if (FALSE) {
  # Mutual information test
  library(mpmi)
  all_sites[,"psp":=sapply(myPSP, isTRUE)]
  tmp <- all_sites[!is.na(acc_s),]
  tmp <- tmp[ 1==rbinom(nrow(tmp),1,prob=0.1)]
  mminjk(tmp[,.(acc_s,IUPred_prob,seq_length,ANCHOR_score)],tmp[,.(psp)])
  #cminjk(tmp)
}

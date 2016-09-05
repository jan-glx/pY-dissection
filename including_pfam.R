source("misc.R")
#loadings all sites
all_sites <- freadlist("../Data/psp_all_sites_all_info.tsv")[!str_detect(SUB_ACC_ID, "-")]
setnames(all_sites, "SUB_ACC_ID", "ac_uniprot")
setkey(all_sites, "ac_uniprot", "res")


#loading surface accesibility
nacc <- fread("zcat ../Data/naccess_human_uniprot.tsv.gz", select=c("ac_uniprot","pos","e_value","acc_s"))
setnames(nacc, "pos", "res")
setkey(nacc,ac_uniprot,res,e_value)[,e_value:=NULL]
setkey(nacc,ac_uniprot,res)
nacc <- unique(nacc)
#jointing allsites into surface accesiblity
all_sites <- nacc[all_sites, on=c("ac_uniprot","res")]

# loading pfam
pfam <- fread("zcat ../Data/pfam_instances_uniprot.tsv.gz")

dt  <- all_sites[,.(ac_uniprot, res,myPSP)][pfam,,allow.cartesian=TRUE,on="ac_uniprot",nomatch=0][res<=end_seq&start_seq<=res][, c("start_seq", "end_seq", "e_value", "name") := NULL]
most_common_100_pfams <- tail(setkey(dt[(myPSP),.N,by=ac_pfam],N),100)
dt[!most_common_100_pfams,ac_pfam:="PF_other",on=("ac_pfam")]
dt=unique(dt, by=c("ac_uniprot", "res", "ac_pfam"))
dt=setkey(dcast(dt[,tmp:=TRUE],ac_uniprot+res+myPSP~ac_pfam,value.var="tmp",fill=FALSE),ac_uniprot, res)[]
dt=dt[all_sites, on= c("ac_uniprot", "res", "myPSP")]
setkey(dt,ac_uniprot,res)

phospho_elm_sub <- fread("../Data/phospho.ELM.pYReader.tsv")
setnames(phospho_elm_sub, c("Accession","Position"), c("ac_uniprot","res"))
phospho_elm_sub=dcast(phospho_elm_sub[,helper:=TRUE],ac_uniprot+res ~ pYReader,value.var="helper",fill=FALSE)
setkey(phospho_elm_sub,ac_uniprot,res)

dt <- phospho_elm_sub[dt]

peptospot <- fread("../Data/sh2_zscores.csv")
peptospot[,neighbourhood:=paste0("_",sequence,"_")][,sequence:=NULL]

# replace nan by FALSE or 0 where appropiate
kincols <- colnames(dt)
kincols <- kincols[(which(kincols=="res")+1) : (which(kincols=="LT_LIT")-1)]
for(col in c("myPSP",kincols,most_common_100_pfams[,ac_pfam],"PF_other",colnames(phospho_elm_sub)[-c(1,2)])){
  dt[is.na(get(col)),(col):=FALSE][]
}
for(col in c("LT_LIT","MS_LIT","MS_CST")){
  dt[is.na(get(col)),(col):=0][]
}




fwritelist(dt,"../Data/psp_all_sites_all_info2.tsv")


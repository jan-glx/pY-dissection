library(readr)
library(stringr)
library(data.table)

seqs = data.table(str_match_all(read_file("../Data/Phosphosite_seq.txt"),
                                ">[^\\|]*\\|human\\|([^\n]+)\n([^>]+)")[[1]][,2:3])
setnames(seqs,c("ACC","sequence"))
seqs[,sequence:=str_replace_all(sequence,"\n","")]
sites=seqs[,.(res=str_locate_all(sequence,'Y')[[1]][,1],
            seq_length=str_length(sequence)),
            by=.(ACC,sequence)][!is.na(res)]

sites[,neighbourhood:=paste0(str_sub("_______",0,pmax(0,8-res)),
                             str_sub(sequence,pmax(1,res-7),pmin(seq_length,res+7)),
                             str_sub("_______",0,pmax(0,7-(seq_length-res))))][,sequence:=NULL]

Rcpp::sourceCpp('anchor.cpp')
setkey(seqs,ACC)
dt = seqs[unique(sites[,ACC]),rbindlist(lapply(sequence,function(seq){fread(anchor(seq))[`One letter code`=="Y"]})),by=ACC]
setnames(dt,c("ACC","res","res_kind","ANCHOR_prob","ANCHOR_binary","IUPred_prob","ANCHOR_score","S","Eint","Egain"))



setkey(sites,ACC,res)
setkey(dt,ACC,res)
sites=dt[sites,mult="first"][,res_kind:=NULL]

fwritelist(sites,"../Data/pps_all_sites_with_iupred.tsv")

# add pspPos info to psp  -----
sites <- freadlist("../Data/pps_all_sites_with_iupred.tsv")

pos_sites = freadlist("../Data/all_sites_pspPos.tsv")[,':='(neighbourhood=NULL,GENE=NULL)]
setkey(pos_sites,SUB_ACC_ID,res)
a=pos_sites[sites]
fwritelist(a,"../Data/psp_all_sites_all_info.tsv")





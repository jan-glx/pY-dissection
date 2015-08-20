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


seqs
sites

Rcpp::sourceCpp('anchor.cpp')
dt = seqs[1:50,rbindlist(lapply(sequence,function(seq){fread(anchor(seq))[`One letter code`=="Y"]})),by=ACC]
dt
setnames(fread(anchor("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")),
              c("res","res_kind","ANCHOR_prob","ANCHOR_binary","IUPred_prob","ANCHOR_score",
                "S","Eint","Egain"))












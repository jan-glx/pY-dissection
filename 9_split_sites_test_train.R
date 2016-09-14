source('misc.R')

all=freadlist("../Data/psp_all_sites_all_info.tsv")
setkey(all,SUB_ACC_ID)
test_ids = fread("../Data/uniprotdata_4all_hits.tsv",sep="\t",
                     select=c("canonic"))
all=all[!str_detect(SUB_ACC_ID,"-")]
fwritelist(all[!test_ids],"../Data/train_sites_for_NN.tsv")
fwritelist(all[test_ids],"../Data/test_sites_for_NN.tsv")

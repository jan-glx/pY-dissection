library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

source("misc.R")

grossmann = fread("../Data/valid_interactions_grossmann.tsv")

# finding sites and uniprot infor for them ----------
uniprot_data = fread("../Data/uniprotdata_4all_hits.tsv",sep="\t",
                     select=c("canonic","Entry.name", "ID","Gene.names...primary..", "Sequence","Phosphotyrosins"))
setnames(uniprot_data,"Gene.names...primary..","gene_name")

uniprot_data[,Phosphotyrosins:=lapply(str_split(Phosphotyrosins,","),as.integer)]
uniprot_data[is.na(Phosphotyrosins),Phosphotyrosins:=lapply(Phosphotyrosins,function(x){integer()})]

sites=uniprot_data[,.(res=str_locate_all(Sequence,'Y')[[1]][,1],
                      seq_length=str_length(Sequence)),
                   by=.(canonic,ID,Entry.name, gene_name,Sequence)][!is.na(res)]

sites[,neighbourhood:=paste0(str_sub("_______",0,pmax(0,8-res)),
                             str_sub(Sequence,pmax(1,res-7),pmin(seq_length,res+7)),
                             str_sub("_______",0,pmax(0,7-(seq_length-res))))][,Sequence:=NULL]


setkey(uniprot_data,ID)
setkey(sites,ID)
setkey(uniprot_data,canonic)
setkey(sites,canonic)

sites[uniprot_data,phosphorylated_by_uniprot:= mapply('%in%',res, Phosphotyrosins) ]

# adding evidence from phophositeplus -------------

psp = fread("../Data/Phosphorylation_site_dataset.tsv",sep = "\t",sep2="; ",#colClasses = c(`CST_CAT#`="character"),
            select=c("ACC_ID","GENE","MOD_RSD","SITE_+/-7_AA","LT_LIT","MS_LIT","MS_CST"))
setnames(psp,"SITE_+/-7_AA","neighbourhood")

psp[,c("res_kind","res"):=data.table(str_match(MOD_RSD,"([YST])(\\d+)-p")[,-1])][,MOD_RSD:=NULL]
psp <- psp[res_kind=="Y"][,':='(res_kind=NULL,res=as.integer(res))][]
psp[is.na(LT_LIT),LT_LIT:=0][is.na(MS_LIT),MS_LIT:=0][is.na(MS_CST),MS_CST:=0]
printSetDiffSizes(psp[,ACC_ID],uniprot_data[,canonic])


setkey(sites,canonic,res)
setkey(psp,ACC_ID,res)


if(sites[psp,any(neighbourhood!=toupper(i.neighbourhood)),nomatch=0])
  dfg
if(sites[psp,any(GENE!=gene_name),nomatch=0])
  dfg


sites[,':='(LT_LIT=0L,MS_LIT=0L,MS_CST=0L)]
sites[psp,':='(LT_LIT=i.LT_LIT,MS_LIT=i.MS_LIT,MS_CST=i.MS_CST),nomatch=0]




# adding kinase specific data from phophositeplus ---------------

pspSub = fread("../Data/Kinase_Substrate_Dataset.tsv",sep = "\t", 
               select=c("KINASE","KIN_ACC_ID","GENE","KIN_ORGANISM","SUBSTRATE","SUB_ACC_ID","SUB_GENE" ,   
                     "SUB_ORGANISM","SUB_MOD_RSD","SITE_+/-7_AA","IN_VIVO_RXN","IN_VITRO_RXN"))
setnames(pspSub,"SITE_+/-7_AA","neighbourhood")
pspSub = pspSub[KIN_ORGANISM=="human"][,':='(KIN_ORGANISM=NULL,
                                             IN_VIVO_RXN=IN_VIVO_RXN=="X",
                                             IN_VITRO_RXN=IN_VITRO_RXN=="X")][]
pspSub[,c("res_kind","res"):=data.table(str_match(SUB_MOD_RSD,"^([YST])(\\d+)$")[,-1])][,SUB_MOD_RSD:=NULL]
pspSub <- pspSub[res_kind=="Y"][,':='(res_kind=NULL,res=as.integer(res))][]


the_nine=c('FYN (P06241)','ABL2 (P42684)','TNK1 (Q13470)','FRK (P42685)','FES (P07332)','PTK2 (Q05397)',
           'SYK (P43405)','BMX (P51813)','JAK2 (O60674)')
the_nine=setnames(data.table(str_match(the_nine,"(.+) \\((.+)\\)")),c("full_string","name","ACC"))

setkey(pspSub,KIN_ACC_ID)
setkey(the_nine,ACC)

pspSub = pspSub[the_nine,nomatch=0]

setkey(pspSub,SUB_ACC_ID,res)
setkey(sites,canonic,res)

sites_pspSub=pspSub[,.(kinases=list(GENE)),by=.(SUB_ACC_ID,res)]
sites[,kinases_psp:=rep(list(character()),.N)][sites_pspSub,kinases_psp:=kinases,nomatch=0]








# adding mechismo output -------------

filename="enFcexLbnj.site_table.tsv"

mechismo_out = fread(file.path(folder,filename),
                     select=c("primary_id_a1", "primary_id_b1", "pos_a1", "mut_a1", 
                              "user input", "mismatch", "conf", "ie","name_b1","idcode","iupred"))[primary_id_a1!=""]

mechismo_out = mechismo_out[name_b1=="" | (name_b1=="[PROT]" ) | (primary_id_b1 %in% unique(uniprot_data$canonic))]
setkey(mechismo_out,primary_id_a1,primary_id_b1,pos_a1,name_b1)
mechismo_out = mechismo_out[mechismo_out[,.(pdb_ids=list(unique(idcode)),
                                            mean_iupred=mean(iupred)),
                                         by=.(primary_id_a1,primary_id_b1,pos_a1,name_b1)
                                         ],
                            mult="first"][,c("idcode","iupred"):=NULL]


dropped = mechismo_out[mismatch!=0,.N,by=primary_id_a1]$primary_id_a1
if (length(dropped) != 0) {
  cat("Waring: dropped Entries with Interactor(s):\n ", paste0(dropped, collapse = ", "),
      "\nbecause mechismo indicates that it used a different sequence as basis. (older uniprot sequence version)\n\n"
  )
  mechismo_out = mechismo_out[!(primary_id_a1 %in% dropped)]
}
mechismo_out[,mismatch:=NULL]

mechismo_out = mechismo_out[mut_a1 == "Yp"][,mut_a1:=NULL]



mechismo_sites = mechismo_out[,.(ie=list(ie[!is.na(ie)]),
                iu_pred=mean(mean_iupred),
                n_pdbs=length(unique(unlist(pdb_ids)))
                ),by=.(primary_id_a1,pos_a1)]
setkey(mechismo_sites,primary_id_a1, pos_a1)


sites=mechismo_sites[sites,nomatch=0]





# unstructured regions are more likeli to have phophorilated Tyrosins (ods ratio: 3.15 (2.27..4.31 95%cI) )
fisher.test(table(sites[,.(iu_pred,phosphorylated_by_uniprot)]))
sites[,mean(phosphorylated_by_uniprot),by=.(iu_pred)]




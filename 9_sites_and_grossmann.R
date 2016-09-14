library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

source("misc.R")
sites= freadlist(file.path("../Data/sites.tsv"),colClasses=c("kinases_psp"="character"))
grossmann = freadlist("../Data/valid_interactions_grossmann.tsv",sep2=", ")

the_nine=c('FYN (P06241)','ABL2 (P42684)','TNK1 (Q13470)','FRK (P42685)','FES (P07332)','PTK2 (Q05397)',
           'SYK (P43405)','BMX (P51813)','JAK2 (O60674)')
the_nine=setnames(data.table(str_match(the_nine,"(.+) \\((.+)\\)")),c("full_string","name","ACC"))

#grossmann[,the_nine[,name]:=data.table(t(sapply(kinase_dep,function(x)the_nine[,full_string]%in%x)))]

grossmann_kinases = grossmann[
  grossmann[,.("kinase"=unlist(kinase_dep)),by=.(Interaction.identifier.s.)][kinase!=""],
  .(confidence,canonic_ID_bait,canonic_ID_prey,kinase),
  on=c("Interaction.identifier.s."="Interaction.identifier.s.")][
    the_nine[,.(full_string,kinase_name=name)],
             on=c("kinase"="full_string")][,kinase:=NULL]
    
sites_kinases = sites[,.(kinase=unlist(kinases_psp)),by=.(primary_id_a1, pos_a1)]




sites_kinases[grossmann_kinases,
              on=c("primary_id_a1"="canonic_ID_bait","kinase"="kinase_name"),
              nomatch=0]
sites_kinases[grossmann_kinases,
              on=c("primary_id_a1"="canonic_ID_prey","kinase"="kinase_name"),
              nomatch=0]

res=rbindlist(list(
grossmann_kinases[sites_kinases,
                  on=c("canonic_ID_bait"="primary_id_a1","kinase_name"="kinase"),
                  nomatch=0][,about:="bait"],
grossmann_kinases[sites_kinases,
                  on=c("canonic_ID_prey"="primary_id_a1","kinase_name"="kinase"),
                  nomatch=0][,about:="prey"]))
setkey(res,canonic_ID_bait,canonic_ID_prey,kinase_name)

fwritelist(res,"../Data/knownPhosphoSitesIninteractions.tsv")


grossmann[sites[,.(primary_id_a1,kinases_psp_prey=kinases_psp)],
          on=c("canonic_ID_prey"="primary_id_a1")
          ][sites[,.(primary_id_a1,kinases_psp_bait=kinases_psp)],
            on=c("canonic_ID_bait"="primary_id_a1")]


proteins_with_info <- unique(sites[sapply(kinases_psp,length)!=0,primary_id_a1])
proteins_with_info <- unique(sites[(phosphorylated_by_uniprot),primary_id_a1])
proteins_with_info <- unique(sites[(MS_CST+MS_LIT+LT_LIT)>0,primary_id_a1])
grossmann[canonic_ID_bait %in% proteins_with_info & canonic_ID_prey %in% proteins_with_info,.N]


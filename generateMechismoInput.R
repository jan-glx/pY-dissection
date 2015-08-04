library(stringr)

rm(list = ls())
df = read.table(
  "../Data/all_intact_data_from_grossmanEtAl.tsv", sep = "\t",header = T,stringsAsFactors = F
)
uniprot_data = read.table(
  "../Data/uniprotdata_4all_hits.tsv", sep = "\t",header = T,colClasses = "character"
)
load(file="../Data/canonical2isoform_idx.RData")


a=uniprot_data[150,]
b= canonical2isoform_idx[[150]]


Y_pos_list = lapply(str_locate_all(uniprot_data$Sequence,"Y"),function(x){x[,1]})
mod_pos_list = lapply(str_match_all(uniprot_data$mechismo_dif_string, "[A-Z](\\d+)[A-Z]+"),function(x){as.numeric(x[,2])})

unmod_Y_pos_list = mapply(setdiff , Y_pos_list, mod_pos_list)

inserted_Y_pos_list = lapply(str_locate_all(uniprot_data$mechismo_dif_string, "Y(?!\\d)"),function(x){x[,1]})


generateAllMechismoStrings = function(ID, mechismo_dif_string, inserted_Y_pos,unmod_Y_pos){
  mechismo_strings= paste(ID,mechismo_dif_string)
  if(length(inserted_Y_pos)>0)
    mechismo_strings = c(mechismo_strings, paste0(ID, " ", 
                              str_sub(mechismo_dif_string,1,inserted_Y_pos),
                              "p",str_sub(mechismo_dif_string,inserted_Y_pos+1,-1)))
  if(length(unmod_Y_pos)>0)
    mechismo_strings = c(mechismo_strings, paste0(ID, " ",mechismo_dif_string," Y",unmod_Y_pos,"Yp"))
  mechismo_strings
}
generateAllMechismoStrings(uniprot_data$canonic[328],uniprot_data$mechismo_dif_string[328],
                           inserted_Y_pos_list[[328]],unmod_Y_pos_list[[328]])
generateAllMechismoStrings(uniprot_data$canonic[1],uniprot_data$mechismo_dif_string[1],
                           inserted_Y_pos_list[[1]],unmod_Y_pos_list[[1]])

unmod_Y_pos_list = mapply(generateAllMechismoStrings, 
                          uniprot_data$canonic,
                          uniprot_data$mechismo_dif_string, 
                          inserted_Y_pos_list,
                          unmod_Y_pos_list)
cat(unlist(unmod_Y_pos_list),file="../Data/mechismoinput.txt",sep="\n")



library(stringr)
library(data.table)

generateMechismoInput <- function(uniprot_data){
  
  Y_pos_list = lapply(str_locate_all(uniprot_data$Sequence,"Y"),function(x){x[,1]})
  mod_pos_list = lapply(str_match_all(uniprot_data$mechismo_dif_string, "[A-Z](\\d+)[A-Z]+"),function(x){as.numeric(x[,2])})
  
  unmod_Y_pos_list = mapply(setdiff , Y_pos_list, mod_pos_list)
  
  generateAllMechismoStringsForInsertion <-
    function(ID,mechismo_dif_string, inserted_Y_pos) {
      if (length(inserted_Y_pos) > 0)
        paste0(
          ID, "/", str_sub(mechismo_dif_string,1,inserted_Y_pos),
          "p",str_sub(mechismo_dif_string,inserted_Y_pos + 1,-1)
        )
    }
  
  generateAllMechismoStrings = function(ID, mechismo_dif_string, unmod_Y_pos){
    if(mechismo_dif_string!=""){
      mechismo_dif_strings = str_split(mechismo_dif_string," ")[[1]]
      mechismo_strings= paste0(ID,'/',mechismo_dif_strings)
      inserted_Y_pos_lists = lapply(str_locate_all(mechismo_dif_strings, "Y(?!\\d)"),function(x){x[,1]})
      mechismo_strings = c(mechismo_strings, unlist(mapply(generateAllMechismoStringsForInsertion,
                                                         ID,mechismo_dif_strings,inserted_Y_pos_lists)))
    }else{
      mechismo_strings= NULL
    }
    if(length(unmod_Y_pos)>0)
      mechismo_strings = c(mechismo_strings, paste0(ID, "/Y",unmod_Y_pos,"Yp"))
    mechismo_strings
  }
  
  mechismo_input_line_list = mapply(generateAllMechismoStrings, 
                            uniprot_data$canonic,
                            uniprot_data$mechismo_dif_string, 
                            unmod_Y_pos_list)
  cat(unlist(mechismo_input_line_list),file="../Data/mechismoinput.txt",sep="\n")
}
uniprot_data = fread("../Data/uniprotdata_4all_hits.tsv")
generateMechismoInput(uniprot_data)

setkey(uniprot_data,ID)
inputs = unique(uniprot_data)[,.(seq_length = str_length(Sequence),
                                 input = paste0(">", ID, "\n", Sequence))
                              ][seq_length>=15&seq_length<=4000][,cum_seq_length := cumsum(seq_length)]
cat(inputs[cum_seq_length<200000][1:.N<=2000,input],
    sep = "\n",
    file= paste0("../Data/NetPhos2_input_",1,".fasta"))


setkey(uniprot_data,ID)
inputs = unique(uniprot_data)[,.(seq_length = str_length(Sequence),
                                 input = paste0(">", ID, "\n", Sequence))
                              ][seq_length>=5][,cum_seq_length := cumsum(seq_length)]
cat(inputs[,input],
    sep = "\n",
    file= paste0("../Data/GPS_input_",1,".fasta"))



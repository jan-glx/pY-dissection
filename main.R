# packages --------------------------------------------------------------------
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)

# read all grossman et al data from intact --------------------------------------------------------------------
rm(list = ls())
df <-
  read.table(
    '../Data/IM-22632.txt',sep = '\t',quote = "",stringsAsFactors = F,header =
      T,comment.char = ""
  )#,stringsAsFactors=F
#discard irrelevant columns
df %<>% select(
  Feature.s..interactor.A,
  Feature.s..interactor.B,
  Identification.method.participant.A,
  Identification.method.participant.B,
  Annotation.s..interactor.A,
  Annotation.s..interactor.B,
  Xref.s..interactor.B,
  Xref.s..interactor.A,
  Experimental.role.s..interactor.A,
  Experimental.role.s..interactor.B,
  ID.s..interactor.A = X.ID.s..interactor.A,
  ID.s..interactor.B,
  Host.organism.s.,
  Interaction.annotation.s.,
  Interaction.identifier.s.,
  Confidence.value.s.,
  Taxid.interactor = Taxid.interactor.B,
  Interaction.type.s.,
  Interaction.detection.method.s.
)

df %<>% gather(key,value, ends_with(".A"), ends_with(".B")) %>%
  extract(key,c("var","Interactor"),"(.*)\\.([AB])$") %>%
  spread(var,value)


df %<>% extract(Confidence.value.s.,c("confidence"),"^intact\\-miscore\\:(.*)$",convert =
                  T) %>%
  extract(
    Experimental.role.s..interactor,c("role"),'^psi\\-mi\\:"MI\\:(?:\\d+)"\\((prey|bait|neutral)(?: component)?\\)$',convert =
      T
  )  %>%
  separate(ID.s..interactor,c("Interactor_ID_db","Interactor_ID"),'\\:',convert =
             F)

# Map some ORFs that indentfied by intact IDs to their corresponding Uniprot IDs
to_map = df %>% filter(Interactor_ID_db == "intact") %>%
  select(Interactor_ID_db, Interactor_ID,Xref.s..interactor) %>%
  extract(
    Xref.s..interactor,c("From"),"refseq:([A-Z]+\\_\\d+(?:.\\d+)?)(?:[(|].*)?$",remove =
      T,convert = F
  ) %>%
  distinct()

url = 'http://www.uniprot.org/mapping/'
params =   c(
  'from', 'REFSEQ_NT_ID',
  'to', 'ACC',
  'format', 'tab',
  'query', paste(to_map$From, collapse = " ")
)
params = lapply(params, URLencode)
query = paste(params[c(T,F)], params[c(F,T)], sep = "=", collapse = "&")
mapping = read.table(paste0(url,'?',query), header = T, stringsAsFactors =
                       F)
mapping %<>% right_join(to_map)
mapping

df %<>% left_join(select(mapping,-From))
df %<>%
  mutate(Interactor_ID = ifelse(Interactor_ID_db == "intact",To, Interactor_ID)) %>%
  select(-Interactor_ID_db,-To)
# done mapping, save relevant data
write.table(
  df, "../Data/all_intact_data_from_grossmanEtAl.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)



# getting sequences for ORFs without Uniprotkb ID
url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta_cds_aa&id='
sequences = data.frame(seq = sapply((mapping %>% filter(is.na(
  To
)))$From,
function(From) {
  paste0(readLines(paste0(url,From))[-1], collapse = "")
}))
sequences$From = rownames(sequences)
rownames(sequences) <- NULL
# and do nothing with it

# Loading sequences and know modifications from uniprot for all ORFst --------------------------------------------------------------------
rm(list = ls())
df = read.table(
  "../Data/all_intact_data_from_grossmanEtAl.tsv", sep = "\t",header = T,stringsAsFactors = F
)

uniprot_ACCs = data.frame(ID = unique(toupper(df$Interactor_ID[!is.na(df$Interactor_ID)])),stringsAsFactors = F) %>%
  extract(ID,c("canonic","isoform"),"^([^\\-]*)(?:\\-(.*))?$", remove = F)
uniprot_ACCs_canonic = unique(uniprot_ACCs$canonic)

url = 'http://www.uniprot.org/uniprot/'
uniseq = NULL
for (uniprot_ACCs_canonic_i in split(uniprot_ACCs_canonic,floor(0:(length(
  uniprot_ACCs_canonic
) - 1) / 100))) {
  params =   c(
    'query', paste0('"accesion+',uniprot_ACCs_canonic_i,'"', collapse = "OR"),
    'columns', 'id,sequence,feature(MODIFIED RESIDUE),feature(ALTERNATIVE SEQUENCE),comment(ALTERNATIVE PRODUCTS)',
    'format', 'tab'
  )
  params = lapply(params, URLencode)
  query = paste(params[c(T,F)], params[c(F,T)], sep = "=", collapse = "&")
  uniseq = rbind(
    uniseq,read.table(
      paste0(url,'?',query), quote = "",comment.char = "", header = T, stringsAsFactors =
        F, sep = "\t"
    )
  )
  print(paste("fetched from uniprot:", nrow(uniseq)))
}

uniprot_data <-
  left_join(uniprot_ACCs, uniseq, by = c("canonic" = "Entry"))
dropped = (uniprot_data  %>% filter(stri_length(uniprot_data$Sequence) ==
                                      0))$ID
if (length(dropped) != 0) {
  warning(
    paste0(
      "dropped Entries with Interactor(s) ",
      paste0(dropped, collapse = ", "),
      " as the UniProt Isoform identifier could not be mapped to a sequence"
    )
  )
  uniprot_data %<>% filter(stri_length(uniprot_data$Sequence) != 0)
}

# done getting data fronm uniprot, save results
write.table(
  uniprot_data, "../Data/all_data_from_uniprot.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)

# parsing data from uniprot --------------------------------------------------------------------
rm(list = ls())
df = read.table(
  "../Data/all_intact_data_from_grossmanEtAl.tsv", sep = "\t",header = T,stringsAsFactors = F
)
uniprot_data = read.table(
  "../Data/all_data_from_uniprot.tsv", sep = "\t",header = T,colClasses = "character"
)

quotemeta <- function(string) {
  gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string)
}

parse_alternative_sequence <-
  function(alt_seq_field, alt_prod_field, isoform_id, canonical_sequence) {
    canonical_sequence <- strsplit(canonical_sequence,"")[[1]]
    pattern <-
      paste0(
        "Name=([^ ;]+)(?: \\{[^\\}]*\\})?;(?: Synonyms=[^;]+;)? IsoId=",quotemeta(isoform_id)
      )
    isoform = str_match(alt_prod_field,pattern)[2]
    if (is.na(isoform)) {
      warning(paste("isoform",isoform_id,"not found"))
      return(NA_character_)
    }
    
    pattern <-
      paste0(
        "VAR_SEQ (\\d+) (\\d+) (?:(Missing)|(\\w+) \\-\\> (\\w+)) \\(in [^\\)]*isoform ",
        quotemeta(isoform),"[^\\)]*\\)."
      )
    df = str_match_all(alt_seq_field,pattern)
    #cat(length(df))
    df = data.frame(df[[1]],stringsAsFactors = F)
    colnames(df) <-
      c("whole_match","start_base","end_base","Missing","From","To")
    if (nrow(df) == 0) {
      warning(paste("isoform",isoform,"(",isoform_id,")","not found"))
      return(NA_character_)
    }
    df %<>%  mutate_each(funs(as.numeric),start_base,end_base)
    
    df %<>% rowwise() %>% dplyr::mutate(mechismo_dif_string = ifelse(
      Missing == "Missing",
      paste0(
        canonical_sequence[start_base:end_base],
        as.character(start_base:end_base),
        "X",
        collapse = " "
      ),
      paste0(
        strsplit(From,"")[[1]],
        as.character(start_base:end_base),
        c(To,rep("X",end_base - start_base)),
        collapse = " "
      )
    ))
    paste(df$mechismo_dif_string, collapse = " ")
  }

uniprot_data =  uniprot_data %>% rowwise() %>%
  dplyr::mutate(mechismo_dif_string =
                  ifelse(
                    isoform == "",
                    "",
                    parse_alternative_sequence(
                      Alternative.sequence,Alternative.products..isoforms.,
                      ID,Sequence
                    )
                  )) %>% ungroup()

canonical2isoform <-
  function(canonical_sequence, ops) {
    sequence =  strsplit(canonical_sequence,"")[[1]]
    colnames(ops) <- c("full","from","position","del","to")
    ops %<>% data.frame(stringsAsFactors = F) %>% mutate(del = del == 'X',
                                                         position = as.numeric(position))
    sequence[ops$pos] = ops$to
    idx_map = cumsum(str_length(sequence))
    list(
      isoform_sequence = paste0(sequence,collapse = ""), canonical2isoform_idx =
        idx_map
    )
  }
canonical_sequence_list = uniprot_data$Sequence
op_list = str_match_all(uniprot_data$mechismo_dif_string, "([A-Z])(\\d+)(X(?![A-Z]))?([A-WYZ]?[A-Z]*)")
tmp = mapply(canonical2isoform , canonical_sequence_list, op_list)
uniprot_data$isoform_sequence = tmp[c(T,F)]
uniprot_data$canonical2isoform_idx = tmp[c(F,T)]

write.table(
  uniprot_data, "../Data/uniprotdata_4all_hits.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)

# extract(df$Interaction.annotation.s.,c(),
#   "[comment:\"\\\"Phosphorylation-dependent interaction. The interaction is only detected in the presence of the following kinases: (.+\\(.+\\),)*\\.\\\"
# \"|comment:\"\\\"The interaction failed to show if kinase-dead versions of ABL2 (P42684) or FYN (P06241) were used.\\\"\"|figure legend:Suppl. table S3, suppl. fig. S4|full coverage:Only protein-protein interactions|curation depth:imex curation)
#
#

#unique(unlist(lapply(df$Interaction.annotation.s.,function(x){strsplit(x,split="\\|")})))

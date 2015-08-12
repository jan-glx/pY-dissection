# packages --------------------------------------------------------------------
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)
library(data.table);options(datatable.showProgress=F)

# read all grossman et al data from intact --------------------------------------------------------------------
rm(list = ls())
dt <- fread('../Data/IM-22632.txt',select=c(
  "Interaction identifier(s)",
  "Interaction annotation(s)",
  "Confidence value(s)",
  "Interaction detection method(s)",
  "Feature(s) interactor A",
  "Feature(s) interactor B",
  "Annotation(s) interactor A" ,
  "Annotation(s) interactor B" ,
  "Xref(s) interactor A",
  "Xref(s) interactor B",
  "#ID(s) interactor A",
  "ID(s) interactor B", 
  "Experimental role(s) interactor A",
  "Experimental role(s) interactor B"
  ))[`Interaction detection method(s)`=='psi-mi:"MI:0397"(two hybrid array)'][, 
     `Interaction detection method(s)`:= NULL]#Since most interaction assays are orthogonal, each detecting its own 
# subset of true interactions, vali- dation of data sets means comparing validation rates rather than discarding 
# pairs that do not bind in the co-IP assay (Braun et al, 2009; Venkatesan et al, 2009). The validation rate of
# ~50% (Fig 3A) is similar for phospho-tyrosine-dependent and phospho-tyrosine- independent interactions. It is 
# lower than what was observed for more stable PPIs such as spliceosomal interactions (Hegele et al, 2012). 
# However, it compares well to validation rates reported for other representative sets of PPIs with this co-IP 
# assay (Braun et al, 2009; Weimann et al, 2013)

setnames(dt,"#ID(s) interactor A","ID(s) interactor A")
setnames(dt, make.names(colnames(dt)))

dt[,
   ':='(kinase_dep = str_match(Interaction.annotation.s.,
                paste0("Phosphorylation\\-dependent interaction\\. The interaction is only detected ",
                       "in the presence of the following kinases\\: ([^\\.]*)"))[,2],
        kinase_dep_binary=TRUE,
        confidence = as.numeric(str_match(Confidence.value.s.,"^intact\\-miscore\\:(.*)$")[,2]))
   ][,
     ':='(Interaction.annotation.s.= NULL,
     Confidence.value.s.=NULL)
     ]
dt[is.na(kinase_dep),
   ':='(kinase_dep="",
        kinase_dep_binary=FALSE)]
df = copy(dt)
df %<>% gather(key,value, ends_with(".A"), ends_with(".B")) 
df %<>% extract(key,c("var","Interactor"),"(.*)\\.([AB])$",remove = T)
df %<>% spread(var,value)
# df = copy(df)
# df[,
#    ':='(
#      role = str_match(
#        Experimental.role.s..interactor,
#        "^psi\\-mi\\:\"MI\\:(?:\\d+)\"\\((prey|bait|neutral)(?: component)?\\)$"
#      )[,2],
#      Experimental.role.s..interactor = NULL
#    )]
df %<>% extract(Experimental.role.s..interactor,c("role"),"^psi\\-mi\\:\"MI\\:(?:\\d+)\"\\((prey|bait|neutral)(?: component)?\\)$")

df %<>% separate(ID.s..interactor,c("Interactor_ID_db","Interactor_ID"),'\\:')


proteins = df[,
              .(ID_DB=Interactor_ID_db[1], ID=Interactor_ID[1], 
                refseq=str_match(Xref.s..interactor[1],"refseq:([A-Z]+\\_\\d+(?:.\\d+)?)(?:[(|].*)?$")[,2], 
                annotations__grossmann=Annotation.s..interactor[1]),
              ,by=.(Interactor_ID_db, Interactor_ID)]

df[,':='(Xref.s..interactor=NULL,
         Annotation.s..interactor=NULL,
         Interactor = NULL)]
backup.DT=copy(df)



df=copy(backup.DT)

df = dcast(df,Interaction.identifier.s. + kinase_dep + kinase_dep_binary + 
        confidence ~ role, 
      value.var = c("Interactor_ID_db", "Interactor_ID", "Feature.s..interactor"))

# done mapping, save relevant data
write.table(
  df, "../Data/interactions_grossmann.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)

# Map some ORFs that indentfied by intact IDs to their corresponding Uniprot IDs

url = 'http://www.uniprot.org/mapping/'
params =   c(
  'from', 'REFSEQ_NT_ID',
  'to', 'ACC',
  'format', 'tab',
                          'query', paste(proteins[Interactor_ID_db == "intact", refseq], collapse = " ")
)
params = lapply(params, URLencode)
query_url = paste0(url,'?',paste(params[c(T,F)], params[c(F,T)], sep = "=", collapse = "&"))
mapping = fread(query_url)
setkey(mapping,From)
mapping=unique(mapping)
setkey(proteins,refseq)
proteins[mapping,':='(ID_DB="uniprotkb",
                      ID=To)]
setkey(proteins,ID_DB,ID)

cat("removed the following entires from the analysis:\n")
show(proteins["intact"])

proteins=proteins["uniprotkb"][,ID_DB:=NULL]

proteins[, c('canonic','isoform') := data.table(str_match(proteins$ID, "^([^\\-]*)(?:\\-(.*))?$")[,2:3])
         ][, isoform=as.integer(isoform)]
setkey(proteins,"canonic")
uniprot_ACCs_canonic = unique(proteins)$canonic

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
uniseq <- as.data.table(uniseq)
# dropped = uniseq[Sequence=="",.(ID=Entry)]
# uniseq <- uniseq[Sequence!=""]
setkey(uniseq,"Entry")
setnames(uniseq,"Entry","canonic")

proteins=uniseq[proteins]
# 
# if (length(dropped) != 0) {
#   warning(
#     paste0(
#       "dropped Entries with Interactor(s) ",
#       paste0(dropped, collapse = ", "),
#       " as the UniProt Isoform identifier could not be mapped to a sequence"
#     )
#   )
# }

# done getting data fronm uniprot, save results
write.table(
  proteins, "../Data/all_data_from_uniprot.tsv", sep = "\t", row.names = F
)

# parsing data from uniprot --------------------------------------------------------------------
rm(list = ls())
df = read.table(
  "../Data/interactions_grossmann.tsv", sep = "\t",header = T,stringsAsFactors = F
)
uniprot_data = fread("../Data/all_data_from_uniprot.tsv",sep="\t")[,isoform:=as.integer(isoform)]


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

uniprot_data[is.na(isoform),
             mechismo_dif_string :=""]

uniprot_data[!is.na(isoform),
             mechismo_dif_string:= parse_alternative_sequence(Alternative.sequence,Alternative.products..isoforms.,
                                                               ID,Sequence),
             by = ID]


canonical2isoform <-
  function(canonical_sequence, ops) {
    sequence =  strsplit(canonical_sequence,"")[[1]]
    colnames(ops) <- c("full","from","position","del","to")
    ops %<>% data.frame(stringsAsFactors = F) %>% mutate(del = del == 'X',
                                                         position = as.numeric(position))
    sequence[ops$pos] = ops$to
    idx_map = cumsum(str_length(sequence))
    data.table(
      isoform_sequence = paste0(sequence,collapse = ""), canonical2isoform_idx =
        list(idx_map)
    )
  }
uniprot_data[, ops := str_match_all(uniprot_data$mechismo_dif_string, "([A-Z])(\\d+)(X(?![A-Z]))?([A-WYZ]?[A-Z]*)")]
uniprot_data[,
             c("isoform_sequence", "canonical2isoform_idx") := 
               do.call(rbind,Map(canonical2isoform , Sequence, ops))][, ops:=NULL]

canonical2isoform_idx = uniprot_data$canonical2isoform_idx
save(canonical2isoform_idx,file="../Data/canonical2isoform_idx.RData")
uniprot_data[,':='(canonical2isoform_idx=NULL,
                   Alternative.products..isoforms.=NULL,
                   Alternative.sequence=NULL)]


uniprot_data[,Modified.residue:=lapply(str_match_all(uniprot_data$Modified.residue,"MOD_RES (\\d+) \\d+ Phosphotyrosine"),
                             function(x){as.integer(x[,2])})
   ]
setnames(uniprot_data,"Modified.residue","Phosphotyrosins")
uniprot_data[,Phosphotyrosins:=sapply(Phosphotyrosins, function(x){paste0(x, collapse=",")})]

write.table(
  uniprot_data, "../Data/uniprotdata_4all_hits.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)

# clean interaction data ---------------------

uniprot_data = fread("../Data/uniprotdata_4all_hits.tsv",sep="\t")
setkey(uniprot_data,Interactor_ID_db,Interactor_ID)
df = fread("../Data/interactions_grossmann.tsv")
setkey(df,Interactor_ID_db_bait,Interactor_ID_bait)
df=df[uniprot_data,canonic_ID_bait:=canonic]
setkey(df,Interactor_ID_db_prey,Interactor_ID_prey)
df=df[uniprot_data,canonic_ID_prey:=canonic]
dropped=df[Interactor_ID_db_prey!="uniprotkb"|Interactor_ID_db_bait!="uniprotkb"]
write.table(
  dropped, file.path("..","Data",paste0("dropped_from_grossman.tsv")), sep = "\t", row.names = F, qmethod =
    "double"
)
df=df[Interactor_ID_db_prey=="uniprotkb"&Interactor_ID_db_bait=="uniprotkb"
      ][,
        ':='(Interactor_ID_db_prey=NULL,
             Interactor_ID_db_bait=NULL)]
setkey(df,canonic_ID_bait,canonic_ID_prey)

write.table(
  df, "../Data/valid_interactions_grossmann.tsv", sep = "\t", row.names = F, qmethod =
    "double"
)


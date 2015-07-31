library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)

df <- read.table('../Data/IM-22632.txt',sep='\t',quote="",stringsAsFactors=F,header=T,comment.char = "")#,stringsAsFactors=F
df %<>% select(Feature.s..interactor.A,
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
               Interaction.detection.method.s.)

df %<>% gather(key,value, ends_with(".A"), ends_with(".B")) %>%
  extract(key,c("var","Interactor"),"(.*)\\.([AB])$") %>%
  spread(var,value)


df %<>% mutate_each(funs(factor), Host.organism.s.:Identification.method.participant,
                    -Annotation.s..interactor,
                    -Interaction.annotation.s.,
                    -Confidence.value.s.)

df %<>% extract(Confidence.value.s.,c("confidence"),"^intact\\-miscore\\:(.*)$",convert=T) %>%
  extract(Experimental.role.s..interactor,c("role"),'^psi\\-mi\\:"MI\\:(?:\\d+)"\\((prey|bait|neutral)(?: component)?\\)$',convert=T)  %>%
  separate(ID.s..interactor,c("Interactor_ID_db","Interactor_ID"),'\\:',convert=F) 
  
to_map = df %>% filter(Interactor_ID_db == "intact") %>% 
  select(Interactor_ID_db, Interactor_ID,Xref.s..interactor) %>% 
  extract(Xref.s..interactor,c("From"),"refseq:([A-Z]+\\_\\d+(?:.\\d+)?)(?:[(|].*)?$",remove=T,convert=F) %>%
  distinct()

url = 'http://www.uniprot.org/mapping/'
params =   c(
  'from', 'REFSEQ_NT_ID',
  'to', 'ACC',
  'format', 'tab',
  'query', paste(to_map$From, collapse=" ")
)
params = lapply(params, URLencode)
query = paste(params[c(T,F)], params[c(F,T)], sep = "=", collapse = "&")
mapping = read.table(paste0(url,'?',query), header = T, stringsAsFactors=F)
mapping %<>% right_join(to_map)
mapping

url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta_cds_aa&id='
sequences = data.frame(seq=sapply((mapping %>% filter(is.na(To)))$From, 
                         function(From){paste0(readLines(paste0(url,From))[-1], collapse="")}))
sequences$From = rownames(sequences)
rownames(sequences) <- NULL













mapping


df %<>% left_join(select(mapping, -From))

df %<>% 
  mutate(Interactor_ID = ifelse(Interactor_ID_db == "intact",To, Interactor_ID)) %>%
  select(-Interactor_ID_db,-To)
summary(df)
unique(df$Interactor_ID)
# extract(df$Interaction.annotation.s.,c(),
#   "[comment:\"\\\"Phosphorylation-dependent interaction. The interaction is only detected in the presence of the following kinases: (.+\\(.+\\),)*\\.\\\"
# \"|comment:\"\\\"The interaction failed to show if kinase-dead versions of ABL2 (P42684) or FYN (P06241) were used.\\\"\"|figure legend:Suppl. table S3, suppl. fig. S4|full coverage:Only protein-protein interactions|curation depth:imex curation)
# 
# 

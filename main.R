if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr,tidyr,dplyr,ggplot2)

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
               Interaction.detection.method.s.) #%>%
#   transform(Xref.s..interactor.A = as.character(Xref.s..interactor.A),
#             Xref.s..interactor.B = as.character(Xref.s..interactor.B),
#             Confidence.value.s. = as.character(Confidence.value.s.),
#             Interaction.annotation.s. = as.character(Interaction.annotation.s.),
#             Annotation.s..interactor.A = as.character(Annotation.s..interactor.A),
#             Annotation.s..interactor.B = as.character(Annotation.s..interactor.B),
#             Interaction.identifier.s. = as.character(Interaction.identifier.s.))
summary(df)

df %<>% gather(key,value, ends_with(".A"), ends_with(".B")) %>%
  extract(key,c("var","Interactor"),"(.*)\\.([AB])$") %>%
  spread(var,value)

summary(df)

df %<>% mutate_each(funs(factor), Host.organism.s.:Identification.method.participant,
                    -Annotation.s..interactor,
                    -Interaction.annotation.s.,
                    -Confidence.value.s.)


summary(df)



df %<>% extract(Confidence.value.s.,c("confidence"),"^intact-miscore\\:(.*)$",convert=T) 

summary(df)

library(data.table)
options(datatable.prettyprint.char = 25L)
library(binom)



source("misc.R")
grossmann  = fread('../Data/IM-22632.txt',
                   select=c(
                     "Interaction identifier(s)",
                     "Interaction annotation(s)",
                     "Interaction detection method(s)",
                     "#ID(s) interactor A",
                     "ID(s) interactor B", 
                     "Experimental role(s) interactor A",
                     "Experimental role(s) interactor B",
                     "Creation date"
                     ), quote = ""
                   )[`Interaction detection method(s)`=='psi-mi:"MI:0397"(two hybrid array)'
                     ][,`Interaction detection method(s)`:= NULL
                       ]
setnames(grossmann,c("#ID(s) interactor A","ID(s) interactor B"),c("id.A","id.B"))
 
# switch interactor A  nd B such that A is always prey:      
grossmann[`Experimental role(s) interactor A`=="psi-mi:\"MI:0496\"(bait)",
          ':='(id.A=id.B,
               id.B=id.A
          )
          ][]

# publication day of grossmann data
day_zero = min(grossmann[, `Creation date`])

grossmann = grossmann[str_detect(`Interaction annotation(s)`,"Phosphorylation-dependent interaction"),
                       .(id.A,id.B)] # only keep pY dep PPIs and only the IDs

grossmann %<>% extract(id.A,"id.A","(intact.*|[^-]*)(?:-[^-]*)*") %>% extract(id.B,"id.B","(intact.*|[^-]*)(?:-[^-]*)*") %>% unique()

preys = unique(grossmann[,id.A])
baits = unique(grossmann[,id.B])
preysAndBaits = unique(c(preys,baits))


intact = 
  fread("../Data/intact/intact.txt", 
        select = c(
          "Interaction identifier(s)",
          "Interaction annotation(s)",
          "Interaction parameter(s)",
          "Interaction detection method(s)",
          "Publication Identifier(s)",
          "#ID(s) interactor A",
          "Creation date",
          "Negative",
          "ID(s) interactor B"
          ), verbose = T, quote = ""
        )[`Creation date`<day_zero & Negative=="false"]
setnames(intact,c("#ID(s) interactor A","ID(s) interactor B"),c("id.A","id.B"))
intact=intact[id.A %in% preysAndBaits][id.B %in% preysAndBaits]
intact=rbind(intact,copy(intact)[,':='(id.A=id.B,id.B=id.A)])[id.A %in% preys][id.B %in% baits]
setkey(intact, id.A,id.B)
intact = unique(intact)

both_hit = intact[grossmann, nomatch=0]

cat("Of the ", nrow(grossmann), " pY dependent PPIs found by Grossmann et al. ", 
    nrow(both_hit), " (",nrow(both_hit)/nrow(grossmann)*100,"%) are found elsewhere as PPI on intact.",
    "\nWhile of the ",length(preys)*length(baits), " PPP that Grossmann et al. could have found as pY PPI ", 
    nrow(intact), " (",nrow(intact)/(length(preys)*length(baits))*100,"%)are found elsewhere as PPI on intact.",
    sep="")
#Of the 294 pY dependent PPIs found by Grossmann et al. 18 (6.122449%) are found elsewhere as PPI on intact.
#While of the 11033 PPP that Grossmann et al. could have found as pY PPI 119 (1.078582%)are found elsewhere as PPI on intact.


preys_in_intact = intact[, .N, keyby=id.A][preys][is.na(N),N:=0][]
tmp <- binom.bayes(mean(preys_in_intact[,N]),length(preys))
preys_in_intact[,p:=binom.bayes(N,length(preys), prior.shape1=tmp$shape1, prior.shape2=tmp$shape2)$mean]

baits_in_intact = intact[, .N, keyby=id.B][baits][is.na(N),N:=0][]
tmp <- binom.bayes(mean(baits_in_intact[,N]),length(baits))
baits_in_intact[,p:=binom.bayes(N,length(baits), prior.shape1=tmp$shape1, prior.shape2=tmp$shape2)$mean]

mean(CJ(preys_in_intact[,p], baits_in_intact[,p])[,V1*V2])

baits_in_intact[preys_in_intact[grossmann,on="id.A",nomatch=0], on="id.B",nomatch=0][,.(E=sum(p*i.p),P=mean(p*i.p))]



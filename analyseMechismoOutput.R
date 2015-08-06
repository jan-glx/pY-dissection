printSetDiffSizes <- function(A,B) {
  nA = deparse(substitute(A))
  nB = deparse(substitute(B))
  uA = unique(A)
  uB = unique(B)
  cat(paste0("In ",nA,": ",length(uA),"\n"))
  cat(paste0("In ",nA, " but not in ", nB ,": ",length(setdiff(uA,uB)),"\n"))
  cat(paste0("In ",nB,": ",length(uB),"\n"))
  cat(paste0("In ",nB, " but not in ", nA ,": ",length(setdiff(uB,uA)),"\n"))
  cat(paste0("In ",nA, " and in ", nB ,": ",length(intersect(uB,uA)),"\n"))
}


analyse_mechismo_output <-
  function (df, uniprot_data, folder, filename) {
    sink(file.path(folder,paste0(filename,".log")))
    
    mechismo_out = read.table(
      file.path(folder,filename), sep = "\t",header = T,stringsAsFactors = F
    )
    mechismo_out %<>% select(primary_id_a1, primary_id_b1,pos_a1:mismatch,conf:ie)
    mechismo_out %<>% distinct(primary_id_a1,primary_id_b1)
    mechismo_out %<>% filter(primary_id_b1 != "")
    
    to_drop = (mechismo_out  %>% filter(as.logical(mismatch)))$primary_id_a1 %>% unique()
    dropped = intersect(uniprot_data$canonic,to_drop)
    
    if (length(dropped) != 0) {
      cat("Waring: dropped Entries with Interactor(s):\n ", paste0(dropped, collapse = ", "),
          "\nbecause mechismo indicates that it used a different sequence as basis. (older uniprot sequence version)\n\n"
      )
      uniprot_data %<>% filter(not(canonic %in% dropped))
    }
    mechismo_out  %<>% filter(not(as.logical(mismatch))) %>% select(-mismatch)
    
    
    to_drop = (
      mechismo_out  %>%
        filter(mut_a1 != "Yp") %>%
        filter(primary_id_b1 %in% uniprot_data$canonic)
    )$primary_id_a1 %>% unique()
    dropped = intersect(uniprot_data$canonic,to_drop)
    
    if (length(dropped) != 0) {
      cat("Warning: WOuld like to have dropped Entries with Interactor(s)\n", paste0(dropped, collapse = ", "),
          "\nbecause mechismo indicatet that isoforms differ in",
          "their interaction partners (subset of the subset of proteins that showed up in the study).\n\n"
      )
      #uniprot_data %<>% filter(not(canonic %in% dropped))
    }
    mechismo_out  %<>% filter(mut_a1 == "Yp") %>% select(-res_a1,-mut_a1)
    
    
    write.table(
      mechismo_out, file.path(folder,paste0(filename,"_clean.tsv")),
      sep = "\t", row.names = F, qmethod =
        "double"
    )
    
    
    
    intersect(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1)
    
    (printSetDiffSizes(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_a1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_b1))
    (printSetDiffSizes(
      uniprot_data$canonic,c(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1)
    ))
    
    # sdf ---------------------------------
    
    
    a = df %>% filter(Interactor_ID %in% uniprot_data$ID)
    a %<>% distinct(Interaction.identifier.s.,Interactor) %>% select(-Interactor)
    a %<>% gather(key,value,Annotation.s..interactor:Xref.s..interactor,-role)
    a %<>% unite(key2,key,role)
    a %<>% spread(key2,value)
    a %<>% filter(not(or(
      is.na(Interactor_ID_bait),is.na(Interactor_ID_prey)
    )))
    a %<>% unite(bp_pair_id,Interactor_ID_bait, Interactor_ID_prey, remove =
                   F)
    
    (printSetDiffSizes(a$Interactor_ID_bait,a$Interactor_ID_prey))
    (printSetDiffSizes(a$Interactor_ID_bait,mechismo_out$primary_id_a1))
    cat(printSetDiffSizes(a$Interactor_ID_bait,mechismo_out$primary_id_b1))
    (printSetDiffSizes(a$Interactor_ID_prey,mechismo_out$primary_id_a1))
    (printSetDiffSizes(a$Interactor_ID_prey,mechismo_out$primary_id_b1))
    
    
    #   very_good= mechismo_out %>%
    #     filter(and(mechismo_out$primary_id_b1 %in% a$Interactor_ID_bait,
    #                mechismo_out$primary_id_a1 %in% a$Interactor_ID_prey))
    #
    #   good= mechismo_out %>%
    #     filter(or(and(mechismo_out$primary_id_b1 %in% a$Interactor_ID_bait,
    #                mechismo_out$primary_id_a1 %in% a$Interactor_ID_prey),
    #               and(mechismo_out$primary_id_b1 %in% a$Interactor_ID_prey,
    #                   mechismo_out$primary_id_a1 %in% a$Interactor_ID_bait)))
    mechismo_out %<>% unite(ab_pair_id,primary_id_a1, primary_id_b1, remove =
                              F) %>%
      unite(ba_pair_id,primary_id_b1, primary_id_a1, remove = F) %>%
      gather(pair_type,pair_id,ab_pair_id,ba_pair_id)# %>%
      #filter(not(and(pair_type == "ba_pair_id", primary_id_a1 == primary_id_b1))) # no reason to filter
    
    joined = inner_join(a,mechismo_out,by = c("bp_pair_id" = "pair_id"))
    write.table(
      joined, file.path(folder,paste0(filename,"_joined.tsv")), sep = "\t", row.names = F, qmethod =
        "double"
    )
    
    
    joined_sum = joined %>% group_by(Interaction.identifier.s.) %>%
      dplyr::summarize(
        ie_sum = sum(ie),
        ie_max = max(ie),
        sites = paste0(user.input,": ",ie,"(",conf,")", collapse =
                         "; "),
        n = n()
      ) %>%
      left_join(a) %>%
      arrange(ie_sum) %>%
      mutate(n_kinases=str_count(kinase.dep,"\\)"),
             kinase.dep_binary=kinase.dep!="") 
    
    write.table(
      joined_sum, file.path(folder,paste0(filename,"_joined_summary.tsv")), sep = "\t", row.names = F, qmethod =
        "double"
    )
    ggplot(joined_sum ,aes(x=n_kinases,y=ie_sum))+geom_point(alpha=0.3,size=5)
    ggsave(file.path(folder,paste0(filename,"_ie_vs_kin.dep.png")))
    print(with(joined_sum, wilcox.test(ie_sum~kinase.dep_binary)))
    sink()
  }


# packages --------------------------------------------------------------------
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

df = read.table(
  "../Data/all_intact_data_from_grossmanEtAl.tsv", sep = "\t",header = T,stringsAsFactors = F
)
uniprot_data = read.table(
  "../Data/uniprotdata_4all_hits.tsv", sep = "\t",header = T,colClasses = "character"
)
load(file = "../Data/canonical2isoform_idx.RData")

"../Data/Jq_IWhrTTn.site_table.tsv"
folder=file.path("..","Data","Out")
ffiles <- list.files( path=folder, pattern = "\\.site\\_table\\.tsv$")
for( file in ffiles){
  analyse_mechismo_output(df, uniprot_data, folder, file)
}
cat("Done Analysing output of ",ffiles, "\n")

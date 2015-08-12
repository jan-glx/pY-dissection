
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
    #  filename="v7BAv7h5Pt.site_table.tsv"
    
    mechismo_out = fread(file.path(folder,filename),
                         select=c("primary_id_a1", "primary_id_b1", "pos_a1", "mut_a1", 
                                  "user input", "mismatch", "conf", "ie","name_b1","idcode","iupred"))
    mechismo_out = mechismo_out[(name_b1=="[PROT]" ) | (primary_id_b1 %in% unique(uniprot_data$canonic))]
    setkey(mechismo_out,primary_id_a1,primary_id_b1,pos_a1,name_b1)
    mechismo_out = mechismo_out[mechismo_out[,.(pdb_ids=paste(idcode, collapse=", "),
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
    
    
    to_drop = mechismo_out[mut_a1 != "Yp" & primary_id_b1 %in% uniprot_data$canonic,unique(primary_id_a1)]
    
    if (length(to_drop) != 0) {
      cat("Warning: Would like to have dropped Entries with Interactor(s)\n", paste0(dropped, collapse = ", "),
          "\nbecause mechismo indicatet that isoforms differ in",
          "their interaction partners (subset of the subset of proteins that showed up in the study).\n\n"
      )
    }
    mechismo_out = mechismo_out[mut_a1 == "Yp"][,mut_a1:=NULL]
    
    
    write.table(
      mechismo_out, file.path(folder,paste0(filename,"_clean.tsv")),
      sep = "\t", row.names = F, qmethod =
        "double"
    )
    
    (printSetDiffSizes(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_a1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_b1))
    (printSetDiffSizes(
      uniprot_data$canonic,c(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1)
    ))
    
    # sdf ---------------------------------
    setkey(mechismo_out,primary_id_a1,primary_id_b1)
    mechismo_interactions=mechismo_out[unique(unique(mechismo_out)[primary_id_b1!="",
                                                          .(primary_id_a1=rbind(primary_id_a1,primary_id_b1),
                                                            primary_id_b1=rbind(primary_id_b1,primary_id_a1))
                                                          ]),
                                       nomatch=0
                                       ]
    
    mechismo_interactions=mechismo_interactions[,
                                       .(ie_sum = sum(ie),
                                         ie_max = max(ie),
                                         sites = paste0(`user input`,": ",ie,"(",conf,")", collapse =
                                                          "; "),
                                         n = .N),
                                       by=.(primary_id_a1,primary_id_b1)]
    df[,n_kinases:=str_count(kinase_dep,"\\)")]
    hits=df[mechismo_interactions,nomatch=0]
    
    
    write.table(
      hits, file.path(folder,paste0(filename,"_joined_summary.tsv")), sep = "\t", row.names = F, qmethod =
        "double"
    )
    ggplot(hits ,aes(x=kinase_dep_binary,y=ie_sum))+geom_point(alpha=0.3,size=5)
    
    ggsave(file.path(folder,paste0(filename,"_ie_vs_kin.dep.png")))
    print(with(hits, wilcox.test(ie_sum~kinase_dep_binary)))
    sink()
  }


# packages --------------------------------------------------------------------
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

uniprot_data = fread("../Data/uniprotdata_4all_hits.tsv",sep="\t")
df = fread("../Data/valid_interactions_grossmann.tsv")
setkey(df,canonic_ID_bait,canonic_ID_prey)

folder=file.path("..","Data","Out")
ffiles <- list.files( path=folder, pattern = "\\.site\\_table\\.tsv$")
for( file in ffiles){
  analyse_mechismo_output(df, uniprot_data, folder, file)
}
cat("Done Analysing output of ",ffiles, "\n")

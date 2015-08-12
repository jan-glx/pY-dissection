rm(list=ls())
analyse_mechismo_output <-
  function (df, uniprot_data, folder, filename) {
    sink(file.path(folder,paste0(filename,".log")))
    #  filename="enFcexLbnj.site_table.tsv"
    
    mechismo_out = fread(file.path(folder,filename),
                         select=c("primary_id_a1", "primary_id_b1", "pos_a1", "mut_a1", 
                                  "user input", "mismatch", "conf", "ie","name_b1","idcode","iupred"))[primary_id_a1!=""]
    
    mechismo_out = mechismo_out[name_b1=="" | (name_b1=="[PROT]" ) | (primary_id_b1 %in% unique(uniprot_data$canonic))]
    setkey(mechismo_out,primary_id_a1,primary_id_b1,pos_a1,name_b1)
    mechismo_out = mechismo_out[mechismo_out[,.(pdb_ids=list(unique(idcode)),
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
      cat("Warning: Would like to have dropped Entries with Interactor(s)\n", paste0(unique(dropped), collapse = ", "),
          "\nbecause mechismo indicatet that isoforms differ in",
          "their interaction partners (subset of the subset of proteins that showed up in the study).\n\n"
      )
    }
    mechismo_out = mechismo_out[mut_a1 == "Yp"][,mut_a1:=NULL]
    
    proteins = mechismo_out[,
                            .(iupred = sum(mean_iupred)),
                            by=primary_id_a1]
    proteins[,.(.N,IDs=list(primary_id_a1)),keyby=iupred]
    has_unstructured_Y = proteins[iupred>0,primary_id_a1]
    n_int_with_unstructured_Y= df[canonic_ID_bait %in% has_unstructured_Y | 
                                    canonic_ID_prey %in% has_unstructured_Y,.N]
    cat("Of the ",nrow(df)," intractions ", n_int_with_unstructured_Y, 
        "contain tyrosine(s) in a region predicted as unstructured and ",
        nrow(df)-n_int_with_unstructured_Y," do not.\n")
    

    
    (printSetDiffSizes(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_a1))
    (printSetDiffSizes(uniprot_data$canonic,mechismo_out$primary_id_b1))
    (printSetDiffSizes(
      uniprot_data$canonic,c(mechismo_out$primary_id_a1,mechismo_out$primary_id_b1)
    ))
    
    # sdf ---------------------------------
    setkey(mechismo_out,primary_id_a1,primary_id_b1)
    helper_idx=unique(mechismo_out)[primary_id_b1!="",
                                               .(primary_id_a1, primary_id_b1)
                                              ][, .(primary_id_a1=rbind(primary_id_a1,primary_id_a1),
                                                    primary_id_b1=rbind(primary_id_b1,primary_id_b1),
                                                    new_primary_id_a1=rbind(primary_id_a1,primary_id_b1),
                                                    new_primary_id_b1=rbind(primary_id_b1,primary_id_a1))
                                                ]
    setkey(helper_idx)
    helper_idx=unique(helper_idx)
    setkey(helper_idx,primary_id_a1,primary_id_b1)
    mechismo_interactions=mechismo_out[helper_idx,nomatch=0][,':='(primary_id_a1=new_primary_id_a1,
                                                         primary_id_b1=new_primary_id_b1,
                                                         new_primary_id_a1=NULL,
                                                         new_primary_id_b1=NULL
                                                         )]
    setkey(mechismo_interactions,primary_id_a1,primary_id_b1)

    
    mechismo_interactions=mechismo_interactions[,
                                       .(ie_sum = sum(ie),
                                         ie_max = max(ie),
                                         ie_mean = mean(ie),
                                         ie_median = median(ie),
                                         sites = paste0(`user input`,": ",ie,"(",conf,")", collapse =
                                                          "; "),
                                         pdb_ids = list(unique(do.call(c,pdb_ids))),
                                         n = .N),
                                       by=.(primary_id_a1,primary_id_b1)]
    
    hits=df[mechismo_interactions,nomatch=0]
    
    
    ggplot(hits ,aes(x=kinase_dep_binary,y=ie_sum))+geom_point(alpha=0.3,size=5)
    ggsave(file.path(folder,paste0(filename,"_ie_vs_kin.dep.png")))
    print(with(hits, wilcox.test(ie_sum~kinase_dep_binary)))
    
    
    hits[,pdb_ids:=sapply(pdb_ids,function(x){paste0(x, collapse=", ")})]
    
    setkey(hits,Interactor_ID_bait)
    uni_colnames = colnames(uniprot_data)
    joined=setkey(hits[setnames(uniprot_data,paste0(uni_colnames,"_bait")),
                       nomatch=0
                       ],
                  Interactor_ID_prey)[setnames(uniprot_data,paste0(uni_colnames,"_prey")),
                                      nomatch=0
                                      ]
    write.table(
      joined, file.path(folder,paste0(filename,"_joined_joined.tsv")), sep = "\t", row.names = F, qmethod =
        "double"
    )
    
    write.table(
      hits, file.path(folder,paste0(filename,"_joined_summary.tsv")), sep = "\t", row.names = F, qmethod =
        "double"
    )
    mechismo_out[,pdb_ids:=sapply(pdb_ids,function(x){paste0(x, collapse=", ")})]
    write.table(
      mechismo_out, file.path(folder,paste0(filename,"_clean.tsv")),
      sep = "\t", row.names = F, qmethod =
        "double"
    )
    sink()
  }


# packages --------------------------------------------------------------------
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

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

uniprot_data = fread("../Data/uniprotdata_4all_hits.tsv",sep="\t")
setkey(uniprot_data,ID)
#uniprot_data[,Phosphotyrosins:=lapply(strsplit(Phosphotyrosins,","),as.integer)]
df = fread("../Data/valid_interactions_grossmann.tsv")
df[,n_kinases:=str_count(kinase_dep,"\\)")]
setkey(df,canonic_ID_bait,canonic_ID_prey)

folder=file.path("..","Data","Out")
ffiles <- list.files( path=folder, pattern = "\\.site\\_table\\.tsv$")
for( file in ffiles){
  analyse_mechismo_output(df, uniprot_data, folder, file)
}
cat("Done Analysing output of ",ffiles, "\n")

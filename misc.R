library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

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

library(jsonlite)
fwriteplus <- function(dt, file) {
  coltypes=sapply(dt, class)
  list_cols=names(coltypes)[coltypes=="list"]
  ii=list_cols[1]
  setnames(dt,list_cols,paste0(".",list_cols))
  invisible(lapply(list_cols,function(ii)dt[,(ii):=sapply(dt[[paste0(".",ii)]],serializeJSON)]))
  write.table(dt[,names(coltypes),with=FALSE],file,sep="\t",quote=F,row.names = F)
  invisible(lapply(list_cols,function(ii)dt[,(ii):=NULL]))
  setnames(dt,paste0(".",list_cols),list_cols)
}

freadplus <- function(...){
  dt <- fread(...)
  coltypes=sapply(dt, class)
  char_cols=names(coltypes)[coltypes=="character"]
  maybeJSON_cols=char_cols[sapply(dt[1,char_cols,with=FALSE], function(x)validate(x))]
  if(length(maybeJSON_cols)!=0){
    JSON_cols=maybeJSON_cols[sapply(dt[,maybeJSON_cols,with=FALSE], 
                               function(x)validate(paste0("[",paste0(x,collapse=","),"]")))]
    invisible(lapply(JSON_cols,function(ii)dt[,(ii):=lapply(dt[[ii]],unserializeJSON)][]))
  }
  dt
}
  
fwritelist <- function(dt, file) {
  coltypes=sapply(dt, class)
  list_cols=names(coltypes)[coltypes=="list"]
  ii=list_cols[1]
  if(length(list_cols)>0) { # rsucks
    setnames(dt,list_cols,paste0(".",list_cols))
    invisible(lapply(list_cols,function(ii)dt[,(ii):=sapply(dt[[paste0(".",ii)]],paste0,collapse=";")])) 
  }
  write.table(dt[,names(coltypes),with=FALSE],file,sep="\t",quote=F,row.names = F)
  if(length(list_cols)>0) {
    invisible(lapply(list_cols,function(ii)dt[,(ii):=NULL]))
    setnames(dt,paste0(".",list_cols),list_cols)
  }
}

freadlist <- function(...,sep2=";"){
  dt <- fread(...)
  coltypes=sapply(dt, class)
  char_cols=names(coltypes)[coltypes=="character"]
  maybelist_cols=char_cols[sapply(dt[1,char_cols,with=FALSE], function(x)substr(x,1,1)!='"')]
  if(length(maybelist_cols)!=0){
    list_cols=maybelist_cols[sapply(dt[,maybelist_cols,with=FALSE], 
                                    function(x)any(str_detect(x,sep2)))]
    invisible(lapply(list_cols,function(ii)dt[,(ii):=strsplit(dt[[ii]],split=sep2)][]))
  }
  dt
}
  
  
  
  
AAs=c("_", "A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
      "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"
)


#sort(unique(sites[,unlist(strsplit(neighbourhood,""))]))
blosum62=structure(c(4L, -1L, -2L, -2L, 0L, -1L, -1L, 0L, -2L, -1L, -1L, 
                     -1L, -1L, -2L, -1L, 1L, 0L, -3L, -2L, 0L, -2L, -1L, 0L, -4L, 
                     -1L, 5L, 0L, -2L, -3L, 1L, 0L, -2L, 0L, -3L, -2L, 2L, -1L, -3L, 
                     -2L, -1L, -1L, -3L, -2L, -3L, -1L, 0L, -1L, -4L, -2L, 0L, 6L, 
                     1L, -3L, 0L, 0L, 0L, 1L, -3L, -3L, 0L, -2L, -3L, -2L, 1L, 0L, 
                     -4L, -2L, -3L, 3L, 0L, -1L, -4L, -2L, -2L, 1L, 6L, -3L, 0L, 2L, 
                     -1L, -1L, -3L, -4L, -1L, -3L, -3L, -1L, 0L, -1L, -4L, -3L, -3L, 
                     4L, 1L, -1L, -4L, 0L, -3L, -3L, -3L, 9L, -3L, -4L, -3L, -3L, 
                     -1L, -1L, -3L, -1L, -2L, -3L, -1L, -1L, -2L, -2L, -1L, -3L, -3L, 
                     -2L, -4L, -1L, 1L, 0L, 0L, -3L, 5L, 2L, -2L, 0L, -3L, -2L, 1L, 
                     0L, -3L, -1L, 0L, -1L, -2L, -1L, -2L, 0L, 3L, -1L, -4L, -1L, 
                     0L, 0L, 2L, -4L, 2L, 5L, -2L, 0L, -3L, -3L, 1L, -2L, -3L, -1L, 
                     0L, -1L, -3L, -2L, -2L, 1L, 4L, -1L, -4L, 0L, -2L, 0L, -1L, -3L, 
                     -2L, -2L, 6L, -2L, -4L, -4L, -2L, -3L, -3L, -2L, 0L, -2L, -2L, 
                     -3L, -3L, -1L, -2L, -1L, -4L, -2L, 0L, 1L, -1L, -3L, 0L, 0L, 
                     -2L, 8L, -3L, -3L, -1L, -2L, -1L, -2L, -1L, -2L, -2L, 2L, -3L, 
                     0L, 0L, -1L, -4L, -1L, -3L, -3L, -3L, -1L, -3L, -3L, -4L, -3L, 
                     4L, 2L, -3L, 1L, 0L, -3L, -2L, -1L, -3L, -1L, 3L, -3L, -3L, -1L, 
                     -4L, -1L, -2L, -3L, -4L, -1L, -2L, -3L, -4L, -3L, 2L, 4L, -2L, 
                     2L, 0L, -3L, -2L, -1L, -2L, -1L, 1L, -4L, -3L, -1L, -4L, -1L, 
                     2L, 0L, -1L, -3L, 1L, 1L, -2L, -1L, -3L, -2L, 5L, -1L, -3L, -1L, 
                     0L, -1L, -3L, -2L, -2L, 0L, 1L, -1L, -4L, -1L, -1L, -2L, -3L, 
                     -1L, 0L, -2L, -3L, -2L, 1L, 2L, -1L, 5L, 0L, -2L, -1L, -1L, -1L, 
                     -1L, 1L, -3L, -1L, -1L, -4L, -2L, -3L, -3L, -3L, -2L, -3L, -3L, 
                     -3L, -1L, 0L, 0L, -3L, 0L, 6L, -4L, -2L, -2L, 1L, 3L, -1L, -3L, 
                     -3L, -1L, -4L, -1L, -2L, -2L, -1L, -3L, -1L, -1L, -2L, -2L, -3L, 
                     -3L, -1L, -2L, -4L, 7L, -1L, -1L, -4L, -3L, -2L, -2L, -1L, -2L, 
                     -4L, 1L, -1L, 1L, 0L, -1L, 0L, 0L, 0L, -1L, -2L, -2L, 0L, -1L, 
                     -2L, -1L, 4L, 1L, -3L, -2L, -2L, 0L, 0L, 0L, -4L, 0L, -1L, 0L, 
                     -1L, -1L, -1L, -1L, -2L, -2L, -1L, -1L, -1L, -1L, -2L, -1L, 1L, 
                     5L, -2L, -2L, 0L, -1L, -1L, 0L, -4L, -3L, -3L, -4L, -4L, -2L, 
                     -2L, -3L, -2L, -2L, -3L, -2L, -3L, -1L, 1L, -4L, -3L, -2L, 11L, 
                     2L, -3L, -4L, -3L, -2L, -4L, -2L, -2L, -2L, -3L, -2L, -1L, -2L, 
                     -3L, 2L, -1L, -1L, -2L, -1L, 3L, -3L, -2L, -2L, 2L, 7L, -1L, 
                     -3L, -2L, -1L, -4L, 0L, -3L, -3L, -3L, -1L, -2L, -2L, -3L, -3L, 
                     3L, 1L, -2L, 1L, -1L, -2L, -2L, 0L, -3L, -1L, 4L, -3L, -2L, -1L, 
                     -4L, -2L, -1L, 3L, 4L, -3L, 0L, 1L, -1L, 0L, -3L, -4L, 0L, -3L, 
                     -3L, -2L, 0L, -1L, -4L, -3L, -3L, 4L, 1L, -1L, -4L, -1L, 0L, 
                     0L, 1L, -3L, 3L, 4L, -2L, 0L, -3L, -3L, 1L, -1L, -3L, -1L, 0L, 
                     -1L, -3L, -2L, -2L, 1L, 4L, -1L, -4L, 0L, -1L, -1L, -1L, -2L, 
                     -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -2L, 0L, 0L, -2L, 
                     -1L, -1L, -1L, -1L, -1L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, 
                     -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, -4L, 
                     -4L, -4L, -4L, 1L), .Dim = c(24L, 24L), .Dimnames = list(c("A", 
                                                                                "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", 
                                                                                "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "GAP"), c("A", "R", 
                                                                                                                                       "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", 
                                                                                                                                       "S", "T", "W", "Y", "V", "B", "Z", "X", "GAP")))
a=prcomp(blosum62[1:20,1:20], scale. = TRUE)$rotation
sapply(2:20,function(x)nrow(unique.matrix(a[,1:x]>0)))
sapply(2:20,function(x)nrow(unique.matrix((a[,1:x]>-0.1)+(a[,1:x]>0.1))))
apply(a,2,function(x)rownames(a)[order(x)])

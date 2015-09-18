library(data.table)
library(stringr)
folder=file.path("..","Data")
filename <- "Out/v7BAv7h5Pt.site_table.tsv"

mechismo_out = fread(file.path(folder,filename))
setkey(mechismo_out,primary_id_a1)

grossmann_dt = fread(file.path(folder, "all_intact_data_from_grossmanEtAl.tsv"))[
  ,':='(Xrefs=str_split(Xref.s..interactor,"\\|"),
        Xref.s..interactor=NULL,
        ID_canonic=str_replace(Interactor_ID,"\\-.*",""))]

Xrefs_dt=data.table(str_match(unique(unlist(grossmann_dt$Xrefs)),
                      "([^\\:]+)\\:([^\\()]+)(?:\\(([^\\)]*)\\))?"))[,ID:=.I]

setnames(Xrefs_dt,c("complete","database","identifier","description","ID"))
setkey(Xrefs_dt,complete)

grossmann_dt[,Xrefs:=lapply(Xrefs,function(Xrefsl){Xrefs_dt[Xrefsl,ID]})]
baits=unique(grossmann_dt[role=="bait",ID_canonic])
preys=unique(grossmann_dt[role=="prey",ID_canonic])



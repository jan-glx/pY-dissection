require(RMySQL)
library(data.table)

dt = fread("../Data/all_intact_data_from_grossmanEtAl.tsv")

dt = dt[, 
        .(ID_prey = Interactor_ID[role == "prey"],
          ID_bait = Interactor_ID[role == "bait"],
          kinase.dep_binary = kinase.dep == "",
          n_dep_kinases = str_count(kinase.dep,"\\)")
        ),
        by = .(Interaction.identifier.s.,kinase.dep,confidence)]


# mechismo_base = fread("../Data/getknownpairs.sql.tsv") # alternative to querying the database
fistdb=dbConnect(RMySQL::MySQL(),"fistdb",user="anonymous")

pair_list=paste0("('",dbEscapeStrings(fistdb,dt$ID_prey),
              "', '",dbEscapeStrings(fistdb,dt$ID_bait),"')",collapse=", ")

part1="
SELECT sa1.id idA,
  sa1.primary_id  acA,
  sb1.id          idB,
  sb1.primary_id  acB,
  GROUP_CONCAT(DISTINCT(fa2.idcode))  idcode
FROM   ContactHit AS ch,
  Contact    AS c,
  FragInst   AS fia2,
  Frag       AS fa2,
  Seq        AS sa1,
  Seq        AS sb1
WHERE  c.id_frag_inst1 = ch.id_frag_inst_a2
AND    c.id_frag_inst2 = ch.id_frag_inst_b2
AND    c.isa_group IS FALSE
AND    c.crystal IS FALSE
AND    sa1.id = ch.id_seq_a1
AND    sb1.id = ch.id_seq_b1
AND    sa1.description NOT REGEXP '^Ig [[:blank:]]+ chain'
AND    sb1.description NOT REGEXP '^Ig [[:blank:]]+ chain'
AND    fia2.id = ch.id_frag_inst_a2
AND    fa2.id = fia2.id_frag
"
part3="
GROUP BY sa1.id, sb1.id
;
"

prey_bait=dbGetQuery(fistdb, paste0(part1,"AND    (sa1.primary_id, sb1.primary_id) IN (\n",pair_list,")",part3))
bait_prey=dbGetQuery(fistdb, paste0(part1,"AND    (sb1.primary_id, sa1.primary_id) IN (\n",pair_list,")",part3))

hits=data.table(rbind(prey_bait,bait_prey))
hits=rbind(hits,hits[,':='(idA=idB,acA=acB,idB=idA,acB=acA)])
setkey(hits,acA,acB)
hits=unique(hits)#[idA>=idB]

setkey(dt,ID_prey,ID_bait)

hits = hits[dt,nomatch=0]


write.table(
  hits, file.path('../Data',paste0("interactions mechismo knows about.tsv")), sep = "\t", row.names = F, qmethod =
    "double"
)

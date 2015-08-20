
Rcpp::sourceCpp('anchor.cpp')
cat(anchor("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"))
dt = setnames(fread(anchor("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")),
              c("res","res_kind","ANCHOR_prob","ANCHOR_binary","IUPred_prob","ANCHOR_score",
                "S","Eint","Egain"))


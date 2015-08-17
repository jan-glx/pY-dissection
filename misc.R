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
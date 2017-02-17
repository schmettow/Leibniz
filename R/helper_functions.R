rlogitnorm  <- function(n,m,s) plogis(rnorm(n,m,s))

## helper function
freq<-function(ms, occ = NULL){
  if(is.null(occ)) occ = unique(ms)
  sapply(occ, function(x) sum(ms == x))
}

# Logit-normal geometric distribution functions ####



#' logit-normal geometric distribution
#' 
#' distribution functions (probability, quantile and hazard) for the
#' logit-normal geometric distribution
#' 
#' The logit-normal geometric is created by imposing a logit-normal prior on
#' the parameter p of the geometric distribution.
#' 
#' @aliases dlngeom plngeom hlngeom
#' @usage dlngeom(x, m, s) plngeom(x, m, s) qlngeom(q, m, s) hlngeom(x, m, s)
#' @param x,q vector of random variable realizations or quantiles
#' @param m distribution parameter (central tendency)
#' @param s distribution parameter (dispersion)
#' @return probability (dlngeom), cumulative probabilty (qlngeom), quantile
#' (qlngeom) or hazard rate (hlngeom)
#' @note %% ~~further notes~~
#' @author Martin Schmettow
#' @seealso dlnbinom dlnbinom_zt dgeom
#' @references Schmettow, M. (2009). Controlling the usability evaluation
#' process under varying defect visibility. In BCS HCI 09: Proceedings of the
#' 23rd British HCI Group Annual Conference on People and Computers:
#' Celebrating People and Technology (pp. 188-197). Swinton, UK: British
#' Computer Society.
#' @examples
#' 
#'   dlngeom(3, 0.5, 2)
#'   qlngeom(3, 0.5, 2)
#'   
#' @export dlngeom
dlngeom<-function(k,m,s) {
  f<-function(P,k,m,s) dgeom(k,P)*dlogitnorm(P,m,s)
  F<-function(k,m,s) integrate(f,0,1,k,m,s)$value
  mapply(F,k,m,s)
}

plngeom<-function(k,m,s) {
  F<-function(k,m,s) 1-dlnbinom(0,k+1,m,s)
  mapply(F,k,m,s)
}

qlngeom<-function(q,m,s){
  f<-function(q,m,s){
    i<-0; 
    while(plngeom(i-1,m,s)<q) {i<-i+1};   
    return(i)
  }
  return(mapply(f,q,m,s))
}

## Hazard function ####

hlngeom<-function(k,m,s) dlngeom(k,m,s)/(1-plngeom(k,m,s))
hgeom<-function(k,p) dgeom(k,p)/(1-pgeom(k,p))


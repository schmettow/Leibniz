.dlnbinom <- function(x, size, m, s,
                      implementation = "pracma",
                      method = "Kronrod",
                      log = F){ ## not used
  x <- as.integer(x)
  size <- as.integer(size)
  fixed = choose(size,x)/(sqrt(2*pi*s^2))
  func = function(P) (1-P)^(size-x-1)*P^(x-1)*exp(-(logit(P)-m)^2/(2*s^2))

  if(x > size | x < 0)  return(0)
  if(s > 0) {
    if(implementation == "builtin")
      I <- integrate(func,
                     0,1,
                     rel.tol=.Machine$double.eps^0.25,
                     subdivisions=100)$value * fixed
    if(implementation == "pracma") {
      I <- pracma::integral(func, 0, 1) * fixed
    }

    return(I)
  }
  if(s == 0) return(dbinom(x, size, plogis(m)))
  if(s < 0) return(0)
}



#' logit-normal binomial distribution
#'
#' distribution functions (probability, random number generation and
#' cumulative) for the logit-normal binomial distribution
#'
#' @usage dlnbinom(x, size, m, s) rlnbinom(n, size, m, s) dlnbinom(p, size, m,
#' s)
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param m distribution parameter (central tendency)
#' @param s distribution parameter (dispersion)
#' @return dlnbinom gives the density, rlnbinom produces random numbers and
#' plnbinom gives the cumulative probability.
#'
#' The logit-normal distribution is created by imposing a logit-normal prior on
#' the binomial parameter p.
#'
#' @author Martin Schmettow
#' @seealso dlnbinom_zt dlngeom
#' @references Schmettow, M. (2009). Controlling the usability evaluation
#' process under varying defect visibility. In BCS HCI 09: Proceedings of the
#' 23rd British HCI Group Annual Conference on People and Computers:
#' Celebrating People and Technology (pp. 188-197). Swinton, UK: British
#' Computer Society.
#' @examples
#'
#'   dlnbinom(3, 10, 0.5, 2)
#'   rlnbinom(10, 10, 0.5, 2)
#'   plnbinom(0.5, 10, 0.5, 2)
#'
#' @export dlnbinom

dlnbinom <-
  Vectorize(.dlnbinom, vectorize.args = c("x","size" ,"m", "s"))



#' @rdname dlnbinom
#' @export

plnbinom <-
  function(x,size,m,s) sum(dlnbinom(c(0:x),size,m,s))

#' @rdname dlnbinom
#' @export
rlnbinom <-
  function(n,size,m,s) rbinom(n,size,rlogitnorm(n,m,s))



mean_lnbinom<-function(n,m,s) sum(dlnbinom(c(0:n),n,m,s)*c(0:n))


var_lnbinom <-
  function(n,m,s) {
  mean <- mean.lnbinom(n,m,s)
  mean((mean-c(0:n))^2*dlnbinom(c(0:n),n,m,s)*n)
}



## Zero-truncated LNB ####


#' zero-truncated logit-normal binomial distribution
#'
#' distribution functions (probability, random number generation and
#' cumulative) for the tzero-ztruncated and zero-one trauncated logit-normal
#' binomial distribution
#'
#' The logit-normal distribution is created by imposing a logit-normal prior on
#' the binomial parameter p. The zero-truncated LNB is derived by restricting
#' the random variable to a strictly positive range and re-scaling the
#' probability mass accordingly. The zero-one truncated distribution is created
#' respectively.
#'
#' @aliases dlnbinom_zot
#' @usage dlnbinom_zt(x, size, m, s) dlnbinom_zot(x, size, m, s)
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param m distribution parameter (central tendency)
#' @param s distribution parameter (dispersion)
#' @return dlnbinom_zt gives zero-truncated probability, dlnbinom_zot returns
#' the zero-one truncated probability
#' @note %% ~~further notes~~
#' @author Martin Schmettow
#' @seealso dlnbinom dlngeom
#' @references Schmettow, M. (2009). Controlling the usability evaluation
#' process under varying defect visibility. In BCS HCI 09: Proceedings of the
#' 23rd British HCI Group Annual Conference on People and Computers:
#' Celebrating People and Technology (pp. 188-197). Swinton, UK: British
#' Computer Society.
#' @examples
#'
#'   dlnbinom_zt(3, 10, 0.5, 2)
#'   dlnbinom_zot(0, 10, 0.5, 2)
#'
#' @export dlnbinom_zt
dlnbinom_zt<-function(x,size,m,s){
  d0<-(1-dlnbinom(0,size,m,s))
  (x>0)*dlnbinom(x,size,m,s)/d0
}

dlnbinom_zot<-function(x,size,m,s){
  d0<-(1-sum(dlnbinom(c(1,2),size,m,s)))
  (x>1)*dlnbinom(x,size,m,s)/d0
}

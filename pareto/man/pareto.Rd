\name{pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{p.dpareto}
\alias{p.ppareto}
\alias{p.qpareto}
\alias{rpareto}
\alias{rcpareto}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function, and pseudorandom number
  generator for the Pareto distribution with parameters 'scale' and 'shape'. 
  Functions with 'p.' allow for parallel processing. rcpareto uses C.
}
\usage{
  dpareto(x, a, b, log = FALSE)
  ppareto(q, a, b, lower.tail = TRUE, log.p = FALSE)
  qpareto(p, a, b, lower.tail = TRUE, log.p = FALSE)
  p.dpareto(x, a, b, log = FALSE, P = 1)
  p.ppareto(q, a, b, lower.tail = TRUE, log.p = FALSE, P = 1)
  p.qpareto(p, a, b, lower.tail = TRUE, log.p = FALSE, P = 1)
  rpareto(n, a, b)
  rcpareto(n, a, b)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{a, b}{scale and shape parameters. Must be positive, 'scale' strictly.}
  \item{n}{
    number of observations. If 'length(n) > 1', the length is taken to be the 
    number required.
  }
  \item{log, log.p}{
    logical; if 'TRUE', probabilities/densities p are returned as log(p).
  }
  \item{lower.tail}{
    logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].
  }
  \item{P}{specifies number of threads in MP functions.}
}
\details{
  The Pareto distribution with scale parameter a > 0 and shape
  parameter 'b' > 0 has density 
  
  f(x|a,b) = (b * a ^ b) / (x ^ (b + 1)) 
  
  at 'a' < 'x' and zero elsewhere.

}
\value{
  'dpareto' gives the density, 'ppareto' gives the distribution function,
  and 'qpareto' gives the quantile function. 
  
  Invalid arguments will give result in return value 'NaN', with a warning.
}
\author{Philip J. Erickson}
\examples{
  dpareto(3, 1, 2)
  p <- (1:5)/2
  dpareto(p, 1, 1, log = TRUE)
  ppareto(1:3, 1, 2, lower.tail = FALSE, log.p = TRUE)
  qpareto(0.5, 1:3, 2)
}
\keyword{distribution}
\keyword{univar}
\keyword{misc}

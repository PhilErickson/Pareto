dpareto <- function(x, a, b, log = FALSE) {
    lengthp <- max(length(x), length(a), length(b))
    result <- .C("paretodens",
                 as.double(x), as.integer(length(x)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)), 
                 as.integer(log), as.integer(lengthp),
                 val = double(lengthp), PACKAGE = "pareto")
	result$val
}

ppareto <- function(q, a, b, lower.tail = TRUE, log.p = FALSE) {
    lengthp <- max(length(q), length(a), length(b))
    result <- .C("paretodist",
                 as.double(q), as.integer(length(q)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)),
                 as.integer(lower.tail), as.integer(log.p), 
                 as.integer(lengthp), val = double(lengthp),
                 PACKAGE = "pareto")
    result$val
}

qpareto <- function(p, a, b, lower.tail = TRUE, log.p = FALSE) {
    lengthp <- max(length(p), length(a), length(b))
    result <- .C("paretoquant",
                 as.double(p), as.integer(length(p)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)),
                 as.integer(lower.tail), as.integer(log.p), 
                 as.integer(lengthp), val = double(lengthp),
                 PACKAGE = "pareto")
    result$val
}

p.dpareto <- function(x, a, b, log = FALSE, P = 1) {
    stopifnot(length(P) == 1)
    lengthp <- max(length(x), length(a), length(b))
    result <- .C("paretodensP",
                 as.double(x), as.integer(length(x)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)), 
                 as.integer(log), as.integer(lengthp),
                 val = double(lengthp), as.integer(P),
                 NAOK = TRUE, PACKAGE = "pareto")
	result$val
}

p.ppareto <- function(q, a, b, lower.tail = TRUE, log.p = FALSE, P = 1) {
    stopifnot(length(P) == 1)
    lengthp <- max(length(q), length(a), length(b))
    result <- .C("paretodistP",
                 as.double(q), as.integer(length(q)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)),
                 as.integer(lower.tail), as.integer(log.p), 
                 as.integer(lengthp), val = double(lengthp),
                 as.integer(P), NAOK = TRUE, PACKAGE = "pareto")
    result$val
}

p.qpareto <- function(p, a, b, lower.tail = TRUE, log.p = FALSE, P = 1) {
    stopifnot(length(P) == 1)
    lengthp <- max(length(p), length(a), length(b))
    result <- .C("paretoquantP",
                 as.double(p), as.integer(length(p)), 
                 as.double(a), as.integer(length(a)), 
                 as.double(b), as.integer(length(b)),
                 as.integer(lower.tail), as.integer(log.p), 
                 as.integer(lengthp), val = double(lengthp),
                 as.integer(P), NAOK = TRUE, PACKAGE = "pareto")
    result$val
}

rpareto <- function(n, a, b) {
    # Generate pseudorandom pareto-distributed vector of length n
    u <- runif(n)
    qpareto(u, a[1], b[1])
}

rcpareto <- function(n, a, b) {
    if(length(n) > 1) n <- length(n) # Just in case vector input
    result <- .C("paretorand", as.integer(n), 
                 as.double(a[1]), as.double(b[1]),
                 val = double(n), PACKAGE = "pareto")
    result$val
}
library(pareto)
stopifnot(tryCatch(is.na(dpareto(3,-2, 1)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(dpareto(3,2, -1)), error = function(e) TRUE))
stopifnot(all.equal(dpareto(3,2,1), 0.2222222222))
stopifnot(all.equal(dpareto(1,2,3), 0.0))
stopifnot(all.equal(dpareto(3:5,2, 1), c(0.2222222222, 0.1250000, 0.0800000)))
stopifnot(all.equal(dpareto(1:5,2, 1), c(0.0, 0.0, 0.2222222222, 0.1250000, 0.0800000)))
stopifnot(all.equal(dpareto(6,2:4, 1), c(0.05555555556, 0.08333333333, 0.11111111111)))
stopifnot(all.equal(log(dpareto(1:5,2, 1)), dpareto(1:5,2, 1, log = TRUE)))
stopifnot(all.equal(dpareto(6,1,2:4), c(0.0092592593, 0.0023148148, 0.0005144033)))

myppareto <- function(x, a, b, lower.tail = TRUE, log.p = FALSE) {
    a <- ifelse(a <= 0, NaN, a)
    b <- ifelse(b <= 0, NaN, b)
    x <- ifelse(x <= a, a, x)
    lp <- b * (log(a) - log(x))
    if (lower.tail) {
        p <- exp(lp)
        if (log.p) log1p(-p)
        else 1 - p
    }
    else {
        if (log.p) lp
        else exp(lp)
    }
}

myqpareto <- function(p, a, b, lower.tail = TRUE, log.p = FALSE) {
    if (log.p) p <- exp(p)
    if (lower.tail) p <- 1 - p
    a <- ifelse(a <= 0, NaN, a)
    b <- ifelse(b <= 0, NaN, b)
    p <- ifelse(p < 0 | p > 1, NaN, p)
    a / p^(1/b)
}
       
stopifnot(tryCatch(is.na(ppareto(3,-2, 1)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(ppareto(3,2, -1)), error = function(e) TRUE))
stopifnot(all.equal(ppareto(3:5,2, 1), myppareto(3:5,2, 1)))
stopifnot(all.equal(ppareto(1:5,2, 1), myppareto(1:5,2, 1)))
stopifnot(all.equal(ppareto(6,2:4, 1), myppareto(6,2:4, 1)))
stopifnot(all.equal(ppareto(6,1,2:4), myppareto(6,1,2:4)))
stopifnot(all.equal(ppareto(6,1,2:4, log.p = TRUE),
                    myppareto(6,1,2:4, log.p = TRUE)))
stopifnot(all.equal(ppareto(6,1,2:4, log.p = TRUE, lower.tail = FALSE),
                    myppareto(6,1,2:4, log.p = TRUE, lower.tail = FALSE)))
stopifnot(all.equal(ppareto(10, 1, 20, log.p = TRUE, lower.tail = FALSE),
                    myppareto(10, 1, 20, log.p = TRUE, lower.tail = FALSE)))

stopifnot(tryCatch(is.na(qpareto(0.5,-2, 1)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(qpareto(0.5,2, -1)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(qpareto(-1,2, 1)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(qpareto(2,2, 1)), error = function(e) TRUE))
stopifnot(all.equal(qpareto((3:5)/6,2, 1), myqpareto((3:5)/6,2, 1)))
stopifnot(all.equal(qpareto((1:5)/6,2, 1), myqpareto((1:5)/6,2, 1)))
stopifnot(all.equal(qpareto(0.25,2:4, 1), myqpareto(0.25,2:4, 1)))
stopifnot(all.equal(qpareto(0.25,1,2:4), myqpareto(0.25,1,2:4)))
stopifnot(all.equal(qpareto(log(0.25),1,2:4, log.p = TRUE),
                    myqpareto(log(0.25),1,2:4, log.p = TRUE)))
stopifnot(all.equal(qpareto(0.25,1,2:4, lower.tail = FALSE),
                    myqpareto(0.25,1,2:4, lower.tail = FALSE)))
stopifnot(all.equal(qpareto(log(0.25),1,2:4, log.p = TRUE,
                            lower.tail = FALSE),
                    myqpareto(log(0.25),1,2:4, log.p = TRUE,
                              lower.tail = FALSE)))
stopifnot(all.equal(ppareto(2:4, 1:3, 5:7), myppareto(2:4, 1:3, 5:7)))

# MP

stopifnot(tryCatch(is.na(p.dpareto(3,-2, 1, P = 2)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.dpareto(3,2, -1, P = 2)), error = function(e) TRUE))
stopifnot(all.equal(p.dpareto(3,2,1, P = 2), 0.2222222222))
stopifnot(all.equal(p.dpareto(1,2,3, P = 2), 0.0))
stopifnot(all.equal(p.dpareto(3:5,2, 1, P = 2), c(0.2222222222, 0.1250000, 0.0800000)))
stopifnot(all.equal(p.dpareto(1:5,2, 1, P = 2), c(0.0, 0.0, 0.2222222222, 0.1250000, 0.0800000)))
stopifnot(all.equal(p.dpareto(6,2:4, 1, P = 2), c(0.05555555556, 0.08333333333, 0.11111111111)))
stopifnot(all.equal(log(p.dpareto(1:5,2, 1, P = 2)), p.dpareto(1:5,2, 1, log = TRUE, P = 2)))
stopifnot(all.equal(p.dpareto(6,1,2:4, P = 2), c(0.0092592593, 0.0023148148, 0.0005144033)))

stopifnot(tryCatch(is.na(p.ppareto(3,-2, 1, P = 2)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.ppareto(3,2, -1, P = 2)), error = function(e) TRUE))
stopifnot(all.equal(p.ppareto(3:5,2, 1, P = 2), myppareto(3:5,2, 1)))
stopifnot(all.equal(p.ppareto(1:5,2, 1, P = 2), myppareto(1:5,2, 1)))
stopifnot(all.equal(p.ppareto(6,2:4, 1, P = 2), myppareto(6,2:4, 1)))
stopifnot(all.equal(p.ppareto(6,1,2:4, P = 2), myppareto(6,1,2:4)))
stopifnot(all.equal(p.ppareto(6,1,2:4, log.p = TRUE, P = 2),
                    myppareto(6,1,2:4, log.p = TRUE)))
stopifnot(all.equal(p.ppareto(6,1,2:4, log.p = TRUE, lower.tail = FALSE, P = 2),
                    myppareto(6,1,2:4, log.p = TRUE, lower.tail = FALSE)))
stopifnot(all.equal(p.ppareto(10, 1, 20, log.p = TRUE, lower.tail = FALSE, P = 2),
                    myppareto(10, 1, 20, log.p = TRUE, lower.tail = FALSE)))

stopifnot(tryCatch(is.na(p.qpareto(0.5,-2, 1, P = 2)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.qpareto(0.5,2, -1, P = 2)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.qpareto(-1,2, 1, P = 2)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.qpareto(2,2, 1, P = 2)), error = function(e) TRUE))
stopifnot(all.equal(p.qpareto((3:5)/6,2, 1, P = 2), myqpareto((3:5)/6,2, 1)))
stopifnot(all.equal(p.qpareto((1:5)/6,2, 1, P = 2), myqpareto((1:5)/6,2, 1)))
stopifnot(all.equal(p.qpareto(0.25,2:4, 1, P = 2), myqpareto(0.25,2:4, 1)))
stopifnot(all.equal(p.qpareto(0.25,1,2:4, P = 2), myqpareto(0.25,1,2:4)))
stopifnot(all.equal(p.qpareto(log(0.25),1,2:4, log.p = TRUE, P = 2),
                    myqpareto(log(0.25),1,2:4, log.p = TRUE)))
stopifnot(all.equal(p.qpareto(0.25,1,2:4, lower.tail = FALSE, P = 2),
                    myqpareto(0.25,1,2:4, lower.tail = FALSE)))
stopifnot(all.equal(p.qpareto(log(0.25),1,2:4, log.p = TRUE,
                            lower.tail = FALSE, P = 2),
                    myqpareto(log(0.25),1,2:4, log.p = TRUE,
                              lower.tail = FALSE)))
stopifnot(all.equal(p.ppareto(2:4, 1:3, 5:7, P = 2), myppareto(2:4, 1:3, 5:7)))

stopifnot(all.equal(p.ppareto(2:4, 1:3, 5:7, P = 2), p.ppareto(2:4, 1:3, 5:7, P = 3)))
stopifnot(all.equal(p.qpareto(0.25,1,2:4, lower.tail = FALSE, P = 2),
                    p.qpareto(0.25,1,2:4, lower.tail = FALSE, P = 3)))
stopifnot(all.equal(p.dpareto(1:5,2, 1, P = 2), p.dpareto(1:5,2, 1, P = 4)))

stopifnot(tryCatch(is.na(p.dpareto(3,2,1, P = 2:3)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.ppareto(3:5,2, 1, P = 2:4)), error = function(e) TRUE))
stopifnot(tryCatch(is.na(p.qpareto((3:5)/6,2, 1, P = 2:3)), error = function(e) TRUE))

set.seed(3212333); a <- rpareto(2, 1, 1)
set.seed(3212333); b <- rcpareto(2, 1, 1)
stopifnot(all.equal(a, b))
set.seed(3212333); a <- rpareto(2, 1:3, 1)
set.seed(3212333); b <- rcpareto(2, 1:2, 1)
stopifnot(all.equal(a, b))
set.seed(3212333); a <- rpareto(2, 1, 1:2)
set.seed(3212333); b <- rcpareto(2, 1, 1:3)
stopifnot(all.equal(a, b))
set.seed(3212333); a <- rpareto(2:3, 1, 1)
set.seed(3212333); b <- rcpareto(3:4, 1, 1)
stopifnot(all.equal(a, b))




library("truncnorm")
library("FNN")

GenRandom <- function(m, s, xlo, xhi, ylo, yhi) {
    dmin <- rtruncnorm(1, 0, mean=m, sd=s)
    x <- runif(1, xlo, xhi)
    y <- runif(1, ylo, yhi)

    return(list(x=x, y=y, dmin=dmin))
}

ValidPoint <- function(rp, points, ngen, dmin) {
    d <- 0
    if (ngen == 0) {
        return(TRUE)
    }

    if (length(points[1:ngen, ]) == 2) {
        data <- matrix(points[1:ngen, ], ncol=2)
        d <- min(knnx.dist(data, rp, k=ngen))
    } else {
        d <- min(knnx.dist(points[1:ngen, ], rp, k=ngen))
    }

    return(d > dmin)
}

dmin2d <- function(n, m, s, xlo, xhi, ylo, yhi) {
    ## Input
    ## n: number of points to simulate
    ## m: mean of Normal distribution
    ## s: s.d. of Normal distribution
    ## xlo, xhi: possible range of X values.
    ## ylo, yhi: possible range of Y values.
    ## Output
    ## Returns a n by 2 matrix

    kMaxTries <- 10
    #matrix to store results
    points <- matrix(data=rep(0, n*2), nrow=n)
    ngen <- 0
    n.tries <- 0

    while(ngen < n &
          n.tries < kMaxTries) {
        r <- GenRandom(m, s, xlo, xhi, ylo, yhi)
        rp <- matrix(c(r$x, r$y), ncol=2, byrow=T)
        if (ValidPoint(rp, points, ngen, r$dmin)) {
            n.tries <- 0
            ngen <- ngen + 1
            points[ngen, ] <- rp
        } else {
            n.tries <- n.tries + 1
        }
    }

    return(points)
}

plotDmin2d <- function(points, ...) {

    plot(points[, 1], points[, 2], ...)
}





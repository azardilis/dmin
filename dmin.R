library("truncnorm")
library("FNN")

genRandom <- function(m, s, xlo, xhi, ylo, yhi) {
    if (m == 0) { dmin = 0 }
    else {dmin <- rtruncnorm(1, 0, mean=m, sd=s)}
    x <- runif(1, xlo, xhi)
    y <- runif(1, ylo, yhi)

    return(list(x=x, y=y, dmin=dmin))
}

validPoint <- function(rp, points, ngen, dmin) {
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

    return(d >= dmin)
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

    kMaxTries <- 100
    #matrix to store results
    points <- matrix(data=rep(0, n*2), nrow=n)
    ngen <- 0
    n.tries <- 0

    while(ngen < n &
          n.tries < kMaxTries) {
        r <- genRandom(m, s, xlo, xhi, ylo, yhi)
        rp <- matrix(c(r$x, r$y), ncol=2, byrow=T)
        if (validPoint(rp, points, ngen, r$dmin)) {
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

ri <- function(points) {
    dists <- knn.dist(points, 1)
    return(mean(dists) / sd(dists))
}

riTest <- function(n) {
    rep.res <- replicate(100, dmin2d(n = n, m = 0, s = 0,
                                     xlo = 0, xhi = 1000,
                                     ylo = 0, yhi = 1000),
                         simplify=FALSE)

    ri.res <- sapply(rep.res, ri)
    return(sort(ri.res, decreasing=T)[50])
}

scoreParams <- function(p) {
    n <- 242
    real.ri <- 4.753291
    m <- p[1]
    s <- p[2]

    rep.res <- replicate(10, dmin2d(n = n, m = m, s = s,
                                     xlo = 0, xhi = 400,
                                     ylo = 0, yhi = 400),
                         simplify=FALSE)
    ri.res <- sapply(rep.res, ri)
    avg.ri <- sum(ri.res) / (length(ri.res) - 1)

    return(abs(real.ri - avg.ri))
}

calcRealRi <- function(fname="spa2.dat") {
    points <- as.matrix(read.table(fname))
    return(ri(points))
}

mcmc <- function() {
    sigma <- 5
    T <- 10
    chain <- matrix(rep(0, T*2), ncol=2)
    chain[1, ] <- c(10, 5)

    for(t in 2:T) {
        print(t)
        ym <- rtruncnorm(1, 0, mean=chain[t-1, 1], sd=sigma)
        ys <- rtruncnorm(1, 0, mean=chain[t-1, 2], sd=sigma)
        cat(sprintf("%f %f \n", ym, ys))
        if (runif(1) > (scoreParams(c(ym, ys)) / scoreParams(chain[t-1, ]))) {
            chain[t, ] <- c(ym, ys)
        } else {
            chain[t, ] <- chain[t-1, ]
        }
    }

    return(chain)
}







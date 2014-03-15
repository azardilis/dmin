library("truncnorm")
library("FNN")


## ---- dmin2d ----
genRandom <- function(m, s, xlo, xhi, ylo, yhi) {
    dmin <- rtruncnorm(1, 0, mean=m, sd=s)
    x <- runif(1, xlo, xhi)
    y <- runif(1, ylo, yhi)

    return(list(x=x, y=y, dmin=dmin))
}

validPoint <- function(rp, points, ngen, dmin) {
    d <- 0
    if (ngen == 0) return(TRUE)

    if (length(points[1:ngen, ]) == 2) {
        data <- matrix(points[1:ngen, ], ncol=2)
        d <- min(knnx.dist(data, rp, k=ngen))
    } else {
        d <- min(knnx.dist(points[1:ngen, ], rp, k=ngen))
    }

    return(d >= dmin)
}

dmin2d <- function(n, m, s, xlo, xhi, ylo, yhi, mt=2000) {
    ## Input
    ## n: number of points to simulate
    ## m: mean of Normal distribution
    ## s: s.d. of Normal distribution
    ## xlo, xhi: possible range of X values.
    ## ylo, yhi: possible range of Y values.
    ## Output
    ## Returns a n by 2 matrix

    #special case where m=0
    if (m == 0) return(cbind(runif(n, xlo, xhi), runif(n, ylo, yhi)))
    kMaxTries <- mt
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

## ---- plotDmin2d ----
plotDmin2d <- function(points, xlo, xhi, ylo, yhi, ...) {
    plot(points, pch=16, xlab="x", ylab="y", ...)
    rect(xlo, ylo, xhi, yhi, lty=2)
}

## ---- ri ----
ri <- function(points) {
    dists <- knn.dist(points, 1)
    return(mean(dists) / sd(dists))
}

## ---- riTest ----
riTest <- function(n, xhi, yhi) {
    print(n)
    rep.res <- replicate(1000, dmin2d(n = n, m = 0, s = 0,
                                      xlo = 0, xhi = xhi,
                                      ylo = 0, yhi = yhi),
                         simplify=FALSE)
    ri.res <- sapply(rep.res, ri)
    return(ri.res)
    #return(sort(ri.res, decreasing=T)[50])
}

## ---- scoreParams ----
scoreParamsPVal <- function(p, n, real.ri) {
    m <- p[1]
    s <- p[2]

    rep.res <- replicate(99, dmin2d(n = n, m = m, s = s,
                                    xlo = 0, xhi = 400,
                                    ylo = 0, yhi = 400),
                         simplify=FALSE)
    ri.res <- c(real.ri, sapply(rep.res, ri))
    uis <- sapply(1:length(ri.res), function(x) { abs(ri.res[x] - mean(ri.res[-x])) })

    pval <- match(1, sort(uis, decreasing=T, index.return=T)$ix) / length(uis)
    return(pval)
}

scoreParams <- function(p, n, real.ri) {
    m <- p[1]
    s <- p[2]

    rep.res <- replicate(99, dmin2d(n = n, m = m, s = s,
                                    xlo = 0, xhi = 400,
                                    ylo = 0, yhi = 400),
                         simplify=FALSE)
    avg.ri <- mean(sapply(rep.res, ri))

    return(abs(real.ri - avg.ri))
}

## ---- calcRealRi ----
calcRealRi <- function(fname="spa2.dat") {
    points <- as.matrix(read.table(fname))
    return(ri(points))
}

## ---- mcmc ----
mcmc <- function(n, real.ri) {
    sigma <- 5
    T <- 200
    chain <- matrix(rep(0, T*2), ncol=2)
    chain[1, ] <- c(5, 5)
    prev.score <- scoreParams(chain[1, ], n, real.ri)
    cat(sprintf("%f %f %f\n", chain[1,1], chain[1,2], prev.score))
    for(t in 2:T) {
        print(t)
        ym <- rtruncnorm(1, 0, mean=chain[t-1, 1], sd=sigma)
        ys <- rtruncnorm(1, 0, mean=chain[t-1, 2], sd=sigma)
        cand.score <- scoreParams(c(ym, ys), n, real.ri)
        cat(sprintf("%f %f %f\n", ym, ys, cand.score))
        if (runif(1) > (cand.score / prev.score)) {
            chain[t, ] <- c(ym, ys)
            prev.score <- cand.score
        } else {
            chain[t, ] <- chain[t-1, ]
        }
    }

    return(chain)
}

## ---- failModel ----
failModel <- function(points) {
    return(any(apply(points, 1, function(x) {any((x == 0) == FALSE)}) == FALSE))
}

## ---- packDensity ----
packDensity <- function() {
    tries <- seq(from=0, to=1000, by=100)
    res <- rep(0, length(tries))
    names(res) <- as.character(tries)
    for (i in 1:length(tries)) {
        points <- dmin2d(n=400, m = 20, s = 0, xlo = 0, xhi = 400, ylo = 0,
                         yhi = 400, mt=tries[i])
        res[i] <- length(which(apply(points, 1, function(x) { any(x == 0)}) == FALSE))
    }

    return(res)
}

packDensityN <- function() {
    n <- seq(10, 400, 10)
    res <- rep(0, length(n))
    for (i in 1:length(n)) {
        print(n[i])
        res[i] <- failModel(dmin2d(n=n[i], 20, 0, 0, 400, 0, 400))
    }
    names(res) <- n

    return(tail(names(which(res == 0)), n=1))
}


testArea <- function(area) {
    side1 <- c(1, 10, 50, 100, 500, 1000, 2500, 5000, 10000, 50000, 100000)
    side2 <- area / side1

    res <- mapply(riTest, rep(200, length(side1)), side1, side2, SIMPLIFY=FALSE)
    ri.m <- sapply(res, mean)
    ri.sd <- sapply(res, sd)
    ri.50 <- sapply(res, function(x) {sort(x, decreasing=T)[50]})

    return(list(side=side1, ri.m=ri.m, ri.sd=ri.sd, ri.50=ri.50))
}

plotRIArea <- function(res.ri) {
    x <- 1:length(res.ri$ri.m)
    plot(x, res.ri$ri.m, pch=19, xlab="side 1", ylab="mean RI", ylim=c(0.9, 2.0), xaxt='n')
    arrows(x, res.ri$ri.m-res.ri$ri.sd, x, res.ri$ri.m+res.ri$ri.sd,
           length=0.05, angle=90, code=3)

    axis(1, at=1:11, labels=as.character(res.ri$side))
}




dens <- function(d, n , xlo, xhi, ylo, yhi) {
    ta <- (xhi - xlo) * (yhi - ylo)
    da <- pi * d^2

    return((1 - (da / ta))^(n-1) * 2*pi*d*n/ta)
}

exp.d <- function(d, pd) {
    return(sum(d * pd))
}

var.d <- function(d, pd, m) {
    return(sqrt(sum(pd * (d - m)^2)))

}

n <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
d <- 1:100
res <- sapply(n, function(x) {dens(d, n=x, 0, 1000, 0, 1000)}, simplify=FALSE)
exp.rep <- sapply(res, function(x) {exp.d(d, x)})
i <- 1:length(res)
var.rep <- sapply(i, function(x) {var.d(d, res[[x]], exp.rep[x])})

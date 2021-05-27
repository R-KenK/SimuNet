# Loading used packages --------------------------------------
library(ggplot2)
library(ggridges)
library(cowplot)
library(data.table)
library(dplyr)
library(extraDistr)

# Draw P0 ---------------------------------------------------------------
## As a function of a gregarious index----
draw_greg.index <- function(n,rng.vec = rnbinom(n,mu = 2,size = .7),min = 0.025,max = 0.975,node_names = NULL) {
  names(rng.vec) <- node_names
  scales::rescale(rng.vec,to = c(min,max))
}

shape_dyads.df <- function(greg.index,n = length(greg.index)) {
  greg.prod <- greg.index %o% greg.index
  diag(greg.prod) <- 0
  expand.grid(i = 1:n,j = 1:n) %>%
    {.$greg.prod <- greg.prod[cbind(.$i,.$j)]; . }
}

draw_dyadic.noise <- function(dyads.df,noise.min = 0.025,mean = 0, sd = 0.25) {
  truncnorm::rtruncnorm(n = nrow(dyads.df),
                        a = -dyads.df$greg.prod + noise.min,
                        b = 1 - dyads.df$greg.prod - noise.min,
                        mean = mean,
                        sd = sd
  )
}

add_dyadic.noise <- function(dyads.df,noise.min = 0.025,mean = 0, sd = 0.25) {
  dyadic.noise <- draw_dyadic.noise(dyads.df,noise.min,mean,sd)
  dyads.df$dyadic.noise <- ifelse(dyads.df$i != dyads.df$j,dyadic.noise,0)
  dyads.df$p.ij <- dyads.df$greg.prod + dyads.df$dyadic.noise
  dyads.df
}

draw_P0_greg <- function(n,rng.vec = rnbinom(n,mu = 2,size = .7),min = 0.025,max = 0.975,noise.min = 0.025,noise.mean = 0,noise.sd = 0.25,node_names = NULL) {
  n %>%
    draw_greg.index(rng.vec = rng.vec,min = min,max = max,node_names = node_names) %>%
    shape_dyads.df() %>%
    add_dyadic.noise(noise.min = noise.min,mean = noise.mean,sd = noise.sd) %>%
    {matrix(.$p.ij,nrow = n,ncol = n,dimnames = list(node_names,node_names))}
}

## As an evenly spaced sequence of probabilities ----
draw_P0_seq <- function(n,min = 0,max = 1,node_names = NULL) {
  P0 <- matrix(0,ncol = n,nrow = n,dimnames = list(node_names,node_names))
  SimuNet::non.diagonal(P0) <- seq(0,1,length.out = n * (n - 1))
  P0
}

## Wrapper ----
draw_P0 <- function(n,method = c("sequence","gregarious.index"),
                    rng.vec = rnbinom(n,mu = 2,size = .7),
                    min = 0.025,
                    max = 0.975,
                    noise.min = 0.025,
                    noise.mean = 0,
                    noise.sd = 0.25,
                    node_names = NULL) {
  method <- match.arg(method)
  switch(method,
         "sequence" = draw_P0_seq(n = n,min = min,max = max,node_names = node_names),
         "gregarious.index" = draw_P0_greg(n = n,
                                           rng.vec = rng.vec,
                                           min = min,
                                           max = max,
                                           noise.min = noise.min,
                                           noise.mean = noise.mean,
                                           noise.sd = noise.sd,
                                           node_names = node_names
         )
  )
}

# Draw P, scans and adjacency matrices ---------------------------------------

## Drawing A0 ----

### Adjacency from binomial probabilities ----
Mat_rbinom <- function(N,p) {
  Adj <- rbinom(length(p),N,p)
  if (is.array(p)) {
    dim(Adj) <- dim(p)
    attr(Adj,"dimnames") <- attr(p,"dimnames")
  }
  Adj
}

## Drawing A.sim ----
### A.sim from BetaBinomial ----
Mat_rbbinom <- function(A.obs,N,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  Adj <-
    extraDistr::rbbinom(
      n = length(A.obs),
      size = N.new,
      alpha = A.obs + alpha.prior,
      beta = N - A.obs + beta.prior
    )
  if (is.array(A.obs)) {
    dim(Adj) <- dim(A.obs)
    attr(Adj,"dimnames") <- attr(A.obs,"dimnames")
    diag(Adj) <- 0
  }
  Adj
}

## Binomial probability from Beta distribution ----
Mat_rbeta <- function(A.obs,N,alpha.prior = 0.5,beta.prior = 0.5) {
  P <- rbeta(n = length(A.obs),A.obs + alpha.prior,N - A.obs + beta.prior)
  if (is.array(A.obs)) {
    dim(P) <- dim(A.obs)
    attr(P,"dimnames") <- attr(A.obs,"dimnames")
    diag(P) <- 0
  }
  P
}

## Draw binary scan list ----
draw_scanList <- function(P = NULL,A.obs = NULL,N,n.scans = N) {
  if (is.null(P) & !is.null(A.obs)) {
    P <- Mat_rbeta(A.obs,N)
  }
  scan.list <-
    list(
      scans = replicate(
        n = n.scans,
        simplify = FALSE,
        Mat_rbinom(P,N = 1)
      ),
      P = P
    )
  class(scan.list) <- "sL"
  scan.list
}

### related methods ----
use_printSpMatrix <- function(M,...,sparse = TRUE, digits = 3, note.dropping.colnames = FALSE,align = "right") {
  M <- Matrix::Matrix(M,sparse = sparse)
  Matrix::printSpMatrix(M,digits = digits,note.dropping.colnames = note.dropping.colnames,align = align,...)
}

print.sL <- function(x,...) {
  n <- length(x$scans)
  if (n >= 10) {
    lapply(1:3,function(i) print_scan_i(x,i))
    cat("\n... (",n - 3," more scans)\n\n",sep = "")
    lapply((n - 2 + 1):n,function(i) print_scan_i(x,i))
  } else {
    lapply(1:n,function(i) print_scan_i(x,i))
  }
}

print_scan_i <- function(x,i) {
  cat(paste0("[[",i,"]]\n"))
  use_printSpMatrix(x$scans[[i]])
}

get_sL <- function(scan.list,attribute) {
  scan.list[[attribute]]
}

## Sum binary scan list into a weighted adjacency matrix ----
sum_scanList <- function(scan.list,attribute = "scans") {
  get_sL(scan.list,attribute) %>%
    Reduce("+",.)
}

# Uncertainty assessment tools -------------------------------------------
# Scan list Bootstrap ----
resample_scanList <- function(scan.list) {
  scan.list$scans <-
    scan.list$scans %>%
    sample(replace = TRUE)
  class(scan.list) <-  "sL"
  scan.list
}

bootstrap_scanList <- function(scan.list,n.boot = 100,attribute = "scans",simplify = FALSE) {
  replicate(n = n.boot,
            simplify = simplify,
            resample_scanList(scan.list) %>%
              sum_scanList(attribute = attribute)
  )
}


# Data extraction tools ---------------------------------------------------

extract_xij <- function(M,N = NULL,x.name = "a",...) {
  n <- nrow(M)
  expand.grid(i = 1:n,j = 1:n) %>% {
    .$N <- N
    .[[x.name]] <- M[cbind(.$i,.$j)]
    .
  } %>%
    cbind(...,.) %>%
    subset(i != j)
}

format_xij.df <- function(A.list,N = NULL,...) {
  n <- nrow(A.list[[1]])
  cbind(
    ...,
    rep = rep(1:length(A.list),each = n^2 - n),
    rbind_lapply(A.list,extract_xij,N = N)
  )
}

merge_with_A0 <- function(xij.df,A0.df) {
  merge(xij.df,A0.df,by = c("i","j","N")) %>%
    arrange(N,rep,j,i) %>%
    relocate(N,n.a0,method,rep,i,j,a0,a)
}

format_and_merge_xij.df <- function(A.list,N,method,A0.df) {
  A.list %>%
    format_xij.df(N = N,method = method) %>%
    merge_with_A0(A0.df)
}

# other tools -------------------------------------------------------------

rbind_lapply <- function(X,FUN,...){
  do.call(rbind,lapply(X = X,FUN = FUN,...))
}

non.diagonal<- function(M,output=c("matrix.logical","vector.values")) {
  output<- match.arg(output)
  if(dim(M)[1]==dim(M)[2]) logicals<- upper.tri(M,diag = FALSE)|lower.tri(M,diag = FALSE) else stop("Matrix provided is not a square matrix.")
  switch(output,
         "matrix.logical" = logicals,
         "vector.values" = M[logicals])
}

`non.diagonal<-` <- function(x,value) {
  dx <- dim(x)
  if (length(dx) != 2L)
    stop("only matrix non-diagonals can be replaced")
  len.i <- prod(dx) - min(dx)
  len.v <- length(value)
  if (len.v != 1L && len.v != len.i)
    stop("replacement non-diagonal has wrong length")
  x[non.diagonal(x)] <- value
  x
}

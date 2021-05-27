# required packages -------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(extraDistr)
devtools::load_all()

# network generating tools ------------------------------------------------

## generating P -----------------------------------------------------------

generate_P.seq <- function(n,node_names = as.character(1:n)) {
  M <- matrix(0,n,n,dimnames = list(node_names,node_names))
  SimuNet::non.diagonal(M) <-
    seq(0,1,length = n^2 - n)
  M
}

generate_P.soc <- function(n,node_names = as.character(1:n)) {
  M <- matrix(0,n,n,dimnames = list(node_names,node_names))
  SimuNet::non.diagonal(M) <-
    seq(0,1,length = n^2 - n)
  M
}

generate_P.wrapper <- function(n,P.method = c("sequence","social"),node_names = as.character(1:n),...) {
  switch(P.method,
         "sequence" = generate_P.seq(n = n, node_names = node_names),
         "social" = generate_P.soc(n = n, node_names = node_names,...)
  )
}


## inferring P from observations ---------------------------------------
mat_rbeta <- function(A,N,alpha.prior = 0.5,beta.prior = 0.5) {
  P <- rbeta(n = length(A),A + alpha.prior,N - A + beta.prior)
  if (is.array(A)) {
    dim(P) <- dim(A)
    attr(P,"dimnames") <- attr(A,"dimnames")
    diag(P) <- 0
  }
  P
}

## draw scan(s) from P ----------------------------------------------------
mat_rbinom <- function(P,N) {
  Adj <- rbinom(length(P),N,P)
  if (is.array(P)) {
    dim(Adj) <- dim(P)
    attr(Adj,"dimnames") <- attr(P,"dimnames")
  }
  Adj
}

draw_scanList <- function(P = NULL,A = NULL,N,n.scans = N) {
  if (is.null(P) & !is.null(A)) {
    P <- mat_rbeta(A,N)
  }
  scan.list <- replicate(n = n.scans,simplify = FALSE,mat_rbinom(P,N = 1))
  class(scan.list) <- "sL"
  scan.list
}

print.sL <- function(x,...) {
  N <- length(x)
  if (N >= 10) {
    lapply(x[1:3],print)
    cat("\n... (",N - 3," more scans)\n\n",sep = "")
    lapply(x[(N - 2 + 1):N],print)
  } else {
    print.default(x)
  }
}


## draw adjacency from Beta-binomial  -------------------------------------
mat_rbbinom <- function(A,N,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  Adj <-
    extraDistr::rbbinom(
      n = length(A),
      size = N.new,
      alpha = A + alpha.prior,
      beta = N - A + beta.prior
    )
  if (is.array(A)) {
    dim(Adj) <- dim(A)
    attr(Adj,"dimnames") <- attr(A,"dimnames")
    diag(Adj) <- 0
  }
  Adj
}

## sum scans --------------------------------------------------------------
sum_scanList <- function(scan.list) {
  Reduce("+",scan.list)
}


## resample scanList -------------------------------------------------------
resample_scanList <- function(scan.list) {
  scan.list <- sample(scan.list,replace = TRUE)
  class(scan.list) <-  "sL"
  scan.list
}

bootstrap_scanList <- function(scan.list,n.boot = 100,simplify = FALSE) {
  replicate(n = n.boot,
            simplify = simplify,
            resample_scanList(scan.list) %>%
              sum_scanList()
  )
}


## wrapper for generating through each method -----------------------------
generate_infered_networks <- function(P,N,methods = c("bbinom","SimuNet","boot"),N.new = N,n.rep = 100,n.samp = 1000) {
  sL0 <- draw_scanList(P = P,n.scans = N)
  A0 <- sum_scan.list(sL0)


  list(
    A.bbinom = ,
    A.sim = ,
    A.boot =
    )
}


# data extraction tools ---------------------------------------------------

## extract xij ------------------------------------------------------------
extract_xij <- function(M,N = NULL,x.name = "x",...,keep.diag = FALSE) {
  n <- nrow(M)
  expand.grid(i = 1:n,j = 1:n) %>% {
    .$N <- N
    .[[x.name]] <- M[cbind(.$i,.$j)]
    .
  } %>%
    cbind(...,.) %>%
    subset(i != j | keep.diag)
}


# calculate social network metrics ----------------------------------------

## whole network metric ---------------------------------------------------
calculate_SNm <- function(Adj,SNm_fun,...,Xapply = lapply) {
  Xapply(Adj,SNm_fun,...)
}

## node network metric ---------------------------------------------------
calculate_SNm <- function(Adj,SNm_fun,transform.to.igraph = FALSE,mode = NULL,weighted = TRUE,...,Xapply = lapply) {
  if (transform.to.igraph) {
    Adj <- lapply(Adj,igraph::graph.adjacency,mode = mode,weighted = weighted)
  }
  Xapply(Adj,SNm_fun,...)
}

## wrapper for node metrics ------------------------------------------------



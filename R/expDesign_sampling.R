# empirical sampling methods ----

## Wrapper to create sampling functions ----

#'  TO WRITE
#'
#' @param sampling TO WRITE
#' @param all.sampled TO WRITE
#' @param method TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
customize_sampling <- function(method = c("group","focal"),
                               sampling =  c("constant","matrix","even","random","function"),
                               all.sampled = TRUE) {

  method <- match.arg(method)
  check_sampling_parameters(method = method,sampling = sampling)

  switch(method,
         "group" = \(scan.list) group_sample(scan.list = scan.list,
                                             sampling = sampling,
                                             all.sampled = all.sampled
         ),
         "focal" = \(scan.list) focal_sample(scan.list = scan.list,
                                             sampling = sampling,
                                             all.sampled = all.sampled
         )
  )
}

#' TO WRITE
#'
#' @param method TO WRITE
#' @param sampling TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @keywords internal
check_sampling_parameters <-
  function(method,
           sampling =  c("constant","matrix","even","random","function")) {
    switch(
      method,
      "group" =
        switch(
          class(sampling)[1],
          "character" = {
            sampling <- match.arg(sampling)
            switch(
              sampling,
              "even" = stop('Invalid `sampling` option for method = "group" (cf. `?customize_sampling`)'),
              "random" =,
              "constant" =, # this is the default choice when no sampling argument is inputted
              "matrix" =,
              "function" = {}
            )
          },
          "function" =,"numeric" =,"matrix" = {}
        ),
      "focal" =
        switch(
          class(sampling)[1],
          "character" = {
            sampling <- match.arg(sampling)
            switch(
              sampling,
              "constant" =, # this is the default choice when no sampling argument is inputted
              "matrix" = stop('Invalid `sampling` option for method = "focal" (cf. `?customize_sampling`)'),
              "random" =,
              "even" =,
              "function" = {}
            )
          },
          "matrix" =,
          "numeric" = stop('Invalid `sampling` option for method = "focal" (cf. `?customize_sampling`)'),
          "function" = {}
        )
    )
  }


## Group-scan sampling ----

#'  TO WRITE
#'
#' @param sampling TO WRITE
#' @param all.sampled TO WRITE
#' @param scan.list TO WRITE
#'
#' @return  TO WRITE
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#' # TO WRITE
group_sample <- function(scan.list,sampling = c("constant","matrix","random","function"),
                         all.sampled = TRUE){

  Adj <- attrs(scan.list,"Adj")
  sf <- attrs(scan.list,"Adj.subfun")
  mode <- attrs(scan.list,"mode")

  obs.P <- determine_obsProb(scan.list = scan.list,sampling = sampling,all.sampled = all.sampled)

  groupSampled <-
    sLapply(scan.list,
            \(s) {
              s[sf(s)] <-
                s[sf(s)] |>
                {\(x) ifelse(stats::rbinom(x,1L,obs.P[sf(obs.P)]) == 1,x,NA)}()
              s |>
                resolve_NA(mode = mode)
            }
    )

  groupSampled <- copy_attrs_to(scan.list,groupSampled)
  attrs(groupSampled,"obs.P") <- obs.P
  groupSampled
}

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param sampling TO WRITE
#' @param all.sampled TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @importFrom scales rescale_max
#' @keywords internal
determine_obsProb <- function(scan.list,sampling = c("constant","matrix","random","function"),all.sampled) {
  Adj <- attrs(scan.list,"Adj")
  sf <- attrs(scan.list,"Adj.subfun")

  if (is.function(sampling)) {
    sampling_fun <- sampling
    sampling <- "function"

    obs.P <- sampling_fun(Adj)
    attr(obs.P,"sampling_fun") <- sampling_fun
  } else {
    switch(typeof(sampling),
           "integer" =,
           "double" = {
             n <- nrow(Adj)
             obs.P <- matrix(sampling,nrow = n,ncol = n,dimnames = dimnames(Adj))
             if (inherits(sampling,"array")) sampling <- "matrix"
             else sampling <- "constant"
           },
           "character" = {
             sampling <- match.arg(sampling)
             switch(sampling,
                    "constant" = stop("Please input a `sampling` variable as a numeric constant/scalar directly (cf. `?group_sample`"),
                    "matrix" = stop("Please input a `sampling` variable as a numeric matrix directly (cf. `?group_sample`"),
                    "function" = stop("Please input a `sampling` variable as a function object directly (cf. `?group_sample`"),
                    "random" = {
                      obs.P <- Adj
                      obs.P[sf(obs.P)] <- runif(length(obs.P[sf(obs.P)]))
                    },
                    # more group-scan sampling "types" can be added here
                    stop("Inputted `sampling` variable not recognized.")
             )
           }
    )
  }
  if (any(obs.P[sf(obs.P)] > 1)) obs.P[sf(obs.P)] <- scales::rescale_max(obs.P[sf(obs.P)])
  obs.P[!sf(obs.P)] <- 0
  attr(obs.P,"sampling") <- sampling
  obs.P
}

## Focal-scan sampling ----

#'  TO WRITE
#'
#' @param scan.list TO WRITE
#' @param sampling TO WRITE
#' @param all.sampled TO WRITE
#'
#' @return  TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
focal_sample <- function(scan.list,sampling = c("even","random","function"),all.sampled = TRUE){
  focalList <- draw_focalList(scan.list = scan.list,
                              sampling = sampling,
                              all.sampled = all.sampled
  )
  focalSampled <- mask_non.focals(scan.list = scan.list,focalList = focalList)
  focalSampled <- copy_attrs_to(scan.list,focalSampled)
  attrs(focalSampled,"focalList") <- focalList
  focalSampled
}

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param sampling TO WRITE
#' @param sampling_fun TO WRITE
#' @param all.sampled TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @keywords internal
draw_focalList<- function(scan.list,sampling = c("even","random","function"),all.sampled = TRUE) {
  if (is.function(sampling)) {
    sampling_fun <- sampling
    sampling <- "function"
  } else {
    sampling_fun <- NULL
    sampling <- match.arg(sampling)
  }

  # retrieve variables from scan.list
  Adj <- attrs(scan.list,"Adj")
  n <- Adj |> nrow()
  node_names <- Adj |> rownames()
  n.scans <- attrs(scan.list,"n.scans")

  # shape future focal.list, filling it with NAs
  focal.list <- rep(NA,n.scans);

  switch(sampling,
         # manage the case of an even focal.list
         "even" = {
           focal.list <-
             c(rep(1:n,n.scans %/% n),sample(1:n,n.scans %% n)) |>
             {\(f) f[sample(seq_along(f))]}()
         },
         # if not even, select at least each each node once, and adjust the rest of the sampling effort needed
         "random" = {
           if (all.sampled){
             if (n > n.scans) stop("n.scans is too small to sample all nodes.")
             focal.list[sample(1:n.scans,n)] <- 1:n; n.scans <- n.scans - n
           }

           focal.list[is.na(focal.list)] <- quick_sample(1:n,n.scans)
           focal.list
         },
         "function" = {
           if (is.null(sampling_fun)) {stop("Please input a `sampling` variable as a function object directly (cf. `?focal_sample`")}
           if (all.sampled) {
             if (n > n.scans) stop("n.scans is too small to sample all nodes.")
             focal.list[sample(1:n.scans,n)] <- 1:n;n.scans <- n.scans - n
           }
           # applies the user-defined function, adjust the minimum probability to be non zero
           P <- sampling_fun(Adj)
           if (any(P == 0)) P <- P + min(P[P > 0])
           # replace remaining NAs for each scan with a node given their probability distribution at each scan
           focal.list[is.na(focal.list)]<- sample(1:n,n.scans,replace = TRUE,prob = P)
           focal.list
         },
         stop("`sampling` class not recognized.")
  )
  names(focal.list) <- node_names[focal.list] # keeps nodes names `NULL` if `Adj` didn't have nodes names
  focal.list
}

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param focalList TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @keywords internal
mask_non.focals <- function(scan.list,focalList) {
  sf <- attrs(scan.list,"Adj.subfun")

  vapply(1:dim(scan.list)[3],
         \(s) {
           scan <- scan.list[,,s]
           foc <- focalList[s]

           scan[-foc,-foc] <- NA
           scan[!sf(scan)] <- 0L
           scan
         },scan.list[,,1]
  )
}

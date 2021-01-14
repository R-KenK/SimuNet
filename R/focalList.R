# focal.list generation and related functions -------------------------------

#' Generator for `focalList` objects
#'
#' @param Adj square integers matrix of occurrences of dyads.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param focal.prob_fun either:
#' \itemize{
#'   \item{Special case `"even"` (default) tries to even out the `focal.list` as much as possible before drawing randomly following a uniform distribution}
#'   \item{`NULL` or `"random"`, pick focals following a uniform distribution}
#'   \item{a user-defined function of (n,Adj) that output a weight of being focal for each node (passed as the `prob` argument to `base::sample` function)}
#' }
#' @param all.sampled logical, should all individuals be sampled before letting them be sampled according to `focal.prob_fun`? Ignored if `focal.prob_fun` is `"even"` (because all nodes will be sampled anyway. Returns an error if `total_scan` is smaller than the number of nodes.
#'
#' @return an `focalList` object (S3 class) containing:
#' \itemize{
#'   \item{`focals`: a vector of focals (as integers)}
#'   \item{`Adj`: inputted `Adj`}
#'   \item{`total_scan`: inputted `total_scan`}
#'   \item{`focal.prob_type`: character scalar, either:
#'     \item{`"even"`: tries to even out the `focal.list` as much as possible before drawing randomly following a uniform distributionall dyad have the same probability of being sampled or not.}
#'     \item{`"random"`: if all node are equiprobable at each scan}
#'     \item{`"user-defined function"`: if the user inputted a function of (n,Adj) to calculate each node probability of being drawn at each scan}
#'   }
#'   \item{`focal.prob_fun`: either:
#'     \item{`"even"`}
#'     \item{`"random"`}
#'     \item{inputted `focal.prob_fun` of (n,Adj)}
#'   }
#'   \item{`all.sampled`: inputted `all.sampled`}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
#' Adj
#'
#' generate_focalList(Adj,total_scan,focal.prob_fun = "even")
#' generate_focalList(Adj,total_scan,focal.prob_fun = "random",all.sampled = FALSE)
#'
#' # using a user-defined function:
#' user_function.n<- function(n,Adj) {1:n} # comparable to a dyad-trait-based bias
#' user_function.n2<- function(n,Adj) {1:n*1:n} # comparable to a dyad-trait-based bias
#' user_function.Adj<- function(n,Adj) {colSums(Adj*Adj)} # comparable to a network-based bias
#'
#' generate_focalList(Adj,total_scan,focal.prob_fun = user_function.n,all.sampled = FALSE)
#' generate_focalList(Adj,total_scan,focal.prob_fun = user_function.n2,all.sampled = TRUE)
#' generate_focalList(Adj,total_scan,focal.prob_fun = user_function.Adj,all.sampled = FALSE)
#'
generate_focalList<- function(Adj,total_scan,
                               focal.prob_fun = "even",all.sampled = TRUE){
  n<- nrow(Adj);nodes_names<- rownames(Adj)

  focal.prob_type<- determine_focal.prob_type(focal.prob_fun = focal.prob_fun)

  focal.list<- list(
    focals = draw_focal.list(focal.prob_fun = focal.prob_fun,all.sampled = all.sampled,focal.prob_type = focal.prob_type,
                             n = n,Adj = Adj,total_scan = total_scan,nodes_names = nodes_names), #perhaps should enrich the type of this element with a vector of names
    Adj = Adj,
    total_scan = total_scan,
    focal.prob_type = focal.prob_type,
    focal.prob_fun = focal.prob_fun,
    all.sampled = all.sampled # Only keep track of the inputted value, so it's possible that all nodes were sampled after the random draw even if this was set to FALSE
  )

  class(focal.list)<- "focalList"
  focal.list
}

#' Print method for `focalList` objects
#' @export
#' @noRd
print.focalList<- function(x,...){
  if(!is.null(names(x$focals))){cat("  node name: ",names(x$focals),"\n")}
  cat(" node index: ",x$focals)
  # print.default(x$focals,...)
}

#' Test if object if a `focalList` object
#'
#' @param x an object to test.
#'
#' @return logical, TRUE if the inputted object is a `focalList` object.
#'
#' @noRd
is.focalList<- function(x){
  inherits(x,"focalList")
}

#' Determine the type of `focalList` objects inputted or to create.
#' from the class of `focal.prob_type` inputted.
#'
#' @param focal.prob_type either:
#' \itemize{
#'   \item{`"even"`: tries to even out the `focal.list` as much as possible before drawing randomly following a uniform distributionall dyad have the same probability of being sampled or not.}
#'   \item{`"random"`: if all node are equiprobable at each scan}
#'   \item{`"user-defined function"`: if the user inputted a function of (n,Adj) to calculate each node probability of being drawn at each scan}
#' }
#'
#' @return a character scalar:
#' \itemize{
#'   \item{`"even"`}
#'   \item{`"random"`}
#'   \item{or `"user-defined function"`}
#' }
#'
#' @noRd
determine_focal.prob_type<- function(focal.prob_fun){
  switch(class(focal.prob_fun),
         "NULL" = "random",
         "character" = switch(focal.prob_fun,
                              "even" = focal.prob_fun,
                              "random" = focal.prob_fun,
                              stop("`focal.prob_fun` character-type not recognized.")
         ),
         "function" = "user-defined function",
         stop("`focal.prob_fun` class not recognized.")
  )
}

#' Draw a vector of focals
#' Internal use. Users should rather use `generate_focalList`. Apply the inputted `focal.prob_fun` to draw a vector of focals.
#'
#' @param focal.prob_fun either:
#' \itemize{
#'   \item{Special case `"even"` tries to even out the `focal.list` as much as possible before drawing randomly following a uniform distribution}
#'   \item{`NULL` or `"random"`, pick focals following a uniform distribution}
#'   \item{a user-defined function of (n,Adj) that output a weight of being focal for each node (passed as the `prob` argument to `base::sample` function)}
#' }
#' @param all.sampled logical, should all individuals be sampled before letting them be sampled according to `focal.prob_fun`? Returns an error if total_scan is smaller than the number of nodes.
#' @param focal.prob_type character scalar:
#' \itemize{
#'   \item{`"even"`}
#'   \item{`"random"`}
#'   \item{or `"user-defined function"`}
#' }
#' @param n integer, number of node in `Adj`.
#' @param Adj square integers matrix of occurrences of dyads.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param nodes_names character vector or `NULL`, names of the nodes in `Adj`
#'
#' @importFrom runif stats
#'
#' @return a named vector of focals (as integers)
#' @noRd
draw_focal.list<- function(focal.prob_fun,all.sampled,focal.prob_type,
                           n,Adj,total_scan,nodes_names) {
  # shape future focal.list, filling it with NAs
  focal.list<- rep(NA,total_scan);

  focal.list<- switch(focal.prob_type,
                      # manage the case of an even focal.list
                      "even" = {
                        focal.list<- c(rep(1:n,total_scan%/%n),sample(1:n,total_scan%%n,replace = FALSE))
                        focal.list[sample(seq_along(focal.list))]
                      },
                      # if not even, select at least each each node once, and adjust the rest of the sampling effort needed
                      "random" = {
                        if(all.sampled){
                          if(n > total_scan){stop("total_scan is too small to sample all nodes.")}
                          focal.list[sample(1:total_scan,n)]<- 1:n;total_scan<- total_scan-n;
                        }

                        focal.list[is.na(focal.list)]<- ceiling(stats::runif(total_scan,0,n))
                        focal.list
                      },
                      "user-defined function" = {
                        if(all.sampled){
                          if(n > total_scan){stop("total_scan is too small to sample all nodes.")}
                          focal.list[sample(1:total_scan,n)]<- 1:n;total_scan<- total_scan-n;
                        }
                        # applies the user-defined function, adjust the minimum probability to be non zero
                        P<- focal.prob_fun(n,Adj);if(any(P==0)){P<-P+min(P[P>0])}
                        # replace remaining NAs for each scan with a node given their probability distribution at each scan
                        focal.list[is.na(focal.list)]<- sample(1:n,total_scan,replace = TRUE,prob = P)
                        focal.list
                      }
                      ,
                      stop("`focal.prob_type` class not recognized.")
  )
  names(focal.list)<- nodes_names[focal.list] # keeps nodes names `NULL` if `Adj` didn't have nodes names
  focal.list
}


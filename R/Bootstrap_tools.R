
# Bootstrap object data management ----------------------------------------

#' Keep bootstrap parameters (and more) as attributes
#' Internal use in Boot_scans(). Add "method", "keep", "mode", "output" attributes to be more easily retrieved by the get function
#'
#' @param Bootstrap Boot_scans() intermediate output
#' @param ... Any named attribute to be included to the Bootstrap object. Mostly determined by the code of Boot_scan(). Will produce error later if no name is attributed to any
#'
#' @return Bootstrap object, with stored attribute for later retrieval through Bootstrap_get.attributes(). A special attribute, attr.list stores those who have been properly inputted, to account for those with special structure.
#' @export
#'
#' @examples
#' #Internal
Bootstrap_add.attributes<- function(Bootstrap,...){
  attr(Bootstrap,"attr.list")<- list(...)
  Bootstrap
}

#' Retrieve bootstrap object's attributes and reassign them in caller frame
#'
#' @param Bootstrap Boot_scans() output
#' @param a character, parameter stored as attribute to retrieve from `Bootstrap`.
#'
#' @return the value of the parameter `a`.
#' @export
#'
#' @examples
#' #Internal
Bootstrap_get.attr<- function(Bootstrap,a){
  attr(Bootstrap,"attr.list")[[a]]
}


#' Retrieve data from specific method from Boot_scans() output
#' Subset rich Bootstrap output choosing what's needed
#'
#' @param Bootstrap Bootstrap output object
#' @param what character scalar, data type requested ("theoretical","group" or "focal")
#' @param get.format character scalar, output format type requested ("list","adjacency", "observed_edges" or "all").
#'
#' @return list which structure depends on chosen data type and Bootstrap attributes
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = FALSE,obs.prob=0.7,
#'                        method = "group",mode = "directed",output = "list")
#' Boot_get.list(Bootstrap,"theoretical")
#' Boot_get.list(Bootstrap,"group")
#' Boot_get.list(Bootstrap,"group","adj")
#' Boot_get.list(Bootstrap,"group","obs")
#'
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = FALSE,obs.prob=0.7,
#'                        method = "group",mode = "directed",output = "adjacency")
#' Boot_get.list(Bootstrap,"theoretical","adj")
#' # Boot_get.list(Bootstrap,"group","list")
#' Boot_get.list(Bootstrap,"group","adj")
#' Boot_get.list(Bootstrap,"group","obs")
#'
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = TRUE,obs.prob=0.7,keep=TRUE,
#'                        method = "both",mode = "directed",output = "all")
#' Boot_get.list(Bootstrap,"focal","all")
#' Boot_get.list(Bootstrap,"group","adjacency")
#' Boot_get.list(Bootstrap,"focal","obs")
#' Boot_get.list(Bootstrap,"group","list")
Boot_get.list<- function(Bootstrap,what=c("theoretical","group","focal"),
                         get.format = c("list","adjacency","observed_edges","all")){
  what<- match.arg(what)
  get.format<- match.arg(get.format)

  method<- Bootstrap_get.attr(Bootstrap,"method");
  output<- Bootstrap_get.attr(Bootstrap,"output");
  total_scan<- Bootstrap_get.attr(Bootstrap,"total_scan");
  n.boot<- Bootstrap_get.attr(Bootstrap,"n.boot");
  scaled<- Bootstrap_get.attr(Bootstrap,"scaled");
  mode<- Bootstrap_get.attr(Bootstrap,"mode");
  use.rare.opti<- Bootstrap_get.attr(Bootstrap,"use.rare.opti");
  focal.list<- Bootstrap_get.attr(Bootstrap,"focal.list");

  if(what!="theoretical" & method!="both" & what!=method){stop("Element requested unavailable in `",substitute(Bootstrap),"`.")}

  switch(output,
         "list" = switch(get.format,
                         "list" = lapply(Bootstrap,function(boot) lapply(boot, function(scan) scan[[what]])),
                         "adjacency" = lapply(Bootstrap,
                                              function(boot){
                                                sum_up.scans(scan_list = boot,focal.list = focal.list,
                                                             scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)[[what]]
                                              }
                         ),
                         "observed_edges" = lapply(Bootstrap,
                                                   function(boot){
                                                     attr(
                                                       sum_up.scans(scan_list = boot,focal.list = focal.list,
                                                                  scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)[[what]],
                                                       "observed_edges"
                                                     )
                                                   }
                         ),
                         "all" = {
                           list(
                             list = lapply(Bootstrap,function(boot) lapply(boot, function(scan) scan[[what]])),
                             adjacency = lapply(Bootstrap,
                                                function(boot){
                                                  sum_up.scans(scan_list = lapply(boot, function(scan) scan[[what]]),
                                                               scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)
                                                }
                             )
                           )
                         }
         ),
         "adjacency" = switch(get.format,
                              "list" = stop("Only summed-up adjacency matrices have been stored in Bootstrap object."),
                              "adjacency" = lapply(Bootstrap,function(boot) boot[[what]]),
                              "observed_edges" = lapply(Bootstrap,function(boot) attr(boot[[what]],"observed_edges")),
                              "all" = {
                                warning("Only summed-up adjacency matrices have been stored in Bootstrap object.")
                                lapply(Bootstrap,function(boot) boot[[what]])
                              }
         ),
         "all" = switch(get.format,
                        "list" = lapply(Bootstrap,function(boot) lapply(boot[["list"]], function(scan) scan[[what]])),
                        "adjacency" = lapply(Bootstrap,function(boot) boot[["adjacency"]][[what]]),
                        "observed_edges" = lapply(Bootstrap,function(boot) attr(boot[["adjacency"]][[what]],"observed_edges")),
                        "all" = {
                          list(
                            list = lapply(Bootstrap,function(boot) lapply(boot[["list"]], function(scan) scan[[what]])),
                            adjancency = lapply(Bootstrap,function(boot) boot[["adjacency"]][[what]])
                          )
                        }
         )
  )
}

#' Scale unscaled bootstrap list
#'
#' @param Bootstrap Bootstrap output object (/!\ WITH output = "adjacency", other output to be implemented if relevant)
#'
#' @return a Bootstrap output object with scaled adjacencies (comparable )
#' @export
#'
#' @examples
#' # Internal use.
Bootstrap.list_post.scale<- function(Bootstrap){
  attr.list<- attr(Bootstrap,"attr.list");
  attr.list[["output"]]<- "adjacency"
  attr.list[["scaled"]]<- TRUE

  Bootstrap.scaled<- lapply(Bootstrap,
                            function(boot){
                              Adj.scaled<- lapply(names(boot),
                                                  function(method){
                                                    observed_edges<- attr(boot[[method]],"observed_edges")
                                                    observed_edges<- ifelse(observed_edges!=0,observed_edges,1) # supposedly boot[[method]]==0 when observed_edges==0, but avoid dividing by 0
                                                    scaled<- boot[[method]]/observed_edges
                                                    diag(scaled)<-0
                                                    scaled
                                                  }
                              )
                              names(Adj.scaled)<- names(boot)
                              Adj.scaled
                            }
  )
  attr(Bootstrap.scaled,"attr.list")<- attr.list
  Bootstrap.scaled
}


#' Retrieve specific simulation parameters of given Bootstrap
#' Internal use. To ease the recollection of a given bootstrap performed through Boot_scans() iterations alongside a parameters.list
#'
#' @param Bootstrap a Bootstrap object (with correct parameter attributes)
#'
#' @return a data frame referencing the simulation parameters
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Boot_get.param<- function(Bootstrap){
  method<- Bootstrap_get.attr(Bootstrap,"method");
  output<- Bootstrap_get.attr(Bootstrap,"output");
  total_scan<- Bootstrap_get.attr(Bootstrap,"total_scan");
  n.boot<- Bootstrap_get.attr(Bootstrap,"n.boot");
  scaled<- Bootstrap_get.attr(Bootstrap,"scaled");
  mode<- Bootstrap_get.attr(Bootstrap,"mode");
  use.rare.opti<- Bootstrap_get.attr(Bootstrap,"use.rare.opti");

  obs.prob<- Bootstrap_get.attr(Bootstrap,"obs.prob");
  focal.list<- Bootstrap_get.attr(Bootstrap,"focal.list");

  data.frame(obs.prob.type = attr(obs.prob,"type"),
             obs.prob.subtype = attr(obs.prob,"subtype"),
             obs.prob.fun = attr(obs.prob,"fun"),
             focal.list.type = attr(focal.list,"type"),
             focal.list.subtype = attr(focal.list,"subtype"),
             total_scan = total_scan,
             mode = mode,
             scaled = scaled,
             use.rare.opti = use.rare.opti
  )
}

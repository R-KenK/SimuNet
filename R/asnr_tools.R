# user-friendly wrappers --------------------------------------------------

#' Import network from asnr package
#' Bansal lab's Animal Social Network Repository
#'
#' @param class character scalar. Phylogenetic class folder of the network to
#'   import (e.g. "Aves","Mammalia",etc.). Supports partial matching.
#' @param species character scalar. species folder of the network to import
#'   (e.g. "Aves","Mammalia",etc.). Needs to match the beginning of the folder's
#'   name Supports partial matching.
#' @param network character scalar. Optional if the folder contains only one
#'   .graphml file. Otherwise the .graphml file name to import. Supports partial
#'   matching.
#' @param url character scalar. Optional if argument `class` and `species` (and
#'   `network` in the case of several graphml in one folder). URL of the
#'   .graphml file to import.
#' @param output character scalar. Either "graph" for an `igraph` graph object,
#'   or "adjacency" for an adjacency matrix
#' @param type character scalar. One of "both", "upper", and "lower". In the
#'   case of undirected network, "upper" or "lower" should probably be
#'   preferred.
#' @param ... additional argument. Useful to pass a `asnr.df` argument to
#'   construct_full.path.
#' @param default_prefix character scalar. URL "prefix" used to retrieve the
#'   graphml _file_ from github.
#' @param full.path character scalar. Optional, can be used to input directly an
#'   URL to a graphml file.
#'
#' @return According to output, either a igraph object or an adjacency matrix.
#' @export
#'
#' @examples
#' import_from_asnr("Aves","cowbird")
#' # partial matching works, but at the moment not case insensitive
#' import_from_asnr("Amph","frog",output = "adj")
#' # Users can also copy-paste a `.graphml`'s URL
#' import_from_asnr(url = "https://github.com/bansallab/asnr/blob/master/Networks/Reptilia/
#' lizard_proximity_weighted/weighted_network_social_T_rugosa.graphml")
#' # To avoid multiple get querries and longer computation time, users can first
#' # load a `asnr.df` data frame by calling once the `asnr_network_df` function
#' asnr.df <- asnr_network_df()
#' import_from_asnr("Amph","frog",output = "adj",asnr.df = asnr.df)
import_from_asnr <- function(class = NULL,
                             species = NULL,
                             network = NULL,
                             url = NULL,
                             output = c("graph","adjacency"),
                             type = c("both", "upper", "lower"),
                             ...,
                             default_prefix = "https://raw.githubusercontent.com/bansallab/asnr/master/Networks/",
                             full.path = NULL) {
  if (is.null(full.path)) {
    if (is.null(url)) {
      full.path <- construct_full.path(class = class,
                                       species = species,
                                       network = network,
                                       default_prefix = default_prefix,
                                       ...
      )
    } else {
      Networks.post <- sub(pattern = ".*Networks/",replacement = "",url)
      full.path <- paste0(default_prefix,Networks.post)
    }
  }
  import_from_graphml(path = full.path,output = output,type = type)
}

#' Generate a data frame of all networks in the Animal Social Network Repository (asnr)
#'
#' @param user character scalar, used to reconstruct the asnr github repository.
#'   Default is "bansallab".
#' @param repo character scalar, used to reconstruct the asnr github repository.
#'   Default is "asnr".
#'
#' @return a character data frame with the column "class", "species", and
#'   "network"
#' @export
#'
#' @examples
#' asnr_network_df()
asnr_network_df <- function(user = "bansallab", repo = "asnr") {
  asnr_graphmls <- retrieve_asnr_graphmls(user = user,repo = repo)

  class.species.network <- substr(asnr_graphmls,nchar("Networks/") + 1,nchar(asnr_graphmls))
  classes <- left_of_slash(class.species.network)
  species.network <- right_of_slash(class.species.network,classes)
  species <- left_of_slash(species.network)
  networks <- right_of_slash(species.network,species)

  data.frame(class = classes,
             species = species,
             network = networks
  )
}


# asnr import internals ---------------------------------------------------

#' Reconstruct an URL to a desired graphml object in the ASNR repository
#'
#' @param class character scalar. Phylogenetic class folder of the network to
#'   import (e.g. "Aves","Mammalia",etc.). Supports partial matching.
#' @param species character scalar. species folder of the network to import
#'   (e.g. "Aves","Mammalia",etc.). Needs to match the beginning of the folder's
#'   name Supports partial matching.
#' @param network character scalar. Optional if the folder contains only one
#'   .graphml file. Otherwise the .graphml file name to import. Supports partial
#'   matching.
#' @param asnr.df output of the `asnr_network_df`function, otherwise a data
#'   frame containing a "class", "species", and "network" columns used to
#'   reconstruct a URL to a graphml file
#' @param default_prefix character scalar. URL "prefix" used to retrieve the
#'   graphml _file_ from github.
#'
#' @return an URL (character scalar) to the desired graphml file
#' @noRd
construct_full.path <- function(class = NULL,
                                species = NULL,
                                network = NULL,
                                asnr.df = asnr_network_df("bansallab","asnr"),
                                default_prefix = "https://raw.githubusercontent.com/bansallab/asnr/master/Networks/") {
  is_unique <- grep(class,unique(asnr.df$class),value = TRUE,fixed = TRUE)
  if (length(is_unique) > 1) {
    stop("Input `class` matches with different class folders.")
  }
  cla <- match.arg(class,unique(asnr.df$class))

  is_unique <- grep(species,unique(asnr.df$species),value = TRUE,fixed = TRUE)
  if (length(is_unique) > 1) {
    warning("Input `species` matches with different species folders. The first one matching has been used.")
  }
  sp <- match.arg(species,unique(asnr.df$species))

  asnr.df.sub <- subset(asnr.df,subset = class == cla & species == sp)

  if (nrow(asnr.df.sub) > 1 & is.null(network)){
    warning(paste0("Several .graphml files in '",cla,"/",sp,"/'.\nFirst one selected. Otherwise specify which .graphml file to choose."))
  }
  network <- match.arg(network,unique(asnr.df.sub$network))

  if(is.na(network)) {
    stop("Network not find in provided class and species repo.")
  }
  paste0(default_prefix,cla,"/",sp,"/",network)
}

#' Retrieve the repository structure of the asnr github repository
#'
#' @param user character scalar, used to reconstruct the asnr github repository.
#'   Default is "bansallab".
#' @param repo character scalar, used to reconstruct the asnr github repository.
#'   Default is "asnr".
#'
#' @return a character vector of all the graphml files in the asnr github
#'   repository
#' @importFrom httr GET
#' @importFrom httr stop_for_status
#' @importFrom httr content
#' @noRd
retrieve_asnr_graphmls <- function(user = "bansallab", repo = "asnr") {
  req <- httr::GET(paste0("https://api.github.com/repos/",user,"/",repo,"/git/trees/master?recursive=1"))
  httr::stop_for_status(req)
  filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)
  filelist <- grep("Networks/", filelist, value = TRUE, fixed = TRUE)
  grep(".graphml", filelist, value = TRUE, fixed = TRUE)
}


# internal wrappers to subset repo's string -------------------------------

#' wrapper to select the first part of a string before a slash
#'
#' @param path a file path containing a slash
#'
#' @return the part at the left of the first slash
#' @noRd
left_of_slash <- function(path) {
  sub("\\/.*", "",path)
}

#' wrapper to select the right part of a string after a slash
#'
#' @param path a file path containing a slash
#' @param left the part at the left of the first slash
#'
#' @return the part at the right of the first slash
#' @noRd
right_of_slash <- function(path,left = left_of_slash(path)) {
  substr(path,nchar(left) + 2,nchar(path))
}

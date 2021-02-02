# user-friendly wrappers --------------------------------------------------

#' Import network from asnr package
#' Bansal lab's Animal Social Network Repository
#'
#' @param class TO WRITE
#' @param species TO WRITE
#' @param network TO WRITE
#' @param url TO WRITE
#' @param output TO WRITE
#' @param type TO WRITE
#' @param ... TO WRITE
#' @param default_prefix TO WRITE
#' @param full.path TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' import_from_asnr("Aves","cowbird")
#' # partial matching works, but at the moment not case insensitive
#' import_from_asnr("Amph","frog",output = "adj")
#' # Users can also copy-paste a `.graphml`'s URL
#' import_from_asnr(url = "https://github.com/bansallab/asnr/blob/master/Networks/Reptilia/
#' lizard_proximity_weighted/weighted_network_social_T_rugosa.graphml")
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

#'  TO WRITE
#'
#' @param user TO WRITE
#' @param repo TO WRITE
#'
#' @return TO WRITE
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

#'  TO WRITE
#'
#' @param class TO WRITE
#' @param species TO WRITE
#' @param network TO WRITE
#' @param asnr.df TO WRITE
#' @param default_prefix TO WRITE
#'
#' @return TO WRITE
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
    stop("Input `species` matches with different species folders.")
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

#' TO WRITE
#'
#' @param user TO WRITE
#' @param repo TO WRITE
#'
#' @return TO WRITE
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

#' TO WRITE
#'
#' @param path TO WRITE
#'
#' @return TO WRITE
#' @noRd
left_of_slash <- function(path) {
  sub("\\/.*", "",path)
}

#' TO WRITE
#'
#' @param path TO WRITE
#' @param left TO WRITE
#'
#' @return TO WRITE
#' @noRd
right_of_slash <- function(path,left = left_of_slash(path)) {
  substr(path,nchar(left) + 2,nchar(path))
}

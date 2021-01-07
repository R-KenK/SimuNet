#' Make the list of Observation subsets
#'
#' Subset the dyadic observations into event/sample-centric observations datasets
#' @param DT data.table of the events/samples of interest, for each and relative to which subsets of data will be produced. Subset between datesdinf and dsup.
#' @param Obs data.table of dyadic observations with their date, and the two nodes in column id and tar (for focal id and target of interaction)
#' @param ID.present logical. If TRUE, provide the subset for a given event/sample only if the individual of intrest from DT is present at least once in a subset of dyadic obseravtions.
#' @param n.core number of (logical) core to use. Default is all minus one.
#'
#' @return a list of subset dyadic observations data table. If no observation were within [dinf,sdup], returns an NA with the reason as an attribute. If data were present but ID.present is TRUE, return an empty data.table with the reason as an attribute.
#' @export
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @examples
#' set.seed(42);ID<- letters[1:10];
#'
#' D<- round(runif(n = 200,min = 1,max = 30));
#' DT<- data.table::data.table(date=as.Date(D,origin = "1991-03-30"),
#'                             id=sample(ID,200,replace = TRUE));
#' DT<- add.prepatent(DT,mean.lag = 0, duration = 15,extension = 0,unit = "day");
#'
#' D.Obs<- round(runif(n = 4000,min = 0,max = 15))
#' Obs<- data.table::data.table(date=as.Date(D.Obs-15,origin = "1991-03-30"),
#'                              id=sample(ID,1000,replace = TRUE),
#'                              tar=sample(ID,1000,replace = TRUE))
#'
#' make.Obs.list(DT,Obs)
#' # multi.core<- system.time(make.Obs.list(DT,Obs,ID.present = FALSE))
#' # multi.core
#' # single.core<- system.time(make.Obs.list(DT,Obs,ID.present = FALSE,n.core = 1))
#' # single.core

make.Obs.list<- function(DT,Obs,ID.present=TRUE,n.core=parallel::detectCores()-1){
  if(is.null(DT$date)) stop("No date column in DT.");if(is.null(Obs$date)) stop("No date column in Obs.")
  i<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...

  data.table::as.data.table(DT);data.table::as.data.table(Obs)
  cl <- snow::makeCluster(n.core);doSNOW::registerDoSNOW(cl); `%dopar%`<- foreach::`%dopar%`;
  on.exit(snow::stopCluster(cl),add = TRUE,after = TRUE)
  Obs.list<- foreach::foreach(i=1:nrow(DT),.packages = c("data.table")) %dopar%{
    di= DT$dinf[i];ds= DT$dsup[i];if(ID.present) ID = as.character(DT$id[i]);
    if(nrow(Obs[date>=di&date<=ds])>0){
      if(nrow(Obs[date>=di&date<=ds][ifelse(ID.present,Obs$id==ID|Obs$tar==ID,TRUE)])!=0){
        return(Obs[date>=di&date<=ds])
      }else{x<-data.table::data.table();attributes(x)$reason<-"individual not in observations subset.";return(x)}
    }else{x<-NA;attributes(x)$reason<-"No observation in subset period.";return(x)}
  }
  Obs.list
}

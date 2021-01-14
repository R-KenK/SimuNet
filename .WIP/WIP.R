# outline algorithm ----


minimize_NA<- function(scan){
  if(!is.empiScan(scan)) {scan}
  switch(scan$mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = {
           resolvable_NA<- is.na(scan) & ((t(scan) == 1) & !is.na(t(scan)))  # scan[i,j] = NA & scan[j,i] = 1 =>  scan[i,j]<- 1
           ifelse(resolvable_NA,1,scan)
         },
         "min" = {
           resolvable_NA<- is.na(scan) & ((t(scan) == 0) & !is.na(t(scan)))  # scan[i,j] = NA & scan[j,i] = 0 =>  scan[i,j]<- 1
           ifelse(resolvable_NA,0,scan)
         },
         "plus" = ,
         "directed" = ,
         "upper" = ,
         "lower" =  ,
         "vector" = raw.scan
  )
}



# manual tests to see if apply_mode does what it should, including on weighted networks ---------------


test.scan.NA<- matrix(c(c(0,1,NA,0,NA),
                        c(NA,0,1,1,NA),
                        c(0,1,0,0,1),
                        c(0,0,1,0,0),
                        c(NA,1,0,0,0)),nrow = 5,byrow = TRUE,dimnames = list(letters[1:5],letters[1:5]))

test<- matrix(c(0,1,0,1,
                0,0,0,0,
                1,0,0,1,
                1,0,1,0),nrow = 4,byrow = TRUE,dimnames = list(letters[1:4],letters[1:4]))

apply_mode(test,"directed")
apply_mode(test,"undirected")
apply_mode(test,"max")
apply_mode(test,"min")
apply_mode(test,"plus")

test.weighted<- matrix(c(0,20,5,3,
                         10,0,2,0,
                         30,2,0,0,
                         3,0,1,0),nrow = 4,byrow = TRUE,dimnames = list(letters[1:4],letters[1:4]))

apply_mode(test.weighted,"directed")
apply_mode(test.weighted,"undirected")
apply_mode(test.weighted,"max")
apply_mode(test.weighted,"min")
apply_mode(test.weighted,"plus")

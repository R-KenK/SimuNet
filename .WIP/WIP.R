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


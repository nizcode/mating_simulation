#---graphing hwe ---

#finding heterozygosity
hetr <- function(data){
  sm <- sum(data)
  no_1<-data[2]
  res <- no_1/sm
  
  return(res)
}






#---Odds ratio ---

# Naisha was calculating was ((mum with gene)/(mum no gene))/((dad with gene)/(dad no gene)) = 1.63
# if genotype is 0(AA) the the score is 0 for that person, ie not a carrier at all
# if genotype is 1(Aa) then the score is 0.5 for that person
# if genotype is 2(aa) then the score is 1 for that person

data <- get(load(paste0('gsp_sim1000_100_',d,'.Rdata')))

carr <- function(gty){
  ca1<-gty[2] * 0.5
  ca2<-gty[3] * 1
  totca <- ca1 + ca2
  noc0<- gty[1] * 1
  noc1<- ca2
  totnoc <- noc0 + noc1
  p <- totca/totnoc
  return(p)
}


for(i in 1:5){
  ac <- apply(data[[i]][,1:3],1,carr)
  ac <- matrix(ac)
  data[[i]]<-cbind(data[[i]],ac)
  
}



g2c<-cbind( matrix(data[[4]][,7], ncol = 1),  matrix(data[[3]][,7], ncol = 1) )
oddr <- apply(g2c,1,function(x) x[1]/x[2])

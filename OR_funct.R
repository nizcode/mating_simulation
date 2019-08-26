#---Odds ratio ---

# Naisha was calculating was ((mum with gene)/(mum no gene))/((dad with gene)/(dad no gene)) = 1.63
# if genotype is 0(AA) the the score is 0 for that person, ie not a carrier at all
# if genotype is 1(Aa) then the score is 0.5 for that person
# if genotype is 2(aa) then the score is 1 for that person

# type in which disassortive you want to look at by writing in argument 1 the degree

odds01 <- function(gty){
  co0 <- gty[1]
  co1 <- gty[2]
  co2 <- gty[3]
  rt01 <- co0/co1
  rt21 <- co2/co1
  return(rt01)
}

odds21 <- function(gty){
  co0 <- gty[1]
  co1 <- gty[2]
  co2 <- gty[3]
  rt01 <- co0/co1
  rt21 <- co2/co1
  return(rt21)
}

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





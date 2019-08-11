#!/usr/bin/env Rscript
source('hardy.R')


#---Odds ratio ---

# Naisha was calculating was ((mum with gene)/(mum no gene))/((dad with gene)/(dad no gene)) = 1.63
# if genotype is 0(AA) the the score is 0 for that person, ie not a carrier at all
# if genotype is 1(Aa) then the score is 0.5 for that person
# if genotype is 2(aa) then the score is 1 for that person

# type in which disassortive you want to look at by writing in argument 1 the degree
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 4 
if(numargs != enumargs) {
  print("sample population number, the number of simuations, amount of disassortive mating you want in decimal, and the number of couples you want in parent generation")
  stop("Stopping right here")
}
N <- args[1]
n <- args[2]
d <- args[3]
psize <- args[4]
data <- get(load(paste0('Res_0112/gsp_sim',N,'_',n,'_',d,'_',psize,'_','.Rdata')))

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
avgMalHet <- mean(data[[3]][,6])
avgFemHet <- mean(data[[4]][,6])

write.table(oddr, file = paste0('oddsRatio_Parents_',N,'_',n,'_','_',d,'_',psize,'_','.txt'))
cat(paste('this is the average of the odds ratio for the parents',mean(oddr),'\n', sep = ' '))
cat(paste('this is the average of fathers heterzygotes',avgMalHet,'and average female heterozygotes',avgFemHet,'\n',sep =' '))
#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
source("samputils.R")
library('HardyWeinberg')
source('hardy.R')
source('dissort_mat.R')
library('doParallel')
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer, cores = 8
registerDoParallel(cl)

cat('write in population size, followed by a space and then followed by the number of simulation you want','dissortative mating, yes or no')

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 3 
if(numargs != enumargs) {
  print("sample population number, the number of simuations")
  stop("Stopping right here")
}
N<-as.numeric(args[1])#pop size
n<-as.numeric(args[2])#number of sim
d <-as.character(args[3])#dissortative function yes or no

options(warn=-1)

#the function the simulates 3 generations
sim<-function(N){
  Nd2 <- N/2 
  #ma is the table of the generation1
  #set.seed(120) for debugging
  ma <- matrix(c(1:N, rep(0,Nd2), rep(1,Nd2)), ncol=2, byrow=FALSE)
  # setting up the matrix, with tha males at the top and females at the bottom, makes the mating process easier
  
  ma <- cbind(ma, rep(c(0,1,1,2), N/4))#allelic frequencies, 0,1,2 in ratio of 1:2:1
  
  attrinames <- c('IID', 'SEX', 'GTY')
  colnames(ma) <- attrinames
  
  
  fshuf <- sample((Nd2+1):N, Nd2) #random shuffle of females
  fshufr <- rep(0,Nd2) 
  fshufr[fshuf-Nd2] <- 1:Nd2 # this is the operation that achieves the reverse of the permutation
  ma <- cbind(ma, c(fshuf, fshufr))# this randomly allocates females to males
  ma <- cbind(ma, c(1:Nd2, fshufr)) 
  
  
  
  attrinames <- c(attrinames, 'MTE', 'FID')
  colnames(ma) <- attrinames
  
  noff <- 1 
  
  offnum <- rep(1:Nd2, noff) 
  #ma2 table for generation2
  ma2 <- matrix(offnum, ncol=1)
  attrinames2 <- 'FID'
  colnames(ma2) <- attrinames2
  
  #put in the ma2 table the number of offspring
  #the population automatically cuts by a half
  g2N <- noff*Nd2
  ma2 <- cbind(ma2, givesx(g2N)) # function givesx, randomly decides if the proband is male or female
  attrinames2 <- c(attrinames2, 'SEX')
  colnames(ma2) <- attrinames2
  
  malegts <-ma[1:Nd2,3]
  femalegts <- ma[ma[1:Nd2,4],3]
  gt2 <- mapply(gengt, malegts, femalegts) #function gengt generates a new genotype assignment given two gt's.
  if(noff>1){
    for(i in 2:noff) {
      gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
    }
  }
  
  ma2 <- cbind(ma2, gt2)
  attrinames2 <- c(attrinames2, 'GTY')
  colnames(ma2) <- attrinames2
  
  #this does the same as ma table, ie males at the top 
  ma2 <- ma2[order(ma2[,2]),]
  fmd <- which(ma2[,2] == 1) 
  fmd <- fmd[1]-1 #fmd = population of males
  num_f<-nrow(ma2)-fmd # num_f = population of females
  
  # the problem now is that there might be unequal female to male
  #this checks which sexe there's more of and which ever there's more of, it shiffles them , ie for random allocation
  if(fmd >= g2N/2) { # majority of males
    pa2_m <-sample(1:fmd, fmd)
    pa2_f <-((fmd+1):g2N)
  } else {
    pa2_f <-sample((fmd+1):g2N, g2N-fmd)
    pa2_m <-(1:fmd)
  }
  
  
  #this part randomily gets rid of males, if there's a male majority
  #or randomly gets rid of females, if females>males
  if(fmd > num_f){
    want_m<-pa2_m[1:num_f]
    evr <- c(want_m, pa2_f)
    ma3<-ma2[evr,]
  } else{
    want_f<-pa2_f[1:fmd]
    evr <- c(want_f, pa2_m)
    ma3<-ma2[evr,]
    
  }
  #---dissort here---
  if(d == 'yes'){
    ma3<- dissort(ma3)
    colnames(ma3) <- attrinames
  } else{
  
  Nma3 <- nrow(ma3)
  dNma3 <- Nma3/2
  #ma4 is generation3, and the next part here follows the same process as generation2
  FID <- ma3[1:Nma3,1]
  ma3 <- cbind(ma3, c(rev(FID),FID))
  attrinames2 <- c(attrinames2, 'MTE')
  colnames(ma3) <- attrinames2
  }
  Nma3 <- nrow(ma3)
  dNma3 <- Nma3/2
  
  offnum <- ma3[1:dNma3,1]
  ma4 <- matrix(offnum, ncol=1)
  attrinames2 <- 'FID'
  colnames(ma4) <- attrinames2
  
  g2N <- noff*dNma3
  ma4 <- cbind(ma4, givesx(g2N))
  attrinames2 <- c(attrinames2, 'SEX')
  colnames(ma4) <- attrinames2
  
  malegts <-ma3[1:dNma3,3]
  femalegts <- rev(ma3[(dNma3+1):Nma3,3])
  gt2 <- mapply(gengt, malegts, femalegts)
  if(noff>2){
    for(i in 2:noff) {
      gt2 <- c(gt2, mapply(gengt, malegts, femalegts)) # 
    }
  }
  
  ma4 <- cbind(ma4, gt2)
  attrinames2 <- c(attrinames2, 'GTY')
  colnames(ma4) <- attrinames2
  
  
  
  gty1<-ma[,3]
  gty2<-ma3[,3]
  gty3<-ma4[,3]
  
  return(list(gty1,gty2,gty3)) #returns the genotypes(0,1 and 2s) of each generation
  
}


#this function sorts out the outpur of sim() into rows, for the input of HWChisq()
sort_rows<-function(gty){
  l0<-length(which(gty==0))
  l1<-length(which(gty==1))
  l2<-length(which(gty==2))
  return(c(l0,l1,l2))
}

#This next part creates a table of MM MN NN pval chisq, the output of this script

Result <- list()
for(gtion in 1:3){
  gt <- foreach(i = 1:n, .combine = rbind)%dopar%{#uses more than one core
    sort_rows(sim(N)[[gtion]])
    
    
  }
  gt <- matrix(gt,ncol=3)
  res <- apply(gt,1,HWChisq)
  resher <- apply(gt,1,hetr)
  
  resy <- foreach(i=1:n, .combine = rbind)%dopar%{
    c(res[[i]]$pval,res[[i]]$chisq,resher[i])
    
  }
  Result[[gtion]]<- cbind(gt,resy)
  res_type <- c("MM","MN","NN","pval","chisq","heterozygosity")
  colnames(Result[[gtion]]) <- res_type
}


stopCluster(cl)


Result<-list("Genration1" = Result[[1]],"Genration2" = Result[[2]] ,"Genration3" = Result[[3]] )

for(i in 1:3) {
  write.table(Result[[i]],file=paste0("Sim_gen",i,'_',N,'_',n,'_',d,'.txt'))
}

save(Result,file = paste0('Sim_result',N,'_',n,'_',d,'.Rdata'))

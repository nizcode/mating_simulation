#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
source("samputils.R")
library('HardyWeinberg')
source('hardy.R')
source('dissort_mat.R')
source('sex_spec_mating.R')
library('doParallel')
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer, cores = 8
registerDoParallel(cl)

cat('write in population size, followed by a space and then followed by the number of simulation you want','fraction (in decimal) of disassortive mating you want')

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 4 
if(numargs != enumargs) {
  print("sample population number, the number of simuations, amount of disassortive mating you want in decimal, and the number of couples you want in parent generation")
  stop("Stopping right here")
}
N<-as.numeric(args[1])#pop size
n<-as.numeric(args[2])#number of sim
psize <- as.numeric(args[4]) #number of couples in parents
diss<-as.numeric(args[3])#dissortative function yes or no

options(warn=-1)

#the function the simulates 3 generations
sim<-function(N){
  Nd2 <- N/2 
  #ma is the table of the generation1
  #set.seed(128) #for debugging
  ma <- matrix(c(1:N, rep(0,Nd2), rep(1,Nd2)), ncol=2, byrow=FALSE)
  # setting up the matrix, with tha males at the top and females at the bottom, makes the mating process easier
  
  ma <- cbind(ma, rep(c(0,1,1,2), N/4))#allelic frequencies, 0,1,2 in ratio of 1:2:1
  
  attrinames <- c('IID', 'SEX', 'GTY')
  colnames(ma) <- attrinames
  if(diss == 1){
    ma <- genDis(ma)
  }else{
  mam <- ma[1:Nd2,]
  maf <- ma[(Nd2+1):N,]
  madm <- mam[1:(Nd2*diss),]
  madf <- maf[1:(Nd2*diss),]
  mad <- rbind(madm,madf)
  
  
  marm <- ma[(Nd2*diss+1):Nd2,] 
  marf <- ma[(Nd2 + (Nd2*diss)+1):N,] 
  
  mar <- rbind(marm,marf)
  mad <- genDis(mad)
  
  Nmar<-nrow(mar)
  Nmar2 <- Nmar/2
  fid <- mar[1:Nmar2,1]
  # so we can cbind() to our matrix, right?
  # Well, no. We need to give the female rows the number of the male they matched with
  # so the reverse of rshuf .. actually this is a tiny bit tricky
  mar <- cbind(mar, c(mar[(Nmar2+1):Nmar,1],fid))
  mar <- cbind(mar, c(mar[1:Nmar2,1],rev(mar[1:Nmar2,1])))
  
  #fshuf <- sample((Nd2+1):N, Nd2) #random shuffle of females
  #fshufr <- rep(0,Nd2) 
  #fshufr[fshuf-Nd2] <- 1:Nd2 # this is the operation that achieves the reverse of the permutation
  #mar <- cbind(mar, c(fshuf, fshufr))# this randomly allocates females to males
  #mar <- cbind(mar, c(1:Nd2, fshufr)) 
  
  ma <- rbind(mad,mar)
  ma <- ma[order(ma[,1],decreasing = FALSE),]
  }
  attrinames <- c(attrinames, 'MTE', 'FID')
  colnames(ma) <- attrinames
  
  N <- nrow(ma)
  Nd2 <- N/2
  
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
  femalegts <- ma[(Nd2+1):N,3]
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
   #p diss
  if(diss == 1){
    ma3 <- genDis(ma3)
    colnames(ma3) <- attrinames
    ma3s<-list(ma3)
  }else{
  sz <- nrow(ma3)
  sz2 <- sz/2
  shum <- sample(1:sz2, sz2)
  shuf <- sample((sz2+1):sz,sz2)
  
  mul <- sz2*diss
  md <- shum[1:mul]
  mf <- shuf[1:mul]
  
  ma2dm <- ma3[md,]
  ma2df <- ma3[mf,]
  
  ma3d <- rbind(ma2dm,ma2df)
  ma3r <- ma2[-c(md,mf),]
  
  
  #--- dissort --- here
  
  
  ma3d<- genDis(ma3d)
  colnames(ma3d) <- attrinames
  
  Nma3 <- nrow(ma3r)
  dNma3 <- Nma3/2
  
  FID <- ma3r[1:Nma3,1]
  # so we can cbind() to our matrix, right?
  # Well, no. We need to give the female rows the number of the male they matched with
  # so the reverse of rshuf .. actually this is a tiny bit tricky
  ma3r <- cbind(ma3r, c(rev(FID),FID))
  attrinames2 <- c(attrinames2,'MTE')
  colnames(ma3r) <- attrinames2
  
  
  ma3s <- list(ma3d,ma3r)
  }
  ma4s <- list()
  
  for(i in 1:length(ma3s)){
    
    Nma3 <- nrow(ma3s[[i]])
    dNma3 <- Nma3/2
    
    offnum <- ma3s[[i]][1:dNma3,1] # FIDs assigned sequentially along
    ma4s[[i]] <- matrix(offnum, ncol=1)
    attrinames2 <- 'FID'
    colnames(ma4s[[i]]) <- attrinames2
    
    # OK, what's next? Well , their sex. That's just random
    # I've turned it into a function, though it hardly needs it
    g2N <- noff*dNma3
    ma4s[[i]] <- cbind(ma4s[[i]], givesx(g2N))
    attrinames2 <- c(attrinames2, 'SEX')
    colnames(ma4s[[i]]) <- attrinames2
    
    malegts <-ma3s[[i]][1:dNma3,3]
    femalegts <- rev(ma3s[[i]][(dNma3+1):Nma3,3])#ma3[ma3[1:dNma3,4],3]
    gt2 <- mapply(gengt, malegts, femalegts)
    # then concatenate the rest
    if(noff>2){
      for(i in 2:noff) {
        gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
      }
    }
    
    ma4s[[i]] <- cbind(ma4s[[i]], gt2)
    attrinames2 <- c(attrinames2, 'GTY')
    colnames(ma4s[[i]]) <- attrinames2
  }  
  
  if(diss != 1){
  ma3 <- rbind(ma3s[[1]][,-5],ma3s[[2]]) 
  ma4 <- rbind(ma4s[[1]],ma4s[[2]]) 
  }else{
    ma3 <- ma3s[[1]]
    ma4 <- ma4s[[1]]
  }
  ma3 <- ma3[order(ma3[,2]),]
  
  
  g2N <- nrow(ma3)
  g2Nd <- g2N/2
  #naishas paper = 1558 males and 1558 females
  nz <- g2Nd - psize
  omi<-sample(ma3[1:g2Nd,1],psize)
  ma3 <- ma3[-(which(ma3[,5]%in%omi)),]
  
  g2N <- nrow(ma3)
  g2Nd <- g2N/2
  
  
  g2m <- ma3[1:g2Nd,3]
  g2f <- ma3[(g2Nd+1):g2N,3]
  
  
  g1m <- ma[1:Nd2,3]
  g1f <- ma[(Nd2+1):N,3]
  
  #gty1<-ma[,3]
  #gty2<-ma3[,3]
  gty3<-ma4[,3]
  
  
  return(list(g1m,g1f,g2m,g2f,gty3,mad)) #returns the genotypes(0,1 and 2s) of each generation
  
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
for(gtion in 1:5){
  gt <- foreach(i = 1:n, .combine = rbind)%dopar%{#uses more than one core

        sort_rows(sim(N)[[1]])
    
    
  }
  
  gt <- matrix(gt,ncol=3)
  res <- apply(gt,1,HWChisq)
  resher <- apply(gt,1,hetr)
  
  resy <- foreach(i=1:n, .combine = rbind)%dopar%{
    c(res[[i]]$pval,res[[i]]$chisq,resher[i])
    
  }
  Result[[gtion]]<- cbind(gt,resy)
  #res_type <- c("MM","MN","NN","pval","chisq","heterozygosity")
  #colnames(Result[[gtion]]) <- res_type
}


stopCluster(cl)


#Result<-list("Gen1_male" = Result[[1]],"Gen1_female" = Result[[2]],"Gen2_male" = Result[[3]],"Gen2_female" = Result[[4]],'Probands' = Result[[5]])

#for(i in 1:5) {
  #write.table(Result[[i]],file=paste0(names(Result[i]),'_',N,'_',n,'_',diss,'_',pszie,'_','.txt'))
#}

save(Result,file = paste0('gsp_sim',N,'_',n,'_',diss,'_',psize,'_','.Rdata'))

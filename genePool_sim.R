#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
source("samputils.R")
library('HardyWeinberg')
source('hardy.R')
source('dissort_mat.R')
source('sex_spec_mating.R')
library('doParallel')
source('dissort3.0.R')
cores=detectCores()
cl <- makeCluster(cores[1]/2) #not to overload your computer, cores = 8
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
  Npd <- N * diss
  Npd2 <- Npd/2
  Npr <- N - Npd
  Npr2 <- Npr/2
  
  #ma is the table of the generation1
  #set.seed(123)# for debugging
  mad <- matrix(c(1:(Npd), rep(0,Npd2), rep(1,Npd2)), ncol=2, byrow=FALSE)
  mar <- matrix(c(1:(Npr), rep(0,Npr2), rep(1,Npr2)), ncol=2, byrow=FALSE)
  
  mad <- cbind(mad, rep(c(0,1,1,2), N/4))#allelic frequencies, 0,1,2 in ratio of 1:2:1
  mar <- cbind(mar, rep(c(0,1,1,2), N/4))#allelic frequencies, 0,1,2 in ratio of 1:2:1
  
  attrinames <- c('IID', 'SEX', 'GTY')
  colnames(mad) <- attrinames
  colnames(mar) <- attrinames
  
  #--- mating starts here --- 
  ma <- list(mad,mar)
  if(diss == 0){
    t  <- 2
  }else if(diss == 1){
    t <- 1
  }else{
    t <- 1:2
  }
  attrinames <- c(attrinames, 'MTE','GtMte', 'FID')
  for(i in t){
  N <- nrow(ma[[i]])
  Nd2 <- N /2
  
  nm <- length(which(ma[[i]][,2]==0))
  nf <- length(which(ma[[i]][,2]==1))
  
  if(nm < nf){
    ma[[i]] <- ma[[i]][-N,]
  }
  
  N <- nrow(ma[[i]])
  Nd2 <- N /2
  
  fshuf <- sample((Nd2+1):N, Nd2) #random shuffle of females
  fshufr <- rep(0,Nd2) 
  fshufr[fshuf-Nd2] <- 1:Nd2 # this is the operation that achieves the reverse of the permutation
  ma[[i]] <- cbind(ma[[i]], c(fshuf, fshufr))
  ma[[i]] <- cbind(ma[[i]], ma[[i]][ma[[i]][(Nd2+1):N,4],3])# this randomly allocates females to males
  ma[[i]] <- cbind(ma[[i]], c(1:Nd2, fshufr)) 
  
  
  
  
  colnames(ma[[i]]) <- attrinames
  
  }
  
  mad <- ma[[1]]
  mar <- ma[[2]]
  
  if(diss > 0){
  N <- nrow(mad)
  Nd2 <- N/2
  TF <- mapply(ranD, mad[1:Nd2,3],mad[1:Nd2,5])
  mad <- cbind(mad , c(TF,rep(NA,Nd2))) 
  #mad[mad[which(mad[1:Nd2,c(4,7)][,2] == 1),4],7] <- 0
  om <- mad[which((mad[,7]) == 0),6]
  #omitF <- mad[c(om),6]
  mad <- mad[-(which(mad[,6]%in%om)),]
  mad<- mad[,-c(7)]
      if( diss == 1){
        ma <- mad
      }else if((diss < 1)&&( diss > 0)){
        #mar <- mar[,-c(5,6)]
        ma <- rbind(mad,mar)
      }
  }else if (diss == 0){
    ma <- mar#[,-c(5,6)]
  }
  

  
  
  
  # --- gen2 starts here ---
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
  #femalegts <- ma[(Nd2+1):N,3]
  femalegts <- ma[1:Nd2,5]
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
  ma3 <- ma3[order(ma3[,2]),]
  #p diss
  # mating of generation 2 starts here 
  
  # this parts split the gen 2 population
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
  ma3r <- ma3[-c(md,mf),]
  
  ma3 <- list(ma3d,ma3r)
  # mating  
  if(diss == 0){
    t  <- 2
    
  }else if(diss == 1){
    t <- 1
  }else{
    t <- 1:2
  }
  
  attrinames <- c(attrinames2, 'MTE','GtMte', 'FID')
  for(i in t){
    N <- nrow(ma3[[i]])
    Nd2 <- N /2
    
    #nm <- length(which(ma[[i]][,2]==0))
    #nf <- length(which(ma[[i]][,2]==1))
    
  #  if(nm < nf){
   #   ma[[i]] <- ma[[i]][-N,]
   # }
    
   
    
    fs <- ma3[[i]][(Nd2+1):N, 1]#random shuffle of females
    fshuf <- sample(fs, Nd2)#random shuffle of females
    yna <- rep(NA,Nd2)
    ma3[[i]] <- cbind(ma3[[i]],c(fshuf,yna),rep(NA,N))
    
    #carries out reverse perm
    for(j in 1:Nd2){
      ma3[[i]][which(ma3[[i]][,1]==ma3[[i]][j,4]),4]<-ma3[[i]][j,1]
    }
    #matches mate with its genotype
    for(j in 1:N){
      ma3[[i]][which(ma3[[i]][,1]==ma3[[i]][j,4]),5] <- ma3[[i]][j,3]
    }
    #family id
    fid<- c(ma3[[i]][1:Nd2,1],ma3[[i]][(Nd2+1):N,4])
    ma3[[i]] <- cbind(ma3[[i]],fid)
    
   
    
    colnames(ma3[[i]]) <- attrinames
    
  }
  
  mad <- ma3[[1]]
  mar <- ma3[[2]]
  
  if(diss > 0){
    N <- nrow(mad)
    Nd2 <- N/2
    TF <- mapply(ranD, mad[1:Nd2,3],mad[1:Nd2,5])
    mad <- cbind(mad , c(TF,rep(NA,Nd2))) 
    #mad[mad[which(mad[1:Nd2,c(4,7)][,2] == 1),4],7] <- 1
    om <- mad[which(mad[,7] == 0),6]
    #omitF <- mad[c(om),6]
    mad <- mad[-(which(mad[,6]%in%om)),]
    mad<- mad[,-c(5,7)]
    if( diss == 1){
      ma3s <- mad
    }else if((diss < 1)&&( diss > 0)){
      mar <- mar[,-c(5)]
      ma3s <- rbind(mad,mar)
    }
  }else if (diss == 0){
    ma3s <- mar[,-c(5)]
  }
  
  ma3s <- ma3s[order(ma3s[,2]),]
  
  
  
  
  
  # generation 3 starts here
  
    
    Nma3 <- nrow(ma3s)
    dNma3 <- Nma3/2
   # fid<- c(ma3s[1:dNma3,1],ma3s[(dNma3+1):Nma3,4])
    #ma3s<-cbind(ma3s,fid)
    #colnames(ma3s)<- c('IID','SEX','GTY','MTE','GTY')
    
    offnum <- ma3s[1:dNma3,1] # FIDs assigned sequentially along
    ma4 <- matrix(offnum, ncol=1)
    attrinames2 <- 'FID'
    colnames(ma4) <- attrinames2
    
    # OK, what's next? Well , their sex. That's just random
    # I've turned it into a function, though it hardly needs it
    g2N <- noff*dNma3
    ma4 <- cbind(ma4, givesx(g2N))
    attrinames2 <- c(attrinames2, 'SEX')
    colnames(ma4) <- attrinames2
    
    malegts <-ma3s[1:dNma3,3]
    femalegts <- rev(ma3s[(dNma3+1):Nma3,3])#ma3[ma3[1:dNma3,4],3]
    gt2 <- mapply(gengt, malegts, femalegts)
    # then concatenate the rest
    if(noff>2){
      for(i in 2:noff) {
        gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
      }
    }
    
    ma4 <- cbind(ma4, gt2)
    attrinames2 <- c(attrinames2, 'GTY')
    colnames(ma4) <- attrinames2
   
  
  
  
  
  g2N <- nrow(ma3s)
  g2Nd <- g2N/2
  #naishas paper = 1558 males and 1558 females
  nz <- g2Nd - psize
  omi<-sample(ma3s[1:g2Nd,1],nz)
  ma3s <- ma3s[-(which(ma3s[,5]%in%omi)),]
  
  g2Nn <- nrow(ma3s)
  g2Nd <- g2Nn/2
  
  
  g2m <- ma3s[1:g2Nd,3]
  g2f <- ma3s[(g2Nd+1):g2Nn,3]
  
  
  g1m <- ma[1:Nd2,3]
  g1f <- ma[(Nd2+1):N,3]
  
  #gty1<-ma[,3]
  #gty2<-ma3[,3]
  gty3<-ma4[,3]
  
  
  return(list(g2m,g2f)) #returns the genotypes(0,1 and 2s) of each generation
  
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
for(gtion in 1:2){
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
  #res_type <- c("MM","MN","NN","pval","chisq","heterozygosity")
  #colnames(Result[[gtion]]) <- res_type
}


stopCluster(cl)


#Result<-list("Gen1_male" = Result[[1]],"Gen1_female" = Result[[2]],"Gen2_male" = Result[[3]],"Gen2_female" = Result[[4]],'Probands' = Result[[5]])

#for(i in 1:5) {
#write.table(Result[[i]],file=paste0(names(Result[i]),'_',N,'_',n,'_',diss,'_',pszie,'_','.txt'))
#}

save(Result,file = paste0('gsp_sim',N,'_',n,'_',diss,'_',psize,'_','.Rdata'))

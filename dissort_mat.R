#---dissortative mating---
#takes in 1 matrix
#0 ->males
#1 ->female


dissort <- function(ma){
N <- nrow(ma)
Nd2 <- N/2

Os<- ma[which(ma[,2]==0),]
Is<- ma[which(ma[,2]==1),]

ma <- rbind(Os,Is) 

mal <- matrix(nrow = N/2, ncol = 3)
fem <- matrix(nrow = N/2, ncol = 3)

#this sets everything up for the loop 
mal <- ma[1:(N/2),]
fem <- ma[(Nd2+1):N,]
nrow(fem)

malshu <- sample(1:nrow(mal), nrow(mal))
mal <- mal[malshu,]
ma <- rbind(mal,fem)

co0<-length(which(fem[,3]==0))
co1<-length(which(fem[,3]==1))
co2<-length(which(fem[,3]==2))
totco<-nrow(fem)




m0<-fem[c(which(fem[,3]==0)),]
m1<-fem[c(which(fem[,3]==1)),]
m2<-fem[c(which(fem[,3]==2)),]

m0<-matrix(m0,ncol = 3)
m1<-matrix(m1,ncol = 3)
m2<-matrix(m2, ncol = 3)

m0mix<-sample(1:co0,co0)
m0mix<-m0[c(m0mix),1]

m1mix<-sample(1:co1,co1)
m1mix<-m1[c(m1mix),1]

m2mix<-sample(1:co2,co2)
m2mix<-m2[c(m2mix),1]

#only use 3 of them in the end, but this is important 
m1m2 <- co1/(co1+co2)
m2m1 <- co2/(co1+co2)

m0m1 <- co0/(co1+co0)
m1m0 <- co1/(co1+co0)

m2m0 <- co2/(co0+co2)
m0m2 <- co0/(co0+co2)

mates<-matrix(nrow=N,ncol=1)
#fmates<-matrix(nrow=Nd2,ncol=1)

# this for loops check if the male geno == 1, if yes then it assigns a 0 or 2 depending on how much of each there's left
# if geno == 0, it can be 1 or 2
# if geno == 2, it can be 1 or 0
# there's a count for geno, and it substracts by when ever it gets taken and assigned to a male
#this changes the ratio and thus changes the odds to make a less favourable to pick than number again

for(i in 1:Nd2){
  #these ratio changes as the count decreases for each
  
  m1m2 <- co1/(co1+co2)
  m1m0 <- co1/(co1+co0)
  m2m0 <- co2/(co0+co2)
    
  
  if((mal[i,3]==0)&&(co1<1)&&(co2<1)) {
    
    mates[i,]<-NA#this is common enough actually whereby there is no negative assortive mate, in this case an NA is assigned and this person will be later omitted
  } else if((mal[i,3]==0)&&(co1>1)&&(co2<1)) {
    mates[i,]<-m1mix[co1]
    co1<-co1-1
  } else if((mal[i,3]==0)&&(co1<1)&&(co2>1)) {
    mates[i,]<-m2mix[co2]
    co2<-co2-1
  
  } else if((mal[i,3]==0)) {
      a <- runif(1)
      if(a<m1m2){
      mates[i,]<-m1mix[co1]
      co1<-co1-1
      }else {
        mates[i,]<-m2mix[co2]
        co2<-co2-1
        }
  
  

  } else if((mal[i,3]==1)&&(co2<1)&&(co0<1)) {
    mates[i,]<-NA
  } else if((mal[i,3]==1)&&(co2>1)&&(co0<1)){
    mates[i,]<-m2mix[co2]
    co2<-co2-1
  } else if((mal[i,3]==1)&&(co2<1)&&(co0>1)) {
    mates[i,]<-m0mix[co0]
    co0<-co0-1
    
  } else if((mal[i,3]==1)) {
    a <- runif(1)
    if(a < m2m0){
      mates[i,]<-m2mix[co2]
      co2<-co2-1
    }else {
      mates[i,]<-m0mix[co0]
      co0<-co0-1
    }
    
  
  }else if((mal[i,3]==2)&&(co1<1)&&(co0<1)) {
    mates[i,]<-NA
  } else if((mal[i,3]==2)&&(co1>1)&&(co0<1)){
    mates[i,]<-m1mix[co1]
    co1<-co1-1
  } else if((mal[i,3]==2)&&(co1<1)&&(co0>1)) {
    mates[i,]<-m0mix[co0]
    co0<-co0-1
    
  } else if((mal[i,3]==2)) {
    a <- runif(1)
    if(a<m1m0){
      mates[i,]<-m1mix[co1]
      co1<-co1-1
    }else {
      mates[i,]<-m0mix[co0]
      co0<-co0-1
    }
    
  }
}
ma <- cbind(ma,mates)


#these would be all the females with no mate
over_0 <- m0mix[0:co0]
over_1 <- m1mix[0:co1]
over_2 <- m2mix[0:co2]
over_mal <- ma[which(is.na(ma[1:Nd2,4])==TRUE),]#all the males with no mate
over_mal <- matrix(over_mal, ncol = 4)
omit <- c(over_0,over_1,over_2,over_mal[,1])


#for now we need to prtend they have a mate in order for the reverse permutation to work
boys <- mal[,1]
girls <- fem[,1]

spl1 <- omit[which(omit %in% boys)]
spl2 <- omit[which(omit %in% girls)]

for(i in 1:length(spl1)){
  ma[which(ma[,1]==spl1[i]),4] <- spl2[i]
}



#reverse permutation, if rows != FID, reason I have an IF statement here is because, if I can avoid doing a loop I will!
if(FALSE %in% (ma[,1]==1:nrow(ma))){
  

for(i in 1:length(boys)){
  ma[which(ma[,1]==girls[i]),4] <- boys[i]

}

fid<- c(ma[1:Nd2,1],ma[(Nd2+1):N,4])
  
ma <- cbind(ma,fid)
#the reverse permutation, if rows = FIDs, reason , if FID == rows, no need for loop = more effecient
}else{
ma <- ma[order(ma[,1],decreasing = FALSE),]
mty <- ma[,4]
maw <- ma[,-4]


fshufr <- rep(0,Nd2)
fshuf <- ma[(Nd2+1):N,]
fshufr[fshuf-Nd2] <- 1:Nd2

ma <- cbind(ma, c(fshuf, fshufr),c(1:Nd2, fshufr))
}





#and here I delete all the mates that are not in assortative mating
for(i in seq(omit)){
  ma<-ma[-(which(ma[,1] == omit[i])),]
}



return(ma)


}


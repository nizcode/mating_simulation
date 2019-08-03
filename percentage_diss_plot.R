#!/usr/bin/env Rscript

#--- plotting ---
library(plotly)
source('hardy.R')
Sys.setenv("plotly_username" = "neilfa")
Sys.setenv("plotly_api_key" = 'tJ5eehh9BxvCwzIuZHhZ')

het2 <- c()
het3 <- c()
pvs2 <- c()
pvs3 <- c()
for(i in 1:10){
  d <- i/10
  
  data <- get(load(paste0('Sim_result10000_1000_',d,'.Rdata')))
  pvs2[i] <- length(which(data$Genration2[,4]<0.05))
  pvs3[i] <- length(which(data$Genration3[,4]<0.05))
  het2[i] <- mean(data$Genration2[,5])
  het3[i] <- mean(data$Genration3[,5])
 
}

pvs <- list(pvs2,pvs3)
het <- list (het2,het3)
ppl <- list()
pht <- list()
for(i in 1:2){
  d <- paste0('Generation',(1+i))
  x <- seq(10,100,by=10)
  yp<- pvs[[i]]
  yh <- het[[i]]
  data3 <- data.frame(x,yp,yh)
  ppl[[i]] <- plot_ly(data3, x = ~x, y= ~yp,name =paste0('pv_number'),type='scatter',mode='lines+markers')%>%
  layout(title = paste0('number of pv < 0.05 in relation to disassortive mating ',d),
         yaxis = list(title= 'pvalue numbers'),
         xaxis = list(title = "% disassortiveness"))
  pht[[i]] <- plot_ly(data3, x = ~x, y= ~yh,name = 'mean heterozygosity',type='scatter',mode='lines+markers')%>%
  layout(title = paste0("mean heterozygosity in relatin to disassortive mating ", d),
         yaxis = list(title='% of mean heterozygotes'),
         xaxis = list (title = "% disassortiveness"))

  
}






plotly_IMAGE(ppl[[2]], format = "png", out_file = paste0("pval_diss",d,".png"))
plotly_IMAGE(pht[[2]], format = "png", out_file = paste0("hetero_diss",d,".png"))



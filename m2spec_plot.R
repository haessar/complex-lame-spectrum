
rm(list=ls())

library(extrafont)
# font_install('fontcm')
# loadfonts()

## Windows
# setwd(dir = "C:/Users/William/Dropbox/Papers_for_publication/") # Home
# setwd(dir = "C:/Users/Will/Dropbox/Papers_for_publication/") # Work
# Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

## Linux
# setwd(dir = "/home/haezar/Dropbox/Papers_for_publication/") # Laptop
setwd('/media/haezar/SAMSUNG/Academia/Research/Papers_for_publication/')
Sys.setenv(R_GSCMD = "/home/haezar/Downloads/ghostscript-9.19-linux-x86_64/gs-919-linux_x86_64")


i = sqrt(as.complex(-1))

load("eigVals_g3p05_2.RData")
load("eigVals_g30.RData")
load("eigVals_g3m05.RData")

sets = list.files(pattern = "eig")

oldparmar = par()$mar

cnt = 0
for (i in sets)
  if (i=="eigVals_g30.RData"){
    plt = 2 # Plot order
    eigVals = rev(eigVals_g30[[1]])
    eigVals = eigVals[1:1000]
    g2=0
    g3=1
  } else if (i=="eigVals_g3m05.RData"){
    plt = 3
    eigVals = rev(eigVals_g3m05[[1]])
    eigVals = eigVals[1:1000]
    g2=-0.5
    g3=1
  } else if (i=="eigVals_g3p05_2.RData"){
    plt = 1
    eigVals = rev(eigVals_g3p05_2[[1]])
    # eigVals = eigVals[1:1000]
    g2=0.5
    g3=1
  }
  
  if (cnt==0){
    pdf(file = "m2spec_all_2.pdf",family = "CM Roman",width = 12,height = 12)
    split.screen(matrix(c(0,0.25,  1, 0.75,  0.5, 0,  1, .5), ncol=4))
    split.screen(figs = c( 1, 2 ), screen=1)
  }
  
  if (plt==1){
    screen(3)
    par(pty="s", mar = c(5.1, 4.6, 4.1, 2.6))
    plot(eigVals,xlim=c(-3,6), ylim=c(-6,6),type = 'n',xlab= "Re",ylab="Im",axes=T,
         main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
    abline(v=0,h=0,lty=2,col="gray")
    cv1 = eigVals[Re(eigVals)<0 & abs(Im(eigVals))>0.1 & Im(eigVals) < 3.5] 
    cv1_tab = data.frame("Re" = Re(cv1), "Im" = Im(cv1))
    cv1_tab = cv1_tab[order(cv1_tab$Im,decreasing = T),]
    cv1 = NULL
    for (j in 1:nrow(cv1_tab)){
      cv1 = c(cv1, cv1_tab$Re[j] + 1i*cv1_tab$Im[j])
    }
    lines(cv1,col="black",lwd=2)
    
    
    cv2 = eigVals[Re(eigVals)<3 & abs(Im(eigVals))<0.001] 
    lines(cv2,col="black",lwd=2)
    
    cv3 = eigVals[Re(eigVals)>=3 & abs(Im(eigVals))<0.001]
    lines(cv3,col="black",lwd=2)
  }
  
  if (plt==2){
    screen(4)
    par(pty="s", mar = c(5.1, 4.6, 4.1, 2.6))
    plot(eigVals,xlim=c(-3,6), ylim=c(-6,6),type = 'n',xlab= "Re",ylab="Im",axes=T,
         main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
    abline(v=0,h=0,lty=2,col="gray")
    eigVals = c(eigVals, 0) # Need to re-add zero to avoid gap at origin
    
    cv1 = eigVals[Re(eigVals) > 3]
    lines(cv1,col="black",lwd=2)
    
    cv2 = eigVals[Re(eigVals)<3 & Im(eigVals)<=0] 
    cv2_tab = data.frame("Re" = Re(cv2), "Im" = Im(cv2))
    cv2_tab = cv2_tab[order(cv2_tab$Im,decreasing = T),]
    cv2 = NULL
    for (j in 1:nrow(cv2_tab)){
      cv2 = c(cv2, cv2_tab$Re[j] + 1i*cv2_tab$Im[j])
    }
    lines(cv2,col="black",lwd=2)
    
    cv3 = eigVals[Re(eigVals)<3 & Im(eigVals)>=0]
    cv3_tab = data.frame("Re" = Re(cv3), "Im" = Im(cv3))
    cv3_tab = cv3_tab[order(cv3_tab$Im,decreasing = T),]
    cv3 = NULL
    for (j in 1:nrow(cv3_tab)){
      cv3 = c(cv3, cv3_tab$Re[j] + 1i*cv3_tab$Im[j])
    }
    lines(cv3,col="black",lwd=2)
  }
  
  if (plt==3){
    screen(2)
    par(pty="s", mar = c(5.1, 4.6, 4.1, 2.6))
    plot(eigVals,xlim=c(-3,6), ylim=c(-6,6),type = 'n',xlab= "Re",ylab="Im",axes=T,
         main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
    abline(v=0,h=0,lty=2,col="gray")
    cv1 = eigVals[Re(eigVals) > 3]
    lines(cv1,col="black",lwd=2)
    
    cv2 = eigVals[Re(eigVals)<3 & Im(eigVals)<=0] 
    lines(cv2,col="black",lwd=2)
    
    cv3 = eigVals[Re(eigVals)<3 & Im(eigVals)>=0]
    lines(cv3,col="black",lwd=2)
    
  
  
  cnt = cnt + 1
  }
  
  if (cnt==length(sets)){
    close.screen(all=TRUE)
    dev.off()
    embed_fonts("m2spec_all_2.pdf", outfile="m2spec_all_2_embed.pdf")
  }
}

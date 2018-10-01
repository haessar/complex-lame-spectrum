# dev.off()


rm(list=ls())
library(elliptic)
library(extrafont)
# setwd(dir = "C:/Users/William/Dropbox/Papers_for_publication/") # Home
setwd(dir = "C:/Users/Will/Dropbox/Papers_for_publication/") # Work
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

i = sqrt(as.complex(-1))


#### cond_solver takes exponentially longer to run for smaller dtx, therefore save results as res2.RData
# load("C:/Users/Will/Dropbox/Papers_for_publication/res4.RData") # Work
# load("C:/Users/William/Dropbox/Papers_for_publication/res4.RData") # Home


cond_solver = function(g2,g3,dtx = 0.1,plot=T,pts=F,lns=F,status=F,equi=F,title=T,equi2=NULL){
  
  rot = 1 #(1-i)/2
  
  if (equi==T){
    om1 <<- Re(half.periods(g=c(g2,g3))[1])
    om3 <<- exp(pi*i/3)*om1
    om = om1
  } else if (equi==F){
    om3 <<- half.periods(g=c(g2,g3))[2]
    om1 <<- half.periods(g=c(g2,g3))[1] + om3  #### half.periods gives om3, om1-om3, respectively
    om = om1
  }
  if (is.null(equi2)==F){
    om1 <<- 2*Re(half.periods(g=c(g2,g3))[1])
    om3 <<- half.periods(g=c(g2,g3))[2]/2 + om1/2
    om = om1
  }
  
  coord1 = (2*om3)*rot
  coord2 = (2*om1 - 2*om3)*rot
  
  fn2 <- function(a,b){ # Function that we want to be equal to zero
    f <- zeta(z=om,g = c(g2,g3))*(a + b)*rot - zeta(z=(a + b)*rot,g = c(g2,g3))*om
    as.numeric(Re(f))
  }
  
  if (plot==T){
    if (title==T){
#       xlim=c(0,2*Re(om1))
#       ylim=c(-2*Im(om3),2*Im(om3))
      xlim=c(0,2*2.662381)
      ylim=c(-2.662381,2.662381)
      
      plot(1,1,type='n',axes=T,bty="L",ylim=ylim,xlim=xlim,xlab="Re",ylab="Im",yaxs="r",xaxs="r",
           cex.lab=1.7,main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.main=2,cex.axis=1.5)
      
#       xyplot(1~1,type="n",ylim=c(-2*Im(om3),2*Im(om3)),xlim=c(0,2*Re(om1)),aspect="iso")
#       ylim=c(-2*Im(om3),2*Im(om3))
#       xlim=c(0,2*Re(om1))
#       rat = (abs(ylim[1]) + abs(ylim[2]) )/ (abs(xlim[1]) + abs(xlim[2]) )
#       rat = (abs(xlim[1]) + abs(xlim[2]) )/ (abs(ylim[1]) + abs(ylim[2]) )
#       
      
#       eqscplot(1,1,ratio = rat,type="n",ylim=c(-2*Im(om3),2*Im(om3)),xlim=c(0,2*Re(om1)),uin=1)
      
    } else {
      xlim=c(0,2*2.662381)
      ylim=c(-2.662381,2.662381)
      plot(1,1,type='n',axes=T,bty="L",ylim=ylim,xlim=xlim,xlab="Re",ylab="Im",yaxs="r",xaxs="r",
           cex.lab=1.7,main = "",cex.main=2,cex.axis=1.5)
    }
    abline(h=0,lty=2,col="gray")
    
    lines(col="gray", (om1 - om3)*seq(-4,4,0.1),lty=2)
    lines(col="gray",2*om1 + (om1 - om3)*seq(-4,4,0.1),lty=2)
    lines(col="gray",2*om3*seq(-4,4,0.1),lty=2)
    lines(col="gray",2*om1 + 2*om3*seq(-4,4,0.1),lty=2)
    
    points(2*om3,pch=4,cex=1.5)
    points(x=Re(2*om1),y=Im(2*om1),pch=4,cex=1.5)
    points(2*om1 - 2*om3,pch=4,cex=1.5)
    
    text(2*om3,labels = bquote(2~omega[3]),pos = 4,cex=1.2,offset = 1)
    text(x=Re(2*om1),y=Im(2*om1),labels = bquote(2~omega[1]),pos = 1,cex=1.2,offset = 1)
    text(2*om1 - 2*om3,labels = bquote(2~omega[1]~"-"~2~omega[3]),pos = 4,cex=1.2,offset = 1)
  }
  
  a = b = seq(0,1,dtx)
  lin1 = a*coord1
  lin2 = b*coord2
  
  ypts=NULL
  xpts=NULL
  
  for (a in 1:length(lin1)){ # Move along x-axis incrementally
    
    if (status==T){
      print(paste0(round(a/length(lin1)*100,1),"%"))
    }
    
    ofv = 0
    for (b in 1:length(lin2)){ # For some x-value, determine fn2 value for all y-axis increments
      ofv[b] = abs(fn2(lin1[a],lin2[b]))
    }
    s_ofv = sort(ofv) # Sort resulting values in ascending order
    ypt = s_ofv[1:2] # For chosen y-values, take the first two
    ypt = lin2[which(ofv %in% ypt)]
    xpt = rep(lin1[a],length(ypt)) # Repeat x-value length(ypt) times to match length of ypt vector.
    if (pts==T) points(xpt + ypt,col="blue") # Plot the resulting points on complex a-b plane.
    if (a==1){
      ypts = ypt
      xpts = xpt
    } else {
      ypts = c(ypts,ypt)
      xpts = c(xpts,xpt)
    }
  }
  
  if (lns==T) lines(xpts,ypts,col="black",lwd=2)
  
  out = data.frame("x" = xpts,"y" = ypts)
  return(out)
}

spec_plotter = function(band1,band2,g2,g3,ylim=c(-5,5),xlim=c(-0.5,5),equi=T,title=T,equi2=NULL){
  
  if (equi==T){
    om1 <<- Re(half.periods(g=c(g2,g3))[1])
    om3 <<- exp(pi*i/3)*om1
    om <<- om1
  } else if (equi==F){
    om3 <<- half.periods(g=c(g2,g3))[2] #### half.periods gives om3, om1-om3, respectively
    om1 <<- half.periods(g=c(g2,g3))[1] + om3 
    om <<- om1
  }
  if (is.null(equi2)==F){
    om1 <<- 2*Re(half.periods(g=c(g2,g3))[1])
    om3 <<- half.periods(g=c(g2,g3))[2]/2 + om1/2
    om = om1
  }
  
  k1 = band1
  k2 = band2
  
  #   plot(-om^2*P(k1,g=c(g2,g3)),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
  #        main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
  #   abline(v=0,h=0,lty=2,col="gray")
  #   
  #   lines(-om^2*P(k1,g=c(g2,g3)),col="black",lwd=2)
  #   lines(-om^2*P(k2,g=c(g2,g3)),col="black",lwd=2)
  
  
  if (title==T){
    plot(-om^2*P(k1,Omega=as.primitive(c(om1,om3))),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
         main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
  } else {
    plot(-om^2*P(k1,Omega=as.primitive(c(om1,om3))),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
         main = "",cex.lab=1.7,cex.main=2,cex.axis=1.5)
  }
  abline(v=0,h=0,lty=2,col="gray")
  
  lines(-om^2*P(k1,Omega=as.primitive(c(om1,om3))),col="black",lwd=2)
  lines(-om^2*P(k2,Omega=as.primitive(c(om1,om3))),col="black",lwd=2)
  
}






## Rhombic k-plane comparison
# pdf(file = "rhombic_comp_sq.pdf",family = "CM Roman",width = 12,height = 12)
pdf(file = "rhombic_comp_sq2.pdf",family = "CM Roman",width = 12,height = 12)
par(mfrow=c(2,2),pty="s")
# par(mfrow=c(2,2))

# 1
# c1 = cond_solver(g2= -1,g3= -1,dtx=0.001,pts=F,equi=F,status=T)
# c5 = cond_solver(g2= 0,g3= -1,dtx=0.001,pts=T,equi=F,status=T,equi2=T)
cond_solver(g2= 0,g3= -1,dtx=0.5,pts=F,equi=F,status=T,equi2=T)

lin = c5[which(abs(Re(c5$x + c5$y) - Re(om1))<0.01),]
curve = c5[which(abs(Re(c5$x + c5$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]
curve2 = curve[which(Re(curve$x + curve$y) <= Re(om1) & Im(curve$x + curve$y) < 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)
curve2 = curve2$x + curve2$y
curve2a = Conj(rev(curve2))
curve2 = c(curve2,curve2a)

lines(lin,lwd=2)
lines(curve1,col="black",lwd=2)
lines(curve2,col="black",lwd=2)

# 2
# options(digits=10)
# c2 = cond_solver(g2= -1,g3= -0.4748432462,dtx=0.001,pts=T,equi=F,status=T)
cond_solver(g2= -1,g3= -0.4748432462,dtx=0.5,pts=F,equi=F,status=T,title=F)
title(main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ "-0.47484..."),cex.main=2)

lin = c2[which(abs(Re(c2$x + c2$y) - Re(om1))<0.01),]
curve = c2[which(abs(Re(c2$x + c2$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]
curve2 = curve[which(Re(curve$x + curve$y) <= Re(om1) & Im(curve$x + curve$y) < 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)
curve2 = curve2$x + curve2$y
curve2a = Conj(rev(curve2))
curve2 = c(curve2,curve2a)

lines(lin,lwd=2)
lines(curve1,col="black",lwd=2)
# lines(curve1a,col="black",lwd=2)
lines(curve2,col="black",lwd=2)
# lines(curve2a,col="black",lwd=2)

# 3
# c3 = cond_solver(g2= -1,g3=0,dtx=0.001,pts=F,equi=F,status=T)
cond_solver(g2= -1,g3= 0,dtx=0.5,pts=F,equi=F,status=T)

lin = c3[which(abs(Re(c3$x + c3$y) - Re(om1))<0.0001),]
curve = c3[which(abs(Re(c3$x + c3$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]
curve2 = curve[which(Im(curve$x + curve$y) <=0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve2 = curve2$x + curve2$y

lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)

# 4
# c4 = cond_solver(g2=-1,g3= 1,dtx=0.005,pts=F,equi=F,status=T)
# c6 = test2
cond_solver(g2= 0,g3= 1,dtx=0.5,pts=F,equi=T,status=T)

lin = c6[which(abs(Re(c6$x + c6$y) - Re(om1))<0.0001),]
curve = c6[which(abs(Re(c6$x + c6$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]
curve2 = curve[which(Im(curve$x + curve$y) <=0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve2 = curve2$x + curve2$y

lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)

dev.off()
# embed_fonts("rhombic_comp_sq.pdf", outfile="rhombic_comp_sq_embed.pdf")
embed_fonts("rhombic_comp_sq2.pdf", outfile="rhombic_comp_sq2_embed.pdf")
##

####

## Rhombic spectrum comparison
pdf(file = "rhombic_spec2.pdf",family = "CM Roman",width = 12,height = 12)
par(mfrow=c(2,2),pty="s")

# 1
cond_solver(g2= 0,g3= -1,dtx=0.5,pts=F,equi=F,status=T,plot=F,equi2 = T)

lin = c5[which(abs(Re(c5$x + c5$y) - Re(om1))<0.01),]
curve = c5[which(abs(Re(c5$x + c5$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)

spec_plotter(band1 = lin,band2 = curve1,g2 = 0,g3 = -1,equi=F,xlim=c(-3,7),ylim=c(-5,5),equi2=T)

# 2
cond_solver(g2= -1,g3= -0.4748432462,dtx=0.5,pts=F,equi=F,status=T,plot=F,title=F)
lin = c2[which(abs(Re(c2$x + c2$y) - Re(om1))<0.01),]
curve = c2[which(abs(Re(c2$x + c2$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)

spec_plotter(band1 = lin,band2 = curve1,g2 = -1,g3 = -0.4748432462,equi=F,xlim=c(-3,7),ylim=c(-5,5),title=F)
title(main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ "-0.47484..."),cex.main=2)


# 3
cond_solver(g2= -1,g3= 0,dtx=0.5,pts=F,equi=F,status=T,plot=F)
lin = c3[which(abs(Re(c3$x + c3$y) - Re(om1))<0.0001),]
curve = c3[which(abs(Re(c3$x + c3$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y

spec_plotter(band1 = lin,band2 = curve1,g2 = -1,g3 = 0,equi=F,xlim=c(-3,7),ylim=c(-5,5))


# 4
cond_solver(g2= 0,g3= 1,dtx=0.5,pts=F,equi=T,status=T,plot=F)
lin = c6[which(abs(Re(c6$x + c6$y) - Re(om1))<0.0001),]
curve = c6[which(abs(Re(c6$x + c6$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y

spec_plotter(band1 = lin,band2 = curve1,g2 = 0,g3 = 1,equi=T,xlim=c(-3,7),ylim=c(-5,5))

dev.off()
embed_fonts("rhombic_spec2.pdf", outfile="rhombic_spec2_embed.pdf")
##

####
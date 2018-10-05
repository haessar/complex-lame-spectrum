rm(list=ls())
library(elliptic)
library(extrafont)
# setwd(dir = "C:/Users/William/Dropbox/Papers_for_publication/") # Home
setwd(dir = "C:/Users/Will/Dropbox/Papers_for_publication/") # Work
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

i = sqrt(as.complex(-1))

### Conditions for Rhombic lattice:
# om1 real and > 0
# Im(om3) > 0
# Re(om3) = 1/2 om1
# Del < 0 

om1 = 1
om3 = i/2 + 1/2*om1 # pseudo-lemniscatic case
om3 = exp(pi*i/3)*om1

g2 = Re(g2.fun(c(om1,om3)))
g3 = Re(g3.fun(c(om1,om3)))


del = g2^3 - 27* g3^2


########################################################################################################
#### FUNCTIONS

cond_solver = function(g2,g3,dtx = 0.1,plot=T,pts=F,lns=F,status=F,equi=F,title=T){

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

  coord1 = (2*om3)*rot
  coord2 = (2*om1 - 2*om3)*rot

  fn2 <- function(a,b){ # Function that we want to be equal to zero
    f <- zeta(z=om,g = c(g2,g3))*(a + b)*rot - zeta(z=(a + b)*rot,g = c(g2,g3))*om
    as.numeric(Re(f))
  }
  
  if (plot==T){
    if (title==T){
      plot(1,1,type='n',axes=T,bty="L",ylim=c(-2*Im(om3),2*Im(om3)),xlim=c(0,2*Re(om1)),xlab="Re",ylab="Im",yaxs="r",xaxs="r",
            cex.lab=1.7,main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.main=2,cex.axis=1.5)
    } else {
      plot(1,1,type='n',axes=T,bty="L",ylim=c(-2*Im(om3),2*Im(om3)),xlim=c(0,2*Re(om1)),xlab="Re",ylab="Im",yaxs="r",xaxs="r",
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

                                                ################

spec_plotter = function(band1,band2,g2,g3,ylim=c(-5,5),xlim=c(-0.5,5),equi=T,title=T){
  
  if (equi==T){
    om1 <<- Re(half.periods(g=c(g2,g3))[1])
    om3 <<- exp(pi*i/3)*om1
    om <<- om1
  } else if (equi==F){
    om3 <<- half.periods(g=c(g2,g3))[2] #### half.periods gives om3, om1-om3, respectively
    om1 <<- half.periods(g=c(g2,g3))[1] + om3 
    om <<- om1
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

####
########################################################################################################



########################################################################################################
#### PSEUDO-LEMNISCATIC CASE

## Pseudo-lemniscatic k-plane alone
pdf(file = "pseudo_lem_plane.pdf",family = "CM Roman",width = 6,height = 6)
par(mfrow=c(1,1),pty="s")

test = cond_solver(g2=-4,g3=0,dtx = 0.001,pts=F)

lin = test[which(abs(Re(test$x + test$y) - Re(om1))<0.0001),]
curve = test[which(abs(Re(test$x + test$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]
curve2 = curve[which(Im(curve$x + curve$y) <=0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve2 = curve2$x + curve2$y

lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)

dev.off()
embed_fonts("pseudo_lem_plane.pdf", outfile="pseudo_lem_plane_embed.pdf")
##

## Pseudo-lemniscatic k-plane and spectrum side-by-side
pdf(file = "pseudo_lem_comp.pdf",family = "CM Roman",width = 12,height = 6)
par(mfrow=c(1,2),pty="s")

cond_solver(g2=-4,g3=0,dtx = 0.05,pts=F)
lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)
spec_plotter(band1 = lin,band2 = curve1,g2 = -4,g3 = 0)

dev.off()
embed_fonts("pseudo_lem_comp.pdf", outfile="pseudo_lem_comp_embed.pdf")
##

####
########################################################################################################


########################################################################################################
#### EQUIANHARMONIC CASE
# 
# test1 = cond_solver(g2 = 0,g3 = 12.82538,dtx = 0.05,pts = T,status = T,equi = T)
# 
# lin1 = i*Im(om3) + seq(Re(om1)/2,3*Re(2*om1)/4,0.001) ### WITH k = lin1, lin2
# lin2 = Re(om1) + i*seq(-2*Im(om3),2*Im(om3),0.001)
# lin3 = -i*Im(om3) + seq(Re(om1)/2,3*Re(2*om1)/4,0.001) ### WITH k = lin1, lin2
# 
# lines(lin1,col="black",lwd=2)
# lines(lin2,col="black",lwd=2)
# lines(lin3,col="black",lwd=2)
# 
# spec_plotter(band1 = lin1,band2 = lin2,g2 = 0,g3 = 12.82538,equi=T,xlim=c(-2,5),ylim=c(-2,2))

## Equianharmonic k-plane alone
pdf(file = "equi_plane.pdf",family = "CM Roman",width = 6,height = 6)
par(mfrow=c(1,1),pty="s")

g2= 0;g3=1
test2 = cond_solver(g2=g2,g3=g3,dtx=0.001,pts=F,equi=T,status=T)

lin = test2[which(abs(Re(test2$x + test2$y) - Re(om1))<0.0001),]
curve = test2[which(abs(Re(test2$x + test2$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]
curve2 = curve[which(Im(curve$x + curve$y) <=0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve2 = curve2$x + curve2$y

lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)

dev.off()
embed_fonts("equi_plane.pdf", outfile="equi_plane_embed.pdf")
##

## Equianharmonic k-plane and spectrum side-by-side
pdf(file = "equi_comp.pdf",family = "CM Roman",width = 12,height = 6)
par(mfrow=c(1,2),pty="s")

cond_solver(g2=0,g3=1,dtx = 0.05,pts=F,equi=T)
lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)
spec_plotter(band1 = lin,band2 = curve1,g2 = g2,g3 = g3,equi=T,xlim=c(-2,5),ylim=c(-2,2))

dev.off()
embed_fonts("equi_comp.pdf", outfile="equi_comp_embed.pdf")
##

#### cond_solver takes exponentially longer to run for smaller dtx, therefore save results as res2.RData
# load("C:/Users/Will/Dropbox/Papers_for_publication/res2.RData") # Work
# load("C:/Users/William/Dropbox/Papers_for_publication/res2.RData") # Home

####
########################################################################################################

########################################################################################################
#### RHOMBIC COMPARISON

## Rhombic k-plane comparison
pdf(file = "rhombic_comp.pdf",family = "CM Roman",width = 12,height = 12)
par(mfrow=c(2,2),pty="s")
# par(mfrow=c(2,2))

# 1
# c1 = cond_solver(g2= -1,g3= -1,dtx=0.001,pts=F,equi=F,status=T)
cond_solver(g2= -1,g3= -1,dtx=0.5,pts=F,equi=F,status=T)

lin = c1[which(abs(Re(c1$x + c1$y) - Re(om1))<0.01),]
curve = c1[which(abs(Re(c1$x + c1$y) - Re(om1))>=0.01),]
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
cond_solver(g2= -1,g3= 1,dtx=0.5,pts=F,equi=F,status=T)

lin = c4[which(abs(Re(c4$x + c4$y) - Re(om1))<0.0001),]
curve = c4[which(abs(Re(c4$x + c4$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]
curve2 = curve[which(Im(curve$x + curve$y) <=0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve2 = curve2$x + curve2$y

lines(lin,lwd=2)
lines(sort(curve1),col="black",lwd=2)
lines(sort(curve2),col="black",lwd=2)

dev.off()
embed_fonts("rhombic_comp.pdf", outfile="rhombic_comp_embed.pdf")
##

####

## Rhombic spectrum comparison
pdf(file = "rhombic_spec.pdf",family = "CM Roman",width = 12,height = 12)
par(mfrow=c(2,2),pty="s")

# 1
cond_solver(g2= -1,g3= -1,dtx=0.5,pts=F,equi=F,status=T,plot=F)

lin = c1[which(abs(Re(c1$x + c1$y) - Re(om1))<0.01),]
curve = c1[which(abs(Re(c1$x + c1$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)

spec_plotter(band1 = lin,band2 = curve1,g2 = -1,g3 = -1,equi=F,xlim=c(-3,7),ylim=c(-5,5))

# 2
cond_solver(g2= -1,g3= -0.4748432462,dtx=0.5,pts=F,equi=F,status=T,plot=F,title=F)
lin = c2[which(abs(Re(c2$x + c2$y) - Re(om1))<0.01),]
curve = c2[which(abs(Re(c2$x + c2$y) - Re(om1))>=0.01),]
curve1 = curve[which(Re(curve$x + curve$y) > Re(om1) & Im(curve$x + curve$y) > 0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y
curve1a = Conj(rev(curve1))
curve1 = c(curve1a,curve1)

spec_plotter(band1 = lin,band2 = curve1,g2 = -1,g3 = -0.4747,equi=F,xlim=c(-3,7),ylim=c(-5,5),title=F)
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
cond_solver(g2= -1,g3= 1,dtx=0.5,pts=F,equi=F,status=T,plot=F)
lin = c4[which(abs(Re(c4$x + c4$y) - Re(om1))<0.0001),]
curve = c4[which(abs(Re(c4$x + c4$y) - Re(om1))>=0.0001),]
curve1 = curve[which(Im(curve$x + curve$y) >0),]

lin = lin$x + lin$y
curve1 = curve1$x + curve1$y

spec_plotter(band1 = lin,band2 = curve1,g2 = -1,g3 = 1,equi=F,xlim=c(-3,7),ylim=c(-5,5))

dev.off()
embed_fonts("rhombic_spec.pdf", outfile="rhombic_spec_embed.pdf")
##

####
########################################################################################################


om1
exp(pi*i/3)*om1

om3

12.82538

pt = gamma(1/3)^6/(2^(14/3)*pi^2) # Proposition 5
points(x=-pt,y=0,col="red",cex=1.3,pch=19)

endpt1 = pt*exp(pi*i/3)
endpt2 = pt*exp(-pi*i/3)
points(endpt1,col="blue",cex=1.3,pch=19)
points(endpt2,col="blue",cex=1.3,pch=19)

tan1 = endpt1 + (2*sqrt(3)*pt/pi - exp(pi*i/3))^2*seq(-1,1,0.1)
tan2 = endpt2 + (2*sqrt(3)*pt/pi - exp(-pi*i/3))^2*seq(-1,1,0.1)
lines(x = Re(tan1),y= Im(tan1))
lines(x = Re(tan2),y= Im(tan2))


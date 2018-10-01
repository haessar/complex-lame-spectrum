x   y    2x-2(x-1)
1   2    2
1.5 1    
2   0    2


x = seq(0.5,1.5,0.1)
y = 2 + (Im(2*om1 - 2*om3))*x
b1 = data.frame("x" = x,"y"=y)
lines(col="red",x,y)

x = seq(1,2,0.1)
y = 4 + (Im(2*om1 - 2*om3))*x
b2 = data.frame("x" = x,"y"=y)
lines(col="red",x,y)

x = seq(0,1,0.1)
y = (Im(2*om1 - 2*om3))*x
b3 = data.frame("x" = x,"y"=y)
lines(col="red",x,y)

y = seq(-2,2,0.1)
x = rep(1,length(y)) 
lines(x,y)
b4 = data.frame("x" = x,"y" = y)

x = seq(0.5,1.5,0.1)
y = rep(1,length(x))
lines(x,y)
b5 = data.frame("x" = x,"y" = y)






spec_plotter(band1 = b4,band2 = b5,om1 = om1,om3=om3,xlim = c(0,10),ylim=c(-2,0))








spec_plotter = function(band1,band2,om1,om3,ylim=c(-5,5),xlim=c(-0.5,5)){
  
  om1 <<- om1#half.periods(g=c(g2,g3))[1] # Set omega_1, omega_3 as global variables
  om3 <<- om3#abs(Re(half.periods(g=c(g2,g3))[2])) + abs(Im(half.periods(g=c(g2,g3))[2]))*i
  om = om1 #+ om3
  
  g2 = Re(g2.fun(c(om1,om3)))
  g3 = Re(g3.fun(c(om1,om3)))
  
  # lines(inf_band,col="red",lwd=1.5);lines(fin_band1,col="red",lwd=1.5);lines(fin_band2,col="red",lwd=1.5) # Plot bands onto a-b plane.
  
  k <<- function(a,b){ # Define the complex spectral parameter k
    (a + i*b)
  }
  
  #### Plot spectral bands
  k1 = k(band1$x,band1$y)
  k2 = k(band2$x,band2$y)
  
  plot(-om^2*P(k1,g=c(g2,g3)),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
       main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
  abline(v=0,h=0,lty=2,col="gray")
  
  lines(-om^2*P(k1,g=c(g2,g3)),col="black",lwd=2)
  lines(-om^2*P(k2,g=c(g2,g3)),col="black",lwd=2)
  

}
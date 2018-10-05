rm(list=ls())
library(elliptic)

i = sqrt(as.complex(-1))

cond_solver = function(g2 = 4, g3 = 1, dtx = 0.1, dty = 0.1,x_bound = c(-1,1),y_bound=c(-1,1),plot=T,pts=F,lns=F){
  
  om1 <<- half.periods(g=c(g2,g3))[1] # Set omega_1, omega_3 as global variables
  om3 <<- abs(Re(half.periods(g=c(g2,g3))[2])) + abs(Im(half.periods(g=c(g2,g3))[2]))*i
  om = om1 + om3
  
  g2 <<- g2 # Set g_2, g_3 as global variables
  g3 <<- g3
  
  fn2 <- function(a,b){ # Function that we want to be equal to zero
    f <- zeta(z=om,g=c(g2,g3))*(a*om1 + b*om3) - zeta(z=(a*om1 + b*om3),g=c(g2,g3))*om
    as.numeric(Re(f))#,Im(f))
  }
  
  if (plot==T){
#         par(mfrow=c(1,1))
#         latplot(p = c(om1,om3),1,xlab="a",ylab="b") # Plot a-b plane.
    plot(1:10,sin(1:10),type="n",xlim=c(-1,1),ylim=c(-1,1),xlab=expression(italic("a")),ylab=expression(italic("b")),bty="L",
         cex.lab=1.7,main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.main=2,cex.axis=1.5)
    abline(v=0,h=0,lty=2,col="gray")
  }
  
  count = 0
  re = seq(x_bound[1],x_bound[2],dtx) # Increments of x-axis
  im = seq(y_bound[1],y_bound[2],dty) # Increments of y-axis
  for (a in 1:length(re)){ # Move along x-axis incrementally
    ofv = 0
    if (abs(re[a])<0.00001 ){#& any(abs(im)<0.00001)){
      count = count + 1
      next # Avoid pole
    }
    
    if (a>1){
      eps = dty*5
      low = ifelse( ypts[a-1-count] - eps < y_bound[1] , y_bound[1], ypts[a-1-count] - eps)
      high = ifelse( ypts[a-1-count] + eps > y_bound[2] , y_bound[2], ypts[a-1-count] + eps)
      im = seq(low,high,dty/10)
    }
    for (b in 1:length(im)){ # For some x-value, determine fn2 value for all y-axis increments
      ofv[b] = abs(fn2(re[a],im[b]))
    }
    s_ofv = sort(ofv) # Sort resulting values in ascending order
    # s_ofv = round(s_ofv,3) # Round to 3 decimals
    # ofv = round(ofv,3)
    ypt = s_ofv[1] # For chosen y-values, take the first
    ypt = im[which(ofv %in% ypt)]
    xpt = rep(re[a],length(ypt)) # Repeat x-value length(ypt) times to match length of ypt vector.
    if (pts==T) points(xpt,ypt,col="blue") # Plot the resulting points on complex a-b plane.
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

spec_plotter = function(quad1,quad2,quad3,quad4,g2,g3,ylim=c(-5,5),xlim=c(-0.5,5)){
  
  om1 <<- half.periods(g=c(g2,g3))[1] # Set omega_1, omega_3 as global variables
  om3 <<- abs(Re(half.periods(g=c(g2,g3))[2])) + abs(Im(half.periods(g=c(g2,g3))[2]))*i
  om = om1 + om3
  
  # lines(inf_band,col="red",lwd=1.5);lines(fin_band1,col="red",lwd=1.5);lines(fin_band2,col="red",lwd=1.5) # Plot bands onto a-b plane.
  
  k <<- function(a,b){ # Define the complex spectral parameter k
    (a*om1 + b*om3)
  }
  
  #### Plot spectral bands
  k1 = k(quad1$x,quad1$y)
  k2 = k(quad2$x,quad2$y)
  k3 = k(quad3$x,quad3$y)
  k4 = k(quad4$x,quad4$y)
  
  plot(-om^2*P(k1,g=c(g2,g3)),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
       main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
  abline(v=0,h=0,lty=2,col="gray")
  
  lines(-om^2*P(k1,g=c(g2,g3)),col="black",lwd=2)
  lines(-om^2*P(k2,g=c(g2,g3)),col="black",lwd=2)
  lines(-om^2*P(k3,g=c(g2,g3)),col="red",lwd=2)
  lines(-om^2*P(k4,g=c(g2,g3)),col="red",lwd=2)
  
}


incr = 0.005
par(pty = "s")
q1=cond_solver(g2 = 0,g3 = 1,dtx = incr,dty = incr,x_bound = c(-1,0),y_bound = c(-1 + incr,0),plot = T,pts=T)
q2=cond_solver(g2 = 0,g3 = 1,dtx = incr,dty = incr,x_bound = c(-1,0),y_bound = c(0,1),plot = F,pts=T)
q3=cond_solver(g2 = 0,g3 = 1,dtx = incr,dty = incr,x_bound = c(incr,1),y_bound = c(incr,1-incr),plot = F,pts=T)
q4=cond_solver(g2 = 0,g3 = 1,dtx = incr,dty = incr,x_bound = c(incr,1),y_bound = c(-1 + incr,0),plot = F,pts=T)

spec_plotter(quad1 = q1,quad2 = q2,quad3 = q3,quad4 = q4,g2 = 0,g3 = 1,xlim = c(-1,30),ylim=c(-1,9))

pnt = 3*gamma(1/3)^6 / (2^(17/3)*pi^2)
pnt
om1*exp(i*pi/3)
points(pnt)

-om^2*e1e2e3(g=c(0,1))["e2"]


points(-pnt*exp(2*pi*i/3))

plot(1:10,xlim=c(-2,2),ylim=c(-2,2),type='n')
lines(om1*(seq(-1,1,0.1)))
lines(om3*(seq(-1,1,0.1)))
lines(om*(seq(-1,1,0.1)))


cond_solver(g2 = 4,g3 = 1,dtx = 0.05,dty = 0.05,x_bound = c(-1,0),y_bound = c(-1 + incr,0),plot = T,pts=T)


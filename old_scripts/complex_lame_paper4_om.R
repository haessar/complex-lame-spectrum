rm(list=ls())
library(elliptic)
library(extrafont)
setwd(dir = "C:/Users/Will/Dropbox/Papers_for_publication/")
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")


### ONE TIME ONLY:
# font_import()
# font_install('fontcm')
# loadfonts()


i = sqrt(as.complex(-1))

### Go with (g2,g3) = (32,1), (16,1), (8,1), (4,1)
### y=0.6 should be boundary

cond_solver = function(om1 = 1, om3 = i, dtx = 0.1, dty = 0.1,x_bound = c(-1,1),y_bound=c(-1,1),plot=T,pts=F,lns=F){
  
  om1 <<- om1
  om3 <<- om3
  om = om1 + om3
  
  
  fn2 <- function(a,b){ # Function that we want to be equal to zero
    f <- zeta(z=om,Omega=c(om1,om3))*(a*om1 + b*om3) - zeta(z=(a*om1 + b*om3),Omega=c(om1,om3))*om
    as.numeric(Re(f))#,Im(f))
  }
  
  if (plot==T){
    #     par(mfrow=c(1,1))
    #     latplot(p = c(1,i),1,xlab="a",ylab="b") # Plot a-b plane.
    plot(1:10,sin(1:10),type="n",xlim=c(-1,1),ylim=c(-1,1),xlab=expression(italic("a")),ylab=expression(italic("b")),bty="L",
         cex.lab=1.7,main = bquote("omega"[1] ~ "=" ~ .(om1) ~ ", " ~ "omega"[3] ~ "=" ~ .(om3)),cex.main=2,cex.axis=1.5)
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


### (g2,g3) = (32,1)
# fin1 = cond_solver(g2 = 32,g3 = 1, dtx = 0.05, dty = 0.05,x_bound = c(-1,0),y_bound = c(0.6,1))
# fin2 = cond_solver(g2 = 32,g3 = 1, dtx = 0.05, dty = 0.05,x_bound = c(0,1),y_bound = c(-1,-0.6),plot=F)
# inf = cond_solver(g2 = 32,g3 = 1, dtx = 0.05, dty = 0.05,x_bound = c(-1,1),y_bound = c(-0.6,0.6),plot=F)


###########################################################################################
###########################################################################################

# inf_band = inf; fin_band1 = fin1; fin_band2 = fin2

spec_plotter = function(fin_band1,fin_band2,inf_band,om1,om3,ylim=c(-5,5),xlim=c(-0.5,5)){
  
  om1 <<- om1
  om3 <<- om3
  om = om1 + om3
  
  # lines(inf_band,col="red",lwd=1.5);lines(fin_band1,col="red",lwd=1.5);lines(fin_band2,col="red",lwd=1.5) # Plot bands onto a-b plane.
  
  k <<- function(a,b){ # Define the complex spectral parameter k
    (a*om1 + b*om3)
  }
  
  #### Plot spectral bands
  k_inf = k(inf_band$x,inf_band$y)
  k_f1 = k(fin_band1$x,fin_band1$y)
  k_f2 = k(fin_band2$x,fin_band2$y)
  
  plot(-om^2*P(k_inf,Omega = c(om1,om3)),ylim=ylim,xlim=xlim,type = 'n',xlab= "Re",ylab="Im",axes=T,
       main = bquote("omega"[1] ~ "=" ~ .(om1) ~ ", " ~ "omega"[3] ~ "=" ~ .(om3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
  abline(v=0,h=0,lty=2,col="gray")
  
  lines(-om^2*P(k_inf,Omega = c(om1,om3)),col="black",lwd=2)
  lines(-om^2*P(k_f1,Omega = c(om1,om3)),col="black",lwd=2)
  lines(-om^2*P(k_f2,Omega = c(om1,om3)),col="black",lwd=2)
  
}

# spec_plotter(fin_band1 = fin1,fin_band2 = fin2,inf_band = inf,g2 = 32,g3=1)



###########################################################################################
###########################################################################################

wrapper = function(om1,om3,incr,plot_which=2,bnd = 0.6){
  
  plt = F
  if (plot_which==1) plt = T
  
  fin1 = cond_solver(om1 = om1, om3 = om3,dtx = incr, dty = incr, x_bound = c(-1,0), y_bound = c(bnd,1),plot = plt,lns = T)
  fin2 = -fin1
  fin2 = fin2[order(fin2$x),]
  lines(fin2,col="black",lwd=2)
  inf = cond_solver(om1 = om1, om3 = om3,dtx = incr, dty = incr,x_bound = c(-1,0), y_bound = c(-bnd,bnd),plot = F,lns = T)
  mirror_inf = -inf
  inf = rbind(inf,mirror_inf)
  inf = inf[order(inf$x),]
  lines(inf,col="black",lwd=2)
  
  if (plot_which==2){
    spec_plotter(fin_band1 = fin1,fin_band2 = fin2,inf_band = inf,om1 = om1, om3 = om3)
  }
  
  out = list()
  out[[1]] = fin1
  out[[2]] = fin2
  out[[3]] = inf
  
  return(out)
}

###########################################################################################
###########################################################################################

par(pty="s")
keep = wrapper(om1 = 1, om3 = 10*i,incr = 0.05,plot_which = 1,bnd = 0.642429)
spec_plotter(keep[[1]],keep[[2]],keep[[3]],om1 = 1,om3 = 100*i,ylim=c(-200,200),xlim=c(-10000,2000))

inf = seq(-1,1,0.1)
fin= seq(-1,0,0.1) + rep(om3,11)

-om^2*P(z = inf,Omega = c(om1,om3))

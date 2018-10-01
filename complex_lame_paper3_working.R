rm(list=ls())
library(elliptic)

i = sqrt(as.complex(-1))

cond_solver = function(g2 = 4, g3 = 1, dtx = 0.1, dty = 0.1,num_pts = 2){
  
  om1 <<- Re(half.periods(g=c(g2,g3))[1]) # Set omega_1, omega_3 as global variables
  om3 <<- Im(half.periods(g=c(g2,g3))[2])
  om = om1 + om3*i
  
  g2 <<- g2 # Set g_2, g_3 as global variables
  g3 <<- g3
  
  fn2 <- function(a,b){ # Function that we want to be equal to zero
    f <- zeta(z=om,g=c(g2,g3))*(a*om1 + i*b*om3) - zeta(z=(a*om1 + i*b*om3),g=c(g2,g3))*om
    as.numeric(Re(f))#,Im(f))
  }
  
  par(mfrow=c(1,1))
  latplot(p = c(1,i),1,xlab="a",ylab="b") # Plot a-b plane.
  
  re = seq(-1,1,dtx) # Increments of x-axis
  im = seq(-1,1,dty) # Increments of y-axis
  for (a in 1:length(re)){ # Move along x-axis incrementally
    ofv = 0
    for (b in 1:length(im)){ # For some x-value, determine fn2 value for all y-axis increments
      ofv[b] = abs(fn2(re[a],im[b]))
    }
    s_ofv = sort(ofv) # Sort resulting values in ascending order
    s_ofv = round(s_ofv,3) # Round to 3 decimals
    ofv = round(ofv,3)
    ypt = s_ofv[1:num_pts] # For chosen y-values, take the first 'num_pts' of them
    ypt = im[which(ofv %in% ypt)]
    xpt = rep(re[a],length(ypt)) # Repeat x-value length(ypt) times to match length of ypt vector.
    points(xpt,ypt,col="blue") # Plot the resulting points on complex a-b plane.
    if (a==1){
      ypts = ypt
      xpts = xpt
    } else {
      ypts = c(ypts,ypt)
      xpts = c(xpts,xpt)
    }
  }
  
  k <<- function(a,b){ # Define the complex spectral parameter k
    (a*om1 + i*b*om3)
  }
  
  out = list() # Output k values, and the corresponding a-b plane values
  out[[1]] = k(xpts,ypts) 
  out[[2]] = data.frame(xpts,ypts)
  return(out)
}

###########################################################################################
###########################################################################################

out = cond_solver(g2 = 4,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out2 = cond_solver(g2 = 1,g3 = 0, dtx = 0.05, dty = 0.05,num_pts = 2)
out3 = cond_solver(g2 = 4,g3 = 1, dtx = 0.001, dty = 0.001,num_pts = 2)
out4 = cond_solver(g2 = 1,g3 = 0, dtx = 0.005, dty = 0.005,num_pts = 2)
out5 = cond_solver(g2 = 24,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out6 = cond_solver(g2 = 12,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out7 = cond_solver(g2 = 8,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out8 = cond_solver(g2 = 16,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out9 = cond_solver(g2 = 32,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out10 = cond_solver(g2 = 4,g3 = 1, dtx = 0.1, dty = 0.1,num_pts = 2)
out11 = cond_solver(g2 = 3.3,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out12 = cond_solver(g2 = 64,g3 = 1, dtx = 0.005, dty = 0.005,num_pts = 2)
out13 = cond_solver(g2 = 0,g3 = 1, dtx = 0.05, dty = 0.05,num_pts = 2)


#### cond_solver takes exponentially longer to run for smaller dtx,dty, therefore save results as res.RData
# load("C:/Users/Will/Dropbox/Papers_for_publication/res.RData") # Work
# load("C:/Users/William/Dropbox/Papers_for_publication/res.RData") # Home

###############################################################
## Set the following to the output parameters that are required:
# k_plane = out6[[1]]
ab_plane = out[[2]]
g2 = 4;g3 = 1
###############################################################

#### RUN ALL THE FOLLOWING TO PLOT SPECTRUM

spec_plotter = function(ab_plane,g2,g3,smooth = F,span=0.075){
  om1 <<- Re(half.periods(g=c(g2,g3))[1]) # Set omega_1, omega_3 as global variables
  om3 <<- Im(half.periods(g=c(g2,g3))[2])
  om = om1 + om3*i
  
  ### RE-plot the a-b lattice and points
#   latplot(p = c(1,i),1,xlab="a",ylab="b") # Plot a-b plane.
#   xpt = ab_plane$xpts; ypt = ab_plane$ypts
#   points(xpt,ypt,col="blue") # Plot the resulting points on complex a-b plane.
  
  #### Establish bands
  x_axis = unique(ab_plane$xpts)
  mat = as.data.frame(matrix(ncol=2,nrow=length(x_axis)))
  colnames(mat) = c("x","y")
  fin_band1 = fin_band2 = inf_band = mat
  f1 = f2 = inf = 0
  for (j in 1:length(x_axis)){
  # for (j in 1:58){
    x = ab_plane$xpts[which(ab_plane$xpts == x_axis[j])]
    y = ab_plane$ypts[which(ab_plane$xpts == x_axis[j])]
    if (j==1){
      fin_band1$x[j] = fin_band2$x[j] = inf_band$x[j] = unique(x)
      fin_band1$y[j] = y[which(y == 1)]
      fin_band2$y[j] = y[which(y == -1)]
      inf_band$y[j] = y[which(y == 0)]
    } else {
      if (any(abs(y - fin_band1$y[j-1-f1])<0.5)){
        fin_band1$x[j-f1] = unique(x)
        fin_band1$y[j-f1] = y[which(abs(y - fin_band1$y[j-1-f1]) == min(abs(y - fin_band1$y[j-1-f1])))]
      } else f1 = f1 + 1
      if (any(abs(y - fin_band2$y[j-1-f2])<0.5)){
        fin_band2$x[j-f2] = unique(x)
        fin_band2$y[j-f2] = y[which(abs(y - fin_band2$y[j-1-f2]) == min(abs(y - fin_band2$y[j-1-f2])))]
      } else f2 = f2 + 1
      #### PROBLEM WITH INF BAND LINE SKIPPING UP TO FINITE BAND.
      # DETERMINE IF +VE OR -VE GRADIENT.
      len = length(which(is.na(inf_band$y)==F))
      if (len>5){
        band = tail(inf_band$y[which(is.na(inf_band$y)==F)],5) # Determine last 3 infinite band y-axis values
        if (all(band == cummin(band))){ # IF negative gradient
          left_max = max(inf_band$y[which(is.na(inf_band$y)==F)]) # Determine max point of inf-band for bound
          if (any(abs(y)<left_max)){
            inf_band$x[j-inf] = unique(x)
            inf_band$y[j-inf] = y[which(abs(y - inf_band$y[j-1-inf]) == min(abs(y - inf_band$y[j-1-inf])))]
          } else inf = inf + 1
        } else {
          if (any(abs(y - inf_band$y[j-1-inf])<0.25)){
            inf_band$x[j-inf] = unique(x)
            inf_band$y[j-inf] = y[which(abs(y - inf_band$y[j-1-inf]) == min(abs(y - inf_band$y[j-1-inf])))]
          } else inf = inf + 1
        }
      } else {
        if (any(abs(y - inf_band$y[j-1-inf])<0.25)){
          inf_band$x[j-inf] = unique(x)
          inf_band$y[j-inf] = y[which(abs(y - inf_band$y[j-1-inf]) == min(abs(y - inf_band$y[j-1-inf])))]
        } else inf = inf + 1
      }
    }
  }
  
  # lines(inf_band,col="red",lwd=1.5);lines(fin_band1,col="red",lwd=1.5);lines(fin_band2,col="red",lwd=1.5) # Plot bands onto a-b plane.
  
  if (smooth == T){ # apply regression line to fit points and create smoother curve?
    
    # Smooth inf_band
    x = inf_band$x[which(is.na(inf_band$x)==F)]
    y = predict(loess(y ~ x, inf_band, span = span))
    inf_band = data.frame(x,y)
    
    # Smooth fin_band1
    fin_band1[which(fin_band1$x > 0.9),] = NA # Prevent erroneous line being drawn
    x = fin_band1$x[which(is.na(fin_band1$x)==F)]
    y = predict(loess(y ~ x, fin_band1, span = span))
    fin_band1 = data.frame(x,y)
    
    # Smooth fin_band2
    fin_band2[which(fin_band2$x < -0.9),] = NA # Prevent erroneous line being drawn
    x = fin_band2$x[which(is.na(fin_band2$x)==F)]
    y = predict(loess(y ~ x, fin_band2, span = span))
    fin_band2 = data.frame(x,y)
    
    # lines(inf_band); lines(fin_band1); lines(fin_band2)
    
    k_inf = k(inf_band$x,inf_band$y)
    k_f1 = k(fin_band1$x,fin_band1$y)
    k_f2 = k(fin_band2$x,fin_band2$y)
    
  } else {
    
    #### Plot spectral bands
    k_inf = k(inf_band$x,inf_band$y)
    fin_band1[which(fin_band1$x > 0.9),] = NA # Prevent erroneous line being drawn
    k_f1 = k(fin_band1$x,fin_band1$y)
    fin_band2[which(fin_band2$x < -0.9),] = NA # Prevent erroneous line being drawn
    k_f2 = k(fin_band2$x,fin_band2$y)
    
  }
  
  # plot(-P(k_inf,g=c(g2,g3)),ylim=c(-3,0.2),xlim=c(-1.2,1),type = 'n',xlab= "Re",ylab="Im",axes=T)
  plot(-P(k_inf,g=c(g2,g3)),ylim=c(-5,0.1),xlim=c(-2.5,2.5),type = 'n',xlab= "Re",ylab="Im",axes=T,
       main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)))
  
  lines(-P(k_inf,g=c(g2,g3)),col="blue",lwd=2)
  lines(-P(k_f1,g=c(g2,g3)),col="blue",lwd=2)
  lines(-P(k_f2,g=c(g2,g3)),col="blue",lwd=2)
  abline(v=0,lty=2);abline(h=0,lty=2)
  
  #### Proposition 2.1
  lines((-100:100)*Conj(om)^2,col="red")
}

par(mfrow=c(2,2))
spec_plotter(ab_plane = out4[[2]],g2 = 1,g3 = 0)
spec_plotter(ab_plane = out5[[2]],g2 = 24,g3 = 1)
spec_plotter(ab_plane = out6[[2]],g2 = 12,g3 = 1)
spec_plotter(ab_plane = out[[2]], g2 = 4,g3 = 1)

par(mfrow=c(1,1))
spec_plotter(ab_plane = out10[[2]],g2 = 4, g3 = 1,smooth=F)

par(mfrow=c(2,2))
span=0.075
spec_plotter(ab_plane = out9[[2]],g2 = 32,g3 = 1,smooth = T,span)
spec_plotter(ab_plane = out8[[2]],g2 = 16,g3 = 1,smooth = T,span)
spec_plotter(ab_plane = out7[[2]],g2 = 8, g3 = 1,smooth = T,span)
spec_plotter(ab_plane = out[[2]], g2 = 4, g3 = 1,smooth = T,span)

par(mfrow=c(1,1))
spec_plotter(ab_plane = out11[[2]],g2 = 3.3, g3 = 1,smooth=F)

par(mfrow=c(2,2))
spec_plotter(ab_plane = out4[[2]],g2 = 1,g3 = 0)
spec_plotter(ab_plane = out12[[2]], g2 = 64, g3 = 1,smooth = T,span)
spec_plotter(ab_plane = out8[[2]],g2 = 16,g3 = 1,smooth = T,span)
spec_plotter(ab_plane = out[[2]], g2 = 4, g3 = 1,smooth = T,span)


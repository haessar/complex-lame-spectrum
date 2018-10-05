
library(elliptic)


spec_plotter = function(g2 = 4, g3 = 1, dtx = 0.1, dty = 0.1,num_pts = 2){
  
  i = sqrt(as.complex(-1))
  
  # g2 = 4
  # g3 = 1
  
  om1 = Re(half.periods(g=c(g2,g3))[1])
  om3 = Im(half.periods(g=c(g2,g3))[2])
  om= om1 + om3*i
  
  fn2 <- function(a,b){
    f <- zeta(z=om,g=c(g2,g3))*(a*om1 + i*b*om3) - zeta(z=(a*om1 + i*b*om3),g=c(g2,g3))*om
    as.numeric(Re(f))#,Im(f))
  }
  
  latplot(p = c(1,i),1)
  
  
  re = seq(-1,1,dtx)
  im = seq(-1,1,dty)
  for (a in 1:length(re)){
    ofv = 0
    for (b in 1:length(im)){
      ofv[b] = abs(fn2(re[a],im[b]))
    }
    s_ofv = sort(ofv)
    s_ofv = round(s_ofv,3)
    ofv = round(ofv,3)
    ypt = s_ofv[1:num_pts]
    ypt = im[which(ofv %in% ypt)]
    # ypt = im[which(ofv == min(ofv))]
    xpt = rep(re[a],length(ypt))
    points(xpt,ypt,col="blue")
    if (a==1){
      ypts = ypt
      xpts = xpt
    } else {
      ypts = c(ypts,ypt)
      xpts = c(xpts,xpt)
    }
  } 
  
  
  k = function(a,b){
    (a*om1 + i*b*om3)
  }
  
  # plot(k(xpts,ypts))
  
  # plot(-P(k(xpts,ypts),g=c(g2,g3)),ylim=c(-10,1),xlim=c(-2,2))

  out = list()
  out[["k"]] = k(xpts,ypts)
  out[["k_plane"]] = data.frame(xpts,ypts)
  return(out)
}

g2 = 4
g3 = 1

out = spec_plotter(g2 = g2,g3 =g3, dtx = 0.04, dty = 0.04,num_pts = 2)

k = out[[1]]
k_plane = out[[2]]

#### Establish 2 bands
x_axis = unique(k_plane$xpts)
mat = as.data.frame(matrix(ncol=2,nrow=length(x_axis)))
colnames(mat) = c("x","y")
fin_band1 = fin_band2 = inf_band = mat
f1 = f2 = inf = 0
for (i in 1:length(x_axis)){
# for (i in 1:20){
  x = k_plane$xpts[which(k_plane$xpts == x_axis[i])]
  y = k_plane$ypts[which(k_plane$xpts == x_axis[i])]
  if (i==1){
    fin_band1$x[i] = fin_band2$x[i] = inf_band$x[i] = unique(x)
    fin_band1$y[i] = y[which(y == 1)]
    fin_band2$y[i] = y[which(y == -1)]
    inf_band$y[i] = y[which(y == 0)]
  } else {
#     if (length(x)==3){
#       fin_band1$x[i-f1] = fin_band2$x[i-f2] = inf_band$x[i-inf] = unique(x)
#       fin_band1$y[i-f1] = y[which(abs(y - fin_band1$y[i-1-f1]) == min(abs(y - fin_band1$y[i-1-f1])))]
#       fin_band2$y[i-f2] = y[which(abs(y - fin_band2$y[i-1-f2]) == min(abs(y - fin_band2$y[i-1-f2])))]
#       inf_band$y[i-inf] = y[which(abs(y - inf_band$y[i-1-inf]) == min(abs(y - inf_band$y[i-1-inf])))]
#     } else if (length(x)==2){
      if (any(abs(y - fin_band1$y[i-1-f1])<0.5)){
        fin_band1$x[i-f1] = unique(x)
        fin_band1$y[i-f1] = y[which(abs(y - fin_band1$y[i-1-f1]) == min(abs(y - fin_band1$y[i-1-f1])))]
      } else f1 = f1 + 1
      if (any(abs(y - fin_band2$y[i-1-f2])<0.5)){
        fin_band2$x[i-f2] = unique(x)
        fin_band2$y[i-f2] = y[which(abs(y - fin_band2$y[i-1-f2]) == min(abs(y - fin_band2$y[i-1-f2])))]
      } else f2 = f2 + 1
      if (any(abs(y - inf_band$y[i-1-inf])<0.5)){
        inf_band$x[i-inf] = unique(x)
        inf_band$y[i-inf] = y[which(abs(y - inf_band$y[i-1-inf]) == min(abs(y - inf_band$y[i-1-inf])))]
      } else inf = inf + 1
    }
  }
# }

lines(inf_band);lines(fin_band1);lines(fin_band2)

k = function(a,b){
  (a*om1 + i*b*om3)
}

k_inf = k(inf_band$x,inf_band$y)
k_f1 = k(fin_band1$x,fin_band1$y)
k_f2 = k(fin_band2$x,fin_band2$y)


plot(-P(k_inf,g=c(g2,g3)),ylim=c(-10,1),xlim=c(-2,2),type = 'n')
lines(-P(k_inf,g=c(g2,g3)))
lines(-P(k_f1,g=c(g2,g3)))
lines(-P(k_f2,g=c(g2,g3)))


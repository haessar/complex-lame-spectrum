# install.packages("elliptic")
# install.packages("BB")

library(elliptic)
library(BB)
library(rootSolve)

i = sqrt(as.complex(-1))

g2 = 4
g3 = 1

om1 = Re(half.periods(g=c(g2,g3))[1])
om3 = Im(half.periods(g=c(g2,g3))[2])
om= om1 + om3*i

eq = function(a,b){
  out = zeta(z=om,g=c(g2,g3))*(a*om1 + i*b*om3) - zeta(z=(a*om1 + i*b*om3),g=c(g2,g3))*om
  return(out)
}

latplot(p = c(1,i),1)

bs=1
fixb = function(a){
  as.complex(eq(a,b=bs))
}
fixb(-0.3)

as=1
fixa = function(b){
  as.complex(eq(a=as,b))
}
fixa(-1)

uniroot(fixb,interval = c(-1,1))



fn1 <- function(x, a){ 
  z <- x[1] + 1i * x[2] 
  f <- exp(z) + a 
  c(Re(f), Im(f)) 
} 

BBsolve(par=c(1,1), fn=fn1, a=1)

fn2 <- function(a,b){
  f <- zeta(z=om,g=c(g2,g3))*(a*om1 + i*b*om3) - zeta(z=(a*om1 + i*b*om3),g=c(g2,g3))*om
  as.numeric(Re(f))#,Im(f))
}

re = seq(-1,1,0.1)
im = seq(-1,1,0.1)
for (a in 1:length(re)){
  ofv = 0
  for (b in 1:length(im)){
    ofv[b] = fn2(re[a],im[b])
  }
  ypt = im[which(abs(ofv) == min(abs(ofv)))]
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

ypts
xpts

k = function(a,b){
  (a*om1 + i*b*om3)
}

k(xpts,ypts)

plot(P(k(xpts,ypts),g=c(g2,g3)))



lines(xpts,ypts)

a = seq(-1,1,0.1)
for (l in 1:length(a)){
  test = BBsolve(par=0,fn = fn2,a=a[l])$par
  BBoptim(par=0,fn = fn2,a=a[l])
  points(x = test[1],y=test[2])
}




###### contourplot

library(lattice)
x <- seq(pi/4, 5 * pi, length.out = 100)
y <- seq(pi/4, 5 * pi, length.out = 100)
r <- as.vector(sqrt(outer(x^2, y^2, "+")))
grid <- expand.grid(x=x, y=y)
grid$z <- cos(r^2) * exp(-r/(pi^3))
levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = FALSE, region = TRUE)






x <- 1:9; names(x) <- x
# Multiplication & Power Tables
x %o% x
y <- 2:8; names(y) <- paste(y,":", sep = "")
outer(y, x, "^")

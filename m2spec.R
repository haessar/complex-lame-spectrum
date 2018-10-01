rm(list=ls())

library(elliptic)
library(Matrix)
library(extrafont)
# setwd(dir = "C:/Users/William/Dropbox/Papers_for_publication/") # Home
setwd(dir = "C:/Users/Will/Dropbox/Papers_for_publication/") # Work
# setwd(dir = "/home/haezar/Dropbox/Papers_for_publication/") # Laptop
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

i = sqrt(as.complex(-1))


m = 2
rootOrd = 1000 # order of roots of unity considered
steps = 8000 # nr of steps used in the discretization
h = 2*rootOrd/steps
g2 = 0.5
g3 = 1
omega1 = half.periods(g=c(g2,g3))[1]
omega3 = half.periods(g=c(g2,g3))[2]

polyroot(c(-g3, -g2, 0, 4))

omega = omega1 #caution is needed here: even if there is a lattice basis with omega1 real Mathematica may not pick such a basis
z0 = omega3


dim = steps + 1
diags = list(rep(-2, dim), rep(1, dim), rep(1, dim), rep(1, dim), rep(1, dim))
K = -(1/h^2)*bandSparse(dim, k = c(0,-1,1,-steps,steps), diagonals = diags)
K = as.matrix(K)

pot = function(x){
  m*(m + 1)*omega^2*P(omega*(2*rootOrd/steps)*x + z0, Omega = c(omega1, omega3))
}  

P(omega*(2*rootOrd/steps)*1 + z0, g = c(g2,g3))

Vdiag = pot(0:(steps))
V = diag(Vdiag)
# View(V)
# View(K)
# View(V+K)
eigVals = eigen(K + V, symmetric = FALSE,only.values = TRUE)
eigVals_g3p05_2 = eigVals

load("eigVals_g3p05_2.RData")

eigVals = eigVals_g30
# g2=0
# g3=1

test = rev(eigVals[[1]])
# test = test[1:1000]
plot(test, xlim=c(-3,5), ylim=c(-4,4))

### PRODUCE PLOTS
# 
# pdf(file = "m2spec_g30.pdf",family = "CM Roman",width = 12,height = 12)
par(mfrow=c(1,1),pty="s")
plot(test,xlim=c(-3,6), ylim=c(-4,4),type = 'n',xlab= "Re",ylab="Im",axes=T,
     main = bquote("g"[2] ~ "=" ~ .(g2) ~ ", " ~ "g"[3] ~ "=" ~ .(g3)),cex.lab=1.7,cex.main=2,cex.axis=1.5)
abline(v=0,h=0,lty=2,col="gray")

cv1 = test[Re(test)<0 & abs(Im(test))>0.1] 
cv1_tab = data.frame("Re" = Re(cv1), "Im" = Im(cv1))
cv1_tab = cv1_tab[order(cv1_tab$Im,decreasing = T),]
cv1 = NULL
for (j in 1:nrow(cv1_tab)){
  cv1 = c(cv1, cv1_tab$Re[j] + i*cv1_tab$Im[j])
}
lines(cv1,col="black",lwd=3)


cv2 = test[Re(test)<3 & abs(Im(test))<0.001] 
lines(cv2,col="black",lwd=3)

cv3 = test[Re(test)>=3 & abs(Im(test))<0.001]
lines(cv3,col="black",lwd=3)


# dev.off()
# embed_fonts("m2spec_g30.pdf", outfile="m2spec_g30_embed.pdf")

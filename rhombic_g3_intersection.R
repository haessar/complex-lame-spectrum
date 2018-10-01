library(elliptic)
options(digits=10)

i = sqrt(as.complex(-1))


cond = function(g3){
  om1 <- 2*Re(half.periods(g=c(-1,g3))[1]) # + om3  #### half.periods gives om3, om1-om3, respectively
  out = zeta(om1,g = c(-1,g3)) + om1*P(om1,g=c(-1,g3))
  return(as.numeric(Re(out)))
}



optimize(cond,interval = c(-0.475,-0.474),tol=0.0000000001)
optimize(cond,interval = c(-0.48,-0.47))
optimize(cond,interval = c(-0.5,-0.4))


uniroot(cond,interval=c(-1,0),tol = 0.0000000001)


cond(-0.4748432462)

library(elliptic)

i = sqrt(as.complex(-1))

### Conditions for Rhombic lattice:
# om1 real and > 0
# Im(om3) > 0
# Re(om3) = 1/2 om1
# Del < 0 


om1 = 1
om3 = ((1 + i)*om1)/2 # pseudo-lemniscatic case
# om2 = om1 + om3


latplot(p = c(2*om1,2*om3),1,xlab="a",ylab="b",ylim=c(-2,2)) # Plot a-b plane.

a = seq(0,1,0.1)
lines(x=a*2*om1,y=rep(0,length(a)),col="red",lwd=2)
lines(a*2*om3,col="red",lwd=2)
# lines(a*2*om2,col="red",lwd=2)


# rot = (1+i)/2
rot = om3

coord1 = (2*om1 - 2*om3)*rot
# coord2 = om2*rot
coord3 = 2*om3*rot

lines(col="blue",coord1*a,lwd=2)
# lines(col="blue",coord2*a)
lines(col="blue",coord3*a,lwd=2)

g2 = Re(g2.fun(c(coord1,coord3)))
g3 = Re(g3.fun(c(coord1,coord3)))

lom1 = half.periods(g=c(-4,0))[1]
lom3 = half.periods(g=c(-4,0))[2]

lines(col="green",lom1*a,lwd=2)
lines(col="green",lom3*a,lwd=2)

coord1 = 2*lom1*rot
# coord2 = om2*rot
coord3 = 2*lom3*rot

lines(col="blue",coord1*a,lwd=2)
# lines(col="blue",coord2*a)
lines(col="blue",coord3*a,lwd=2)

half.periods(g=c(1,0))


####### Determine the pseudo-lemniscatic g2,g3 values that are analagous to lemniscatic g2=1,g3=0:

lom1 = half.periods(g=c(1,0))[1]
lom3 = half.periods(g=c(1,0))[2]

om3 = lom3/(2*rot)
om1 = lom1/(2*rot) + om3

g2 = Re(g2.fun(c(om1,om3))) # g2 = -4
g3 = Re(g3.fun(c(om1,om3))) # g3 = 0






rot = (1-i)/2

coord1 = (2*om3)*rot
coord2 = (2*om1 - 2*om3)*rot

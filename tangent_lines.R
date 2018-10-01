lam = i*(gamma(1/4)^4)/(16*pi)

points(lam)
points(-lam)

tan = function(t) +lam + (gamma(1/4)^4/(8*pi^2) - i)^2*t
lines(tan(-1:1),col="blue")

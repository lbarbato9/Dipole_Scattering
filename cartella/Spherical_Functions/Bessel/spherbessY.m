
function spherbessY=spherbessY(nu,z)
spherbessY=(sqrt(pi/2)).*( sqrt(1./z)) .*( bessely((nu + 1/2), z));
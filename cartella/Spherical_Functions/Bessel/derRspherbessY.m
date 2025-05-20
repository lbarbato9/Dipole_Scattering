
function RderY = derRspherbessY(nu, z)

RderY = ((-z.*spherbessY(nu + 1, z) + (nu + 1).*spherbessY(nu, z))./z);

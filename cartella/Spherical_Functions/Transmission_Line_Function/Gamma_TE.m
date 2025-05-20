

function S=Gamma_TE(n,w,mu,k,pho,Z)

    Zp=(1i*w*mu/k).*(spherbessH(n,1,k*pho)./(derRspherbessH(n,1,k*pho)));   % Z progressive
    Zr=(1i*w*mu/k).*(spherbessH(n,2,k*pho)./(derRspherbessH(n,2,k*pho)));   % Z regressive
    R=(Zp./Zr);
    S=(Z-Zp)./(Zp-(Z.*R));

end 

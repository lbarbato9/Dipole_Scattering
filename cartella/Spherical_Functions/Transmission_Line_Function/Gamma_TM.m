

function S=Gamma_TM(n,w,mu,k,pho,Z)

    Zp=(1i*w*mu/k).*(derRspherbessH(n,1,k*pho)./(spherbessH(n,1,k*pho)));   % Z progressive
    Zr=(1i*w*mu/k).*(derRspherbessH(n,2,k*pho)./(spherbessH(n,2,k*pho)));   % Z regressive
    R=(Zp./Zr);
    S=(Z-Zp)./(Zp-(Z.*R));

end 

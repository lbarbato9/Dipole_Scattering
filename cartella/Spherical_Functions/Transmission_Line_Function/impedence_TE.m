function [Z,A]=impedence_TE(n,w,mu,k,a_before,a_current,Z_before_boundary)

        Zj=(1i*w*mu/k)*(spherbessJ(n,k*a_before)/(derRspherbessJ(n,k*a_before))); 
        t=(spherbessY(n,k*a_before))/(spherbessJ(n,k*a_before));             
        tp=(derRspherbessY(n,k*a_before))/(derRspherbessJ(n,k*a_before));   
        A=((1i*t*Zj)-(1i*tp*Z_before_boundary))./(Z_before_boundary-Zj);
        A=1e-16*(1+1i)+A;

        Zj=(1i*w*mu/k)*(spherbessJ(n,k*a_current)/(derRspherbessJ(n,k*a_current))); 
        t=(spherbessY(n,k*a_current))/(spherbessJ(n,k*a_current));             
        tp=(derRspherbessY(n,k*a_current))/(derRspherbessJ(n,k*a_current));
        Z=(Zj*((A+(1i*t))/(A+(1i*tp))));

end
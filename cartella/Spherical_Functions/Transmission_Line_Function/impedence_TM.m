function [Z,A]=impedence_TM(mode,w,mu,k,a_before,a_current,Z_before_boundary)

% Evaluation of A, Transport coefficient at current boundary (beginning of a new medium) 
Zj=(1i*w*mu/k)*(derRspherbessJ(mode,k*a_before)/(spherbessJ(mode,k*a_before))); 
t=(spherbessY(mode,k*a_before))/(spherbessJ(mode,k*a_before));             
tp=(derRspherbessY(mode,k*a_before))/(derRspherbessJ(mode,k*a_before));   
A=((1i*tp*Zj)-(1i*t*Z_before_boundary))./(Z_before_boundary-Zj);
A=1e-16*(1+1i)+A;


% Evaluation of the impedence at the next boundary of the current medium (end of
% the medium)
Zj=(1i*w*mu/k)*(derRspherbessJ(mode,k*a_current)/(spherbessJ(mode,k*a_current))); 
t=(spherbessY(mode,k*a_current))/(spherbessJ(mode,k*a_current));             
tp=(derRspherbessY(mode,k*a_current))/(derRspherbessJ(mode,k*a_current));
Z=(Zj*((A+(1i*tp))/(A+(1i*t))));
end



%% SPHERICAL SCATTERING PROBLEM MULTILAYER FIELD COEFFICIENT EXPANSION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir     = './';
addpath([home_dir 'Spherical_Functions/Bessel/'])
addpath([home_dir 'Spherical_Functions/Spherical Armonics/'])
addpath([home_dir 'Spherical_Functions/Trasmission line functions/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocated space
Z_TE=zeros(mode_number,layer_number-1);
S_in_TE=zeros(mode_number,layer_number-1);
S_out_TE=zeros(mode_number,layer_number-1);
A_TE=zeros(mode_number,layer_number);
Ep_TE=zeros(mode_number,layer_number);
Er_TE=zeros(mode_number,layer_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=1:mode_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ep_TE(r,end)=Ep1_TE(r);

% Impedences at each boundary evaluation (from inside to outside medium):
% q=1
Z_TE(r,1)=(1i*w*mu(1)/k(1)).*(spherbessJ(r,k(1)*a(1))./(derRspherbessJ(r,k(1)*a(1))));

for q=2:layer_number-1
    [Z_TE(r,q),A_TE(r,q)]=impedence_TE(r,w,mu(q),k(q),a(q-1),a(q),Z_TE(r,q-1));   
    
end

% Reflection coefficients at each boundary in both upper and lower medium 
for q=1:layer_number-1 
S_in_TE(r,q)=Gamma_TE(r,w,mu(q),k(q),a(q),Z_TE(r,q));
S_out_TE(r,q)=Gamma_TE(r,w,mu(q+1),k(q+1),a(q),Z_TE(r,q));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D1=((spherbessH(r,1,k(end)*a(end-1)))./(spherbessH(r,2,k(end)*a(end-1))));
Er_TE(r,end)=Ep_TE(r,end).*(S_out_TE(r,end).*D1);

for q=layer_number-1:-1:2
G=(1+S_out_TE(r,q))/(1+S_in_TE(r,q));
D2=(spherbessH(r,1,k(q+1)*a(q))./(spherbessH(r,1,k(q)*a(q))));
Ep_TE(r,q)=Ep_TE(r,q+1).*D2.*G;
Er_TE(r,q)=((A_TE(r,q)-1)/(A_TE(r,q)+1)).*Ep_TE(r,q);

end

D3=spherbessH(r,1,k(2)*a(1))./(spherbessH(r,1,k(1)*a(1)));
Ep_TE(r,1)=Ep_TE(r,2).*D3.*((1+S_out_TE(r,1))/(1+S_in_TE(r,1)));
Er_TE(r,1)=Ep_TE(r,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end %multimodal for 
clear q mask mask1 mask11 H G E D1 D2 D3 
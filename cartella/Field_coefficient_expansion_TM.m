%% SPHERICAL SCATTMRING PROBLEM MULTILAYER FIELD COEFFICIENT EXPANSION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir     = './';
addpath([home_dir 'Spherical_Functions/Bessel/'])
addpath([home_dir 'Spherical_Functions/Spherical Armonics/'])
addpath([home_dir 'Spherical_Functions/Trasmission line functions/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocaTMd space
Z_TM=zeros(mode_number,layer_number-1);
S_in_TM=zeros(mode_number,layer_number-1);
S_out_TM=zeros(mode_number,layer_number-1);
A_TM=zeros(mode_number,layer_number);
Ep_TM=zeros(mode_number,layer_number);
Er_TM=zeros(mode_number,layer_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=1:mode_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Ep_TM(r,end)=Ep1_TM(1,r); % ultimo coeff corrisponde al primo mezzo
    
% Impedences at each boundary evaluation (from inside to outside medium):
Z_TM(r,1)=(1i*w*mu(1)/k(1)).*(derRspherbessJ(r,k(1)*a(1)))./(spherbessJ(r,k(1)*a(1)));

for q=2:layer_number-1
    [Z_TM(r,q),A_TM(r,q)]=impedence_TM(r,w,mu(q),k(q),a(q-1),a(q),Z_TM(r,q-1));   
end

% Reflection coefficients at each boundary in both upper and lower medium 
for q=1:layer_number-1 
S_in_TM(r,q)=Gamma_TM(r,w,mu(q),k(q),a(q),Z_TM(r,q));
S_out_TM(r,q)=Gamma_TM(r,w,mu(q+1),k(q+1),a(q),Z_TM(r,q));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D1=((derRspherbessH(r,1,k(end)*a(end-1)))./(derRspherbessH(r,2,k(end)*a(end-1))));
Er_TM(r,end)=Ep_TM(r,end).*(S_out_TM(r,end).*D1);

for q=layer_number-1:-1:2
G=(1+S_out_TM(r,q))/(1+S_in_TM(r,q));
D2=(derRspherbessH(r,1,k(q+1)*a(q))./(derRspherbessH(r,1,k(q)*a(q))));
Ep_TM(r,q)=Ep_TM(r,q+1).*D2.*G;
Er_TM(r,q)=((A_TM(r,q)-1)/(A_TM(r,q)+1)).*Ep_TM(r,q);

end

D3=derRspherbessH(r,1,k(2)*a(1))./(derRspherbessH(r,1,k(1)*a(1)));
Ep_TM(r,1)=Ep_TM(r,2).*D3.*((1+S_out_TM(r,1))/(1+S_in_TM(r,1)));
Er_TM(r,1)=Ep_TM(r,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear q mask mask1 mask11 H G E D1 D2 D3 


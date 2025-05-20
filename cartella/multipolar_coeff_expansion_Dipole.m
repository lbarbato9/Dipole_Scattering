%%                 DIPOLE EXPANSION COEFFICIENTS SCRIPT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir = './';
addpath([home_dir 'Spherical_Functions/Bessel/'])
addpath([home_dir 'Spherical_Functions/Spherical Armonics/'])
addpath([home_dir 'Spherical_Functions/Trasmission line functions/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocated space for variable
an = zeros(1,mode_number);
bn = zeros(1,mode_number);
pv=zeros(1,mode_number);
qv=zeros(1,mode_number);
Dv=zeros(1,mode_number);
uv=zeros(1,mode_number);
vv=zeros(1,mode_number);
piG=zeros(1,mode_number);
Pnn=zeros(1,mode_number);
tautau=zeros(1,mode_number);
M3_theta=zeros(1,mode_number);
M3_phi=zeros(1,mode_number);
N3_theta=zeros(1,mode_number);
N3_phi=zeros(1,mode_number);
N3_pho=zeros(1,mode_number);

M3_theta_odd=zeros(1,mode_number);
N3_theta_even=zeros(1,mode_number);

Ep1_TE = zeros(1,mode_number);
Ep1_TM = zeros(1,mode_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scattering Coefficients on the sphere Boundary

for n = 1 : mode_number

num1 = ric_besselj_derivative(n,k(end)*a(end-1),1) * spherbessJ(n,k(end-1)*a(end-1)) - spherbessJ(n,k(end)*a(end-1)) *  ric_besselj_derivative(n,k(end-1)*a(end-1),1);
den1 = ric_besselj_derivative(n,k(end-1)*a(end-1),1) * spherbessH(n,1,k(end)*a(end-1)) - spherbessJ(n,k(end-1)*a(end-1)) *  ric_besselh_derivative(n,1,k(end)*a(end-1),1);

an(n) = num1/den1;

num2 = ric_besselj_derivative(n,k(end-1)*a(end-1),1) * spherbessJ(n,k(end)*a(end-1)) - eps(end-1) * ric_besselj_derivative(n,k(end)*a(end-1),1) * spherbessJ(n,k(end-1)*a(end-1));
den2 = eps(end-1) * ric_besselh_derivative(n,1,k(end)*a(end-1),1) * spherbessJ(n,k(end-1)*a(end-1)) - ric_besselj_derivative(n,k(end-1)*a(end-1),1) * spherbessH(n,1,k(end)*a(end-1));

bn(n) = num2/den2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radial Oriented Dipole Coefficient Expansion

if pol==0

m=0;
eps_m = 1;

Pnn(1,n)=Pnm(n,m,theta_d);
tautau(1,n)=Tau(n,m,theta_d);
piG(1,n)=Pigreco(n,m,theta_d);
Dv(1,n) =eps_m *((2*n+1)*factorial(n-m))/(4*n*(n+1)*factorial(n+m));
N3_pho(1,n)=   (spherbessH(n,1,k(end)*r_d)./(k(end).*r_d)).*(exp(1i*m*phi_d)).*(n.*(n+1)).*Pnn(1,n);
qv(1,n) = (1i*k(end)^3/(eps_r(end)*pi))*N3_pho(1,n);

uv(1,n) = an(1,n) * pv(m+1,n) ;
vv(1,n) = bn(1,n) * qv(m+1,n) ;

Ep1_TE(1,n) =  (Dv(1,n)*(pv(1,n)/ 2) + uv(1,n));
Ep1_TM(1,n) =  (Dv(1,n)*(qv(1,n)/ 2) + vv(1,n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radial Oriented Dipole Coefficient Expansion

elseif pol==1
m=1;
eps_m = 2;
Dv(1,n) =eps_m *((2*n+1)*factorial(n-m))/(4*n*(n+1)*factorial(n+m));
Pnn(1,n)=Pnm(n,m,theta_d);
tautau(1,n)=Tau(n,m,theta_d);
piG(1,n)=Pigreco(n,m,theta_d);

% M3_theta_odd(1,n)= (spherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*((piG(1,n)))*(1/2);
% N3_theta_even(1,n)=-((spherbessH(n,1,k(end))/(k(end)*r_d)) +  derRspherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*(tautau(1,n))*(1/2);   
% N3_theta_even(1,n)=-(derRspherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*(tautau(1,n))*(1/2);   



M3_theta_odd(1,n)= (spherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*((piG(1,n)))*(1/2);
N3_theta_even(1,n)=-((spherbessH(n,1,k(end))/(k(end)*r_d)) +  derRspherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*(tautau(1,n))*(1/2);   
N3_theta_even(1,n)=-(derRspherbessH(n,1,k(end)*r_d)).*(cos(m*phi_d)).*(tautau(1,n))*(1/2); 




pv(1,n) = (1i*k(end)^3/(eps_r(end)*pi))*p0*M3_theta_odd(1,n);
qv(1,n) = (1i*k(end)^3/(eps_r(end)*pi))*p0*N3_theta_even(1,n);

uv(1,n) = an(1,n) * pv(1,n) ;
vv(1,n) = bn(1,n) * qv(1,n) ;

Ep1_TE(1,n) =  (Dv(1,n)*(pv(1,n)/ 2) + uv(1,n));
Ep1_TM(1,n) =  (Dv(1,n)*(qv(1,n)/ 2) + vv(1,n));

end


end %sum over n

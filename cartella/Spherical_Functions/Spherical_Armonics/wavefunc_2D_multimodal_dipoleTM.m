%                VECTOR HARMONICS SPHERICAL CASE 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function describe vector harmonics in Spherical case, evaluating
% fields component in a tidimensional space for a Loop Coil antenna.
%
%  INPUT:         
%                n        in analysis mode 
%                zita     Medium impedence (sqrt(mu/eps)) [Ohm]
%                k        vector wave components [m^-1]
%                pho      radial variable [m]
%                Ep,Er    Progressive and regressive wave coefficient for each mode, they are vectors 
%
%
% OUTPUT:       Vectors filled with fields component:
%
%               H_phi, H_pho, H_theta
%               E_phi, E_pho, E_theta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:      Vincenzo Miranda
% Last Update: 20/05/204
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H_phi,H_pho,H_theta,E_phi,E_pho,E_theta] = wavefunc_2D_multimodal_dipoleTM(pol,n_tot,n_points,w,mu,k,pho,theta,phi,Ep,Er)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_phi=zeros(n_points,n_points);                                     %%%
H_pho=zeros(n_points,n_points);                                     %%%
H_theta=zeros(n_points,n_points);                                   %%%
E_phi=zeros(n_points,n_points);                                     %%%
E_pho=zeros(n_points,n_points);                                     %%%
E_theta=zeros(n_points,n_points);                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pol==0
m=0; % Radial Oriented Dipole
elseif pol==1
m=1; % Tangentially oriented Dipole 
end

for n=1:n_tot

Pnn=Pnm(n,m,theta);
tautau=Tau(n,m,theta);
piG=Pigreco(n,m,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pol==0
M3_phi=   (-spherbessH(n,1,k*pho)).*(tautau).*(exp(1i*m*phi));  
M4_phi=   (-spherbessH(n,2,k*pho)).*(tautau).*(exp(1i*m*phi)); 

M3_theta= (spherbessH(n,1,k*pho)).*((piG)).*(exp(1i*m*phi)); 
M4_theta= (spherbessH(n,2,k*pho)).*((piG)).*(exp(1i*m*phi)); 

N3_pho=   (spherbessH(n,1,k*pho)./(k.*pho)).*(n.*(n+1)).*Pnn.*(exp(1i*m*phi)); 
N3_theta= (derRspherbessH(n,1,k*pho)).*(tautau).*(exp(1i*m*phi));   
N3_phi=   (derRspherbessH(n,1,k*pho)).*(piG).*(exp(1i*m*phi));                

N4_pho=   (spherbessH(n,2,k*pho)./(k.*pho)).*(n.*(n+1)).*Pnn.*(exp(1i*m*phi)); 
N4_theta= (derRspherbessH(n,2,k*pho)).*(tautau).*(exp(1i*m*phi)); 
N4_phi=   (derRspherbessH(n,2,k*pho)).*(piG).*(exp(1i*m*phi)); 

elseif pol==1

M3_phi=   (-spherbessH(n,1,k*pho)).*(tautau).*(sin(m*phi));  
M4_phi=   (-spherbessH(n,2,k*pho)).*(tautau).*(sin(m*phi)); 

M3_theta= (spherbessH(n,1,k*pho)).*((piG)).*cos(m*phi);
M4_theta= (spherbessH(n,2,k*pho)).*((piG)).*cos(m*phi); 

N3_pho=   (spherbessH(n,1,k*pho)./(k.*pho)).*(n.*(n+1)).*Pnn.*cos(m*phi);
N3_theta= (derRspherbessH(n,1,k*pho)).*(tautau).*cos(m*phi);   
N3_phi=   (-derRspherbessH(n,1,k*pho)).*(piG).*(sin(m*phi));                

N4_pho=   (spherbessH(n,2,k*pho)./(k.*pho)).*(n.*(n+1)).*Pnn.*cos(m*phi); 
N4_theta= (derRspherbessH(n,2,k*pho)).*(tautau).*cos(m*phi);  
N4_phi=   (-derRspherbessH(n,2,k*pho)).*(piG).*(sin(m*phi));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TM POLARIZATION 

Htheta_field= (Ep(n)*M3_theta + Er(n)*M4_theta); 
H_theta= H_theta + Htheta_field;

Hphi_field=(Ep(n)*M3_phi + Er(n)*M4_phi); 
H_phi= H_phi + Hphi_field;

Ephi_field= Ep(n)*N3_phi + Er(n)*N4_phi;
E_phi= E_phi+ Ephi_field;

Epho_field= Ep(n)*N3_pho + Er(n)*N4_pho;
E_pho= E_pho +Epho_field;

Etheta_field= Ep(n)*N3_theta + Er(n)*N4_theta;
E_theta= E_theta+ Etheta_field;

end

H_theta=(k/(1i*w*mu))*H_theta;
H_phi=(k/(1i*w*mu))*H_phi;

end % end function
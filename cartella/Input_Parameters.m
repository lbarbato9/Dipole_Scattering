%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          INPUT PARAMETER SPHERICAL SCATTERING LOOP COIL - MULTILAYER
%
% This file create the input parameters for the solution of the scattering
% problem from a layered sphere, with a loop coil as input source. 
% 
% Input General parameters--> Frequency, mode numbers, layer numbers
%
% Input Electromagnetic parameter --> K and Zita and geometrical parameters 
% for N-medium. This script allow to use metamaterials also.
%
% Electromagnetic Parameters are summed up as follow:
%
% eps=[eps n,eps n-1,eps n-2 ,...eps1] 
%
% from the internal to external medium, be careful!
%
% Input loop coil --> Radius, distance from origin and Current are required
% the source is placed algong z-axis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Author: Vincenzo Miranda
%Last Update: 15/09/24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_0=8.8540*1e-12; %vaccum permittivity value
mu_0=4*pi*1e-7;   %vacuum permeability value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Input parameters

    prompt={'Frequency [Hz]','Modes Number','Layers Number','Sampling period [m]'};
    name='Incident Wave Parameters';
    defaultanswer={'297.2*1e6','5','2','0.008'};
    numlines=1;
    options.Windowstyle='normal';
    options.resize=10;
    data0=inputdlg(prompt,name,numlines,defaultanswer,options); 

f=str2num(data0{1,1});         %[Hz]
w=2*pi*f;                      %[Rad/s]
mode_number=str2num(data0{2,1});
layer_number=str2num(data0{3,1});
dl=str2num(data0{4,1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coil Input Parameters
% 
prompt={'Dipole Moment p','Dipole distance','Dipole Orientation 0) Radial 1) Tangential'};
name='Dipole Parameters';
defaultanswer={'1','0.12','0'};
numlines=1;
options.Resize=2;
options.Windowstyle='normal';
data_dipole=inputdlg(prompt,name,numlines,defaultanswer,options);

p0=str2num(data_dipole{1,1});
r_d=str2num(data_dipole{2,1});
pol=str2num(data_dipole{3,1});
theta_d=0;
phi_d=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electromagnetic Mediums Parameters 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocated space for variable storadge:
eps_r=zeros(layer_number,1);            %%%
eps=zeros(layer_number,1);              %%%
mu_r=zeros(layer_number,1);             %%%
mu=zeros(layer_number,1);               %%%
sigma=zeros(layer_number,1);            %%%
d=zeros(layer_number,1);                %%%
k=zeros(layer_number,1);                %%%
zita=zeros(layer_number,1);             %%%
a=zeros(layer_number,1);                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for q=1:layer_number 

    prompt1={"Relative Permittivity medium " + q ,"Conductivity [S/m] medium " + q,"Relative Permeability " + q,"Layer thickness medium [m] "  + q };
    name='Material"s Parameters';
    defaultanswer={'50','0','1','0.09'};
    numlines=1;
    options.Resize=2;
    options.Windowstyle='normal';
    data1=inputdlg(prompt1,name,numlines,defaultanswer,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting input data from prompt table:

    eps_r(q,1)=str2num(data1{1,1});
    eps(q,1)=eps_r(q,1)*eps_0;
    sigma(q,1)=str2num(data1{2,1});
    if sigma(q,1)==0
        sigma(q,1)=1e-9;
    end
    mu_r(q,1)=str2num(data1{3,1});
    mu(q,1)=mu_r(q,1)*mu_0;
    d(q,1)=str2num(data1{4,1});

    if q==1 
    a(q,1)=d(q,1);
    else
    a(q,1)=a(q-1,1)+d(q,1); 
    end

%%%%%%%%%%% Wave vectors Components (Metamaterials Added) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluating all roots of wave-number for metamaterials correct parameter input   
    polynomial_k=w^2*eps(q,1)*mu(q,1)-1i*w*sigma(q,1)*mu(q,1);
    polynomial_kroots=roots([1 0 -polynomial_k]);
    polynomial_zita=mu(q,1)/eps(q,1);
    polynomial_zitaroots=roots([1 0 -polynomial_zita]);

% DPS (Re{eps}>0, Re{mu}>0) --->       [zita>0]                      [k>0]
% 
% DNG (Re{eps}<0, Re{mu}<0) --->       [zita>0]                      [k>0]
% 
% ENG (Re{eps}<0, Re{mu}>0) --->       [zita<0 ,Complex pure]        [k>0, Complex pure]
% 
% MNG (Re{eps}>0, Re{mu}<0) --->       [zita>0 ,Complex pure]        [k>0, Complex pure]

% DPS 
if real(eps(q,1))>0 && real(mu(q,1))>0
   k(q,1)=polynomial_kroots(1,1);
   zita(q,1)=polynomial_zitaroots(1,1);
%DNG
elseif real(eps(q,1))<0 && real(mu(q,1))<0 
   k(q,1)=polynomial_kroots(2,1);
   zita(q,1)=polynomial_zitaroots(1,1);
%ENG
elseif real(eps(q,1))<0 && real(mu(q,1))>0 
   k(q,1)=polynomial_kroots(1,1);
   zita(q,1)=polynomial_zitaroots(2,1);
%MNG
elseif real(eps(q,1))>0 && real(mu(q,1))<0 
   k(q,1)=polynomial_kroots(1,1);
   zita(q,1)=polynomial_zitaroots(1,1);

end %if material parameter

end % main for 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical 3D Geometry Definition
FOV = 2 * a(end);
N_points = round(FOV / dl);

% Defining X,Y and Z-axis
ll = linspace(-a(end), a(end), N_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XY plane
[XX1, YY1] = meshgrid(ll, ll);

% Calcolo delle coordinate sferiche
pho_xy = sqrt(YY1.^2 + XX1.^2);
phi_xy = atan2(YY1, XX1);
phi_xy = wrapTo2Pi(phi_xy); % Porta phi tra [0, 2*pi]
theta_xy = pi/2 * ones(size(XX1)); % Equatore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YZ PLANE
[YY2,ZZ2] = meshgrid(ll, ll);

% Calcolo delle coordinate sferiche
pho_yz = sqrt(YY2.^2 + ZZ2.^2);
theta_yz = acos(YY2./pho_yz); % Correzione nella formula per latitudine
phi_yz = pi/2 * ones(size(YY2)); % Longitudine fissa a pi/2
phi_yz(ZZ2 < 0) = 3 * pi / 2; % Aggiusta per Y negativo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZX plane
[XX3, ZZ3] = meshgrid(ll, ll);

% Calcolo delle coordinate sferiche
YY3 = zeros(size(XX3)); % Piano sagittale a y = 0
pho_xz = sqrt(XX3.^2 + ZZ3.^2);
theta_xz = acos(XX3./pho_xz); % Correzione della formula
phi_xz = zeros(size(XX3)); % Longitudine inizialmente 0
phi_xz(ZZ3 < 0) = pi; % Aggiusta per X negativo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ro=zeros(3,length(XX3),length(XX3));
theta=zeros(3,length(XX3),length(XX3));
phi=zeros(3,length(XX3),length(XX3));

ro(1,:,:)=pho_xz;
ro(2,:,:)=pho_xy;
ro(3,:,:)=pho_yz;

pho(1,:,:)=phi_xz;
phi(2,:,:)=phi_xy;
phi(3,:,:)=phi_yz;

theta(1,:,:)=theta_xz;
theta(2,:,:)=theta_xy;
theta(3,:,:)=theta_yz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control plot for spherical coordinate system
% 
% Primo indice nella matrice è il piano X-Z 
% Secondo indice nella matrice è il piano X-Y
% terzo indice nella matrice è il piano Y-Z 

% figure
% subplot(3,2,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],theta_xz)
% colorbar 
% title('\theta_X_Z')
% ylabel('X')
% xlabel('Z')
% pbaspect([1,1,1])
% 
% subplot(3,2,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],phi_xz)
% colorbar
% title('\phi_X_Z')
% ylabel('X')
% xlabel('Z')
% pbaspect([1,1,1])
% 
% subplot(3,2,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],theta_yz)
% colorbar 
% title('\theta_X_Z')
% ylabel('Y')
% xlabel('Z')
% pbaspect([1,1,1])
% 
% subplot(3,2,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],phi_yz)
% colorbar 
% title('\phi_Y_Z')
% ylabel('Y')
% xlabel('Z')
% pbaspect([1,1,1])
% 
% subplot(3,2,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],theta_xy)
% colorbar 
% title('\theta_X_Y')
% ylabel('Y')
% xlabel('X')
% pbaspect([1,1,1])
% 
% subplot(3,2,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],phi_xy)
% colorbar 
% title('\phi_X_Y')
% ylabel('Y')
% xlabel('X')
% pbaspect([1,1,1])

sgtitle('Control plot for spherical coordinate system')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear data0 defaultanswer name numlines options prompt0 data1 prompt1 data_loopcoil prompt dims
clear polynomial_zitaroots polynomial_k polynomial_kroots polynomial_zita aus1 p t
clear data_input_table1 data_input_table0 input_table1 input_table0 home_dir index_field_plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EM-FIELD EVALUATION FROM EXPANSION COEFFICIENTS:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir     = './';
addpath([home_dir 'Spherical_Functions/Bessel/'])
addpath([home_dir 'Spherical_Functions/Spherical Armonics/'])
addpath([home_dir 'Spherical_Functions/Trasmission line functions/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocaTMd space TM

% Space for the multimodal field in all medium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_pho_TM_acc=zeros(3,layer_number,N_points,N_points);
E_theta_TM_acc=zeros(3,layer_number,N_points,N_points);
E_phi_TM_acc=zeros(3,layer_number,N_points,N_points);
H_theta_TM_acc=zeros(3,layer_number,N_points,N_points);
H_phi_TM_acc=zeros(3,layer_number,N_points,N_points);
H_pho_TM_acc=zeros(3,layer_number,N_points,N_points);

acc1=zeros(3,N_points,N_points);
acc2=zeros(3,N_points,N_points);
acc3=zeros(3,N_points,N_points);
acc4=zeros(3,N_points,N_points);
acc5=zeros(3,N_points,N_points);
acc6=zeros(3,N_points,N_points);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TM 
% H_phi, H_pho=0, H_theta
% E_phi, E_pho, E_theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0, 'Progress...');  % Crea una barra di avanzamento
figure

for plane=1:3
for q=1:layer_number                  %medium for 

[acc1(plane,:,:),acc2(plane,:,:),acc3(plane,:,:),acc4(plane,:,:),acc5(plane,:,:),acc6(plane,:,:)] = wavefunc_2D_multimodal_dipoleTM(pol,mode_number,N_points,w,mu(q),k(q),squeeze(ro(plane,:,:)),squeeze(theta(plane,:,:)),squeeze(phi(plane,:,:)),Ep_TM(:,q),Er_TM(:,q));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q==1 %internal medium mask
ro_mask= squeeze(ro(plane,:,:))<=a(q);

elseif q==layer_number %external medium mask
mask11= squeeze(ro(plane,:,:))>a(q-1);
mask1=squeeze(ro(plane,:,:))<=r_d;
ro_mask=mask1.*mask11;

else % intermediate medium mask
mask11=squeeze(ro(plane,:,:))>=a(q-1,1);
mask1=squeeze(ro(plane,:,:))<a(q,1);
ro_mask=mask1.*mask11;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(layer_number,1,q)
imagesc([-FOV/2 FOV/2], [-FOV/2 FOV/2], ro_mask);
title(sprintf('Layer: %d, eps = %.2f, sigma = %.2f S/m, a= %.2f m', q, eps_r(q), sigma(q), a(q)), 'FontSize', 12); % Inserisce i valori delle variabili nel titolo
pbaspect([1 1 1]);
xlabel('[m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('[m]', 'FontSize', 10, 'FontWeight', 'bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_phi_TM_acc(plane,q,:,:)=squeeze(acc1(plane,:,:)).*ro_mask;
H_pho_TM_acc(plane,q,:,:)=squeeze(acc2(plane,:,:)).*ro_mask;
H_theta_TM_acc(plane,q,:,:)=squeeze(acc3(plane,:,:)).*ro_mask;

E_phi_TM_acc(plane,q,:,:)=squeeze(acc4(plane,:,:)).*ro_mask;
E_pho_TM_acc(plane,q,:,:)=squeeze(acc5(plane,:,:)).*ro_mask;
E_theta_TM_acc(plane,q,:,:)=squeeze(acc6(plane,:,:)).*ro_mask;

 progress = q/layer_number;
 waitbar(progress, h, sprintf('Simulation in Progress TM case: %d%%', round(progress*100)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Total fields TM sum of the selected modes in the medium between the sphere
% and the antenna
E_theta_TM(:,:,:)=sum(E_theta_TM_acc,2);
E_phi_TM(:,:,:)=sum(E_phi_TM_acc,2);
E_pho_TM(:,:,:)=sum(E_pho_TM_acc,2);
H_theta_TM(:,:,:)=sum(H_theta_TM_acc,2);
H_phi_TM(:,:,:)=sum(H_phi_TM_acc,2);
H_pho_TM(:,:,:)=sum(H_pho_TM_acc,2);

coil_mask= (ro(plane,:,:))<=r_d;

E_theta_TM(:,:,:)=E_theta_TM.*coil_mask;
E_phi_TM(:,:,:)=E_phi_TM.*coil_mask;
E_pho_TM(:,:,:)=E_pho_TM.*coil_mask;
H_theta_TM(:,:,:)=H_theta_TM.*coil_mask;
H_phi_TM(:,:,:)=H_phi_TM.*coil_mask;
H_pho_TM(:,:,:)=H_pho_TM.*coil_mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluation of absolut value of total fields multimodal and for each layer:
H_TM= abs(H_theta_TM).^2 + abs(H_pho_TM).^2 +abs(H_phi_TM).^2;
H_TM= sqrt(H_TM);

E_TM= abs(E_theta_TM).^2  +abs(E_phi_TM).^2 +abs(E_pho_TM).^2;
E_TM= sqrt(E_TM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EM-FIELD EVALUATION FROM EXPANSION COEFFICIENTS:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home_dir     = './';
addpath([home_dir 'Spherical_Functions/Bessel/'])
addpath([home_dir 'Spherical_Functions/Spherical Armonics/'])
addpath([home_dir 'Spherical_Functions/Trasmission line functions/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-allocated space TE

% Space for the multimodal field in all medium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_pho_TE_acc=zeros(3,layer_number,N_points,N_points);
E_theta_TE_acc=zeros(3,layer_number,N_points,N_points);
E_phi_TE_acc=zeros(3,layer_number,N_points,N_points);
H_theta_TE_acc=zeros(3,layer_number,N_points,N_points);
H_phi_TE_acc=zeros(3,layer_number,N_points,N_points);
H_pho_TE_acc=zeros(3,layer_number,N_points,N_points);

acc1=zeros(3,N_points,N_points);
acc2=zeros(3,N_points,N_points);
acc3=zeros(3,N_points,N_points);
acc4=zeros(3,N_points,N_points);
acc5=zeros(3,N_points,N_points);
acc6=zeros(3,N_points,N_points);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TE 
% H_phi, H_pho, H_theta
% E_phi, E_pho=0, E_theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for plane=1:3

for q=1:layer_number                  %medium for 

[acc1(plane,:,:),acc2(plane,:,:),acc3(plane,:,:),acc4(plane,:,:),acc5(plane,:,:),acc6(plane,:,:)] = wavefunc_2D_multimodal_dipoleTE(pol,mode_number,N_points,w,mu(q),k(q),squeeze(ro(plane,:,:)),squeeze(theta(plane,:,:)),squeeze(phi(plane,:,:)),Ep_TE(:,q),Er_TE(:,q));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q==1 %internal medium mask
ro_mask= squeeze(ro(plane,:,:))<=a(q);

elseif q==layer_number %external medium mask
ro_mask= squeeze(ro(plane,:,:))>a(q-1);

else %intermediate medium mask
mask11=squeeze(ro(plane,:,:))>=a(q-1,1);
mask1=squeeze(ro(plane,:,:))<a(q,1);
ro_mask=mask1.*mask11;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_phi_TE_acc(plane,q,:,:)=squeeze(acc1(plane,:,:)).*ro_mask;
H_pho_TE_acc(plane,q,:,:)=squeeze(acc2(plane,:,:)).*ro_mask;
H_theta_TE_acc(plane,q,:,:)=squeeze(acc3(plane,:,:)).*ro_mask;

E_phi_TE_acc(plane,q,:,:)=squeeze(acc4(plane,:,:)).*ro_mask;
E_pho_TE_acc(plane,q,:,:)=squeeze(acc5(plane,:,:)).*ro_mask;
E_theta_TE_acc(plane,q,:,:)=squeeze(acc6(plane,:,:)).*ro_mask;

 % 
 % progress = q / layer_number;
 % waitbar(progress, h, sprintf('Simulation in Progress TE case: %d%%', round(progress*100)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total fields TE sum of the selected modes in the medium q

E_theta_TE(:,:,:)=sum(E_theta_TE_acc,2);
E_phi_TE(:,:,:)=sum(E_phi_TE_acc,2);
E_pho_TE(:,:,:)=sum(E_pho_TE_acc,2);
H_theta_TE(:,:,:)=sum(H_theta_TE_acc,2);
H_phi_TE(:,:,:)=sum(H_phi_TE_acc,2);
H_pho_TE(:,:,:)=sum(H_pho_TE_acc,2);

coil_mask= (ro(plane,:,:))<=r_d;

E_theta_TE(:,:,:)=E_theta_TE.*coil_mask;
E_phi_TE(:,:,:)=E_phi_TE.*coil_mask;
E_pho_TE(:,:,:)=E_pho_TE.*coil_mask;
H_theta_TE(:,:,:)=H_theta_TE.*coil_mask;
H_phi_TE(:,:,:)=H_phi_TE.*coil_mask;
H_pho_TE(:,:,:)=H_pho_TE.*coil_mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluation of absolut value of total fields multimodal and for each layer:

H_TE= abs(H_theta_TE).^2 + abs(H_pho_TE).^2 +abs(H_phi_TE).^2;
H_TE= sqrt(H_TE);

E_TE= abs(E_theta_TE).^2  +abs(E_phi_TE).^2 +abs(E_pho_TE).^2;
E_TE= sqrt(E_TE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(h);
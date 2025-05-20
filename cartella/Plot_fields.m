%% PLOT SECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th=linspace(0,2*pi,1000);
x_a=zeros(layer_number,length(th));
y_a=zeros(layer_number,length(th));

for q=1:layer_number
x_a(q,:)= a(q)*cos(th);
y_a(q,:)= a(q)*sin(th);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TE and TM FIELDS (Absolute Value)                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Primo indice nella matrice è il piano X-Z 
% Secondo indice nella matrice è il piano X-Y
% terzo indice nella matrice è il piano Y-Z 

pE_TE_xz=squeeze(E_TE(1,:,:));
pH_TE_xz=squeeze(H_TE(1,:,:))*mu_0*1e6;
pE_TE_yz=squeeze(E_TE(3,:,:));
pH_TE_yz=squeeze(H_TE(3,:,:))*mu_0*1e6;
pE_TE_xy=squeeze(E_TE(2,:,:));
pH_TE_xy=squeeze(H_TE(2,:,:))*mu_0*1e6;

pE_TM_xz=squeeze(E_TM(1,:,:));
pH_TM_xz=squeeze(H_TM(1,:,:))*mu_0*1e6;
pE_TM_yz=squeeze(E_TM(3,:,:));
pH_TM_yz=squeeze(H_TM(3,:,:))*mu_0*1e6;
pE_TM_xy=squeeze(E_TM(2,:,:));
pH_TM_xy=squeeze(H_TM(2,:,:))*mu_0*1e6;

ro_mask_int= squeeze(ro(1,:,:))<=a(1);
ro_mask_dipole= squeeze(ro(1,:,:))<=r_d;
%ro_mask_dipole= ro(2,:,:)<=a(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TM + TE FIELDS                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Magnetic Field 
H_theta= H_theta_TE + H_theta_TM;
H_phi= H_phi_TE + H_phi_TM;
H_pho= H_pho_TE;

H= abs(H_theta).^2 + abs(H_pho).^2 + abs(H_phi).^2 ;    %|B|^2
H= sqrt(H);                                             %|B|

%Electric Field:
E_theta= E_theta_TE + E_theta_TM;
E_phi= E_phi_TE + E_phi_TM;
E_pho= E_pho_TM;
E= abs(E_theta).^2 + abs(E_pho).^2 + abs(E_phi).^2 ;    %|B|^2
E= sqrt(E);    

% Primo indice nella matrice è il piano X-Z 
% Secondo indice nella matrice è il piano X-Y
% terzo indice nella matrice è il piano Y-Z 

pE_xz=squeeze(E(1,:,:));
pH_xz=squeeze(H(1,:,:))*mu_0*1e6;
pE_yz=squeeze(E(3,:,:));
pH_yz=squeeze(H(3,:,:))*mu_0*1e6;
pE_xy=squeeze(E(2,:,:));
pH_xy=squeeze(H(2,:,:))*mu_0*1e6;
% pH_xy_rotata = imrotate(pH_xy, 90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
% subplot(2,3,1)
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pH_xz.*ro_mask_int)
colormap jet
colorbar 
pbaspect([1 1 1]) 
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Z [m]','FontSize',10,'FontWeight','bold')
ylabel('X [m]','FontSize',10,'FontWeight','bold')
title('|B_z_x| [μT]')

figure
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],rot90(pH_xy.*ro_mask_int./(0.99*(max(pH_xy.*ro_mask_int,[],'all')))))
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],rot90(pH_xy.*ro_mask_int))

colormap jet
colorbar 
pbaspect([1 1 1]) 
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Y [m]','FontSize',10,'FontWeight','bold')
ylabel('X [m]','FontSize',10,'FontWeight','bold')
title('|B_x_y| [μT]')

figure
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pH_yz.*ro_mask_int)
colormap jet
colorbar 
pbaspect([1 1 1])
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Y [m]','FontSize',10,'FontWeight','bold')
ylabel('Z [m]','FontSize',10,'FontWeight','bold')
title('|B_zy| [μT]')

figure
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_xz.*ro_mask_int)
colormap jet
colorbar 
pbaspect([1 1 1])
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Z [m]','FontSize',10,'FontWeight','bold')
ylabel('X [m]','FontSize',10,'FontWeight','bold')
title('|E_z_x| [V/m]')

figure
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_xy.*ro_mask_int)
colormap jet
colorbar
pbaspect([1 1 1])
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Y [m]','FontSize',10,'FontWeight','bold')
ylabel('X [m]','FontSize',10,'FontWeight','bold')
title('|E_x_y| [V/m]')

figure
imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_yz.*ro_mask_int)
colormap jet
colorbar
pbaspect([1 1 1])
hold on 
plot(x_a(1,:),y_a(1,:),'k--',LineWidth=1)
hold off
xlabel('Y [m]','FontSize',10,'FontWeight','bold')
ylabel('Z [m]','FontSize',10,'FontWeight','bold')
title('|E_z_y| [V/m]')

sgtitle('TM + TE Fields in the inner medium','Fontsize',12,'FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pH_xz.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_z_x| [μT]')
% 
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pH_xy.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_x_y| [μT]')
% 
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pH_yz.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_zy| [μT]')
% 
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_xz.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_z_x| [V/m]')
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_xy.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_x_y| [V/m]')
% 
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],pE_yz.*ro_mask_dipole)
% colormap jet
% colorbar
% pbaspect([1 1 1]) 
% hold on 
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_z_y| [V/m]')
% 
% sgtitle('TM + TE Fields in the whole space','Fontsize',12,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                      TE FIELDS COMPONENTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Primo indice nella matrice è il piano X-Z 
% % Secondo indice nella matrice è il piano Y-Z 
% % terzo indice nella matrice è il piano X-Y 
% 
% pE_phi_TE_yz=squeeze(E_phi_TE(3,:,:));
% pH_phi_TE_yz=squeeze(H_phi_TE(3,:,:))*mu_0*1e6;
% pE_theta_TE_yz=squeeze(E_theta_TE(3,:,:));
% pH_theta_TE_yz=squeeze(H_theta_TE(3,:,:))*mu_0*1e6;
% pE_pho_TE_yz=squeeze(E_pho_TE(3,:,:));
% pH_pho_TE_yz=squeeze(H_pho_TE(3,:,:))*mu_0*1e6;
% 
% pE_phi_TE_xz=squeeze(E_phi_TE(1,:,:));
% pH_phi_TE_xz=squeeze(H_phi_TE(1,:,:))*mu_0*1e6;
% pE_theta_TE_xz=squeeze(E_theta_TE(1,:,:));
% pH_theta_TE_xz=squeeze(H_theta_TE(1,:,:))*mu_0*1e6;
% pE_pho_TE_xz=squeeze(E_pho_TE(1,:,:));
% pH_pho_TE_xz=squeeze(H_pho_TE(1,:,:))*mu_0*1e6;
% 
% pE_phi_TE_xy=squeeze(E_phi_TE(2,:,:));
% pH_phi_TE_xy=squeeze(H_phi_TE(2,:,:))*mu_0*1e6;
% pE_theta_TE_xy=squeeze(E_theta_TE(2,:,:));
% pH_theta_TE_xy=squeeze(H_theta_TE(2,:,:))*mu_0*1e6;
% pE_pho_TE_xy=squeeze(E_pho_TE(2,:,:));
% pH_pho_TE_xy=squeeze(H_pho_TE(2,:,:))*mu_0*1e6;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TE_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TE_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TE_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TE_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TE_xy).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TE_xy).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TE Fields XY Plane','Fontsize',12,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TE_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TE_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TE_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TE_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TE_yz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TE_yz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TE Fields YZ Plane','Fontsize',12,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TE_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TE_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TE_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TE_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TE_xz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TE_xz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TE Fields ZX Plane','Fontsize',12,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                      TM FIELDS COMPONENTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Primo indice nella matrice è il piano X-Z 
% % Secondo indice nella matrice è il piano Y-Z 
% % terzo indice nella matrice è il piano X-Y 
% 
% pE_phi_TM_yz=squeeze(E_phi_TM(3,:,:));
% pH_phi_TM_yz=squeeze(H_phi_TM(3,:,:))*mu_0*1e6;
% pE_theta_TM_yz=squeeze(E_theta_TM(3,:,:));
% pH_theta_TM_yz=squeeze(H_theta_TM(3,:,:))*mu_0*1e6;
% pE_pho_TM_yz=squeeze(E_pho_TM(3,:,:));
% pH_pho_TM_yz=squeeze(H_pho_TM(3,:,:))*mu_0*1e6;
% 
% pE_phi_TM_xz=squeeze(E_phi_TM(1,:,:));
% pH_phi_TM_xz=squeeze(H_phi_TM(1,:,:))*mu_0*1e6;
% pE_theta_TM_xz=squeeze(E_theta_TM(1,:,:));
% pH_theta_TM_xz=squeeze(H_theta_TM(1,:,:))*mu_0*1e6;
% pE_pho_TM_xz=squeeze(E_pho_TM(1,:,:));
% pH_pho_TM_xz=squeeze(H_pho_TM(1,:,:))*mu_0*1e6;
% 
% pE_phi_TM_xy=squeeze(E_phi_TM(2,:,:));
% pH_phi_TM_xy=squeeze(H_phi_TM(2,:,:))*mu_0*1e6;
% pE_theta_TM_xy=squeeze(E_theta_TM(2,:,:));
% pH_theta_TM_xy=squeeze(H_theta_TM(2,:,:))*mu_0*1e6;
% pE_pho_TM_xy=squeeze(E_pho_TM(2,:,:));
% pH_pho_TM_xy=squeeze(H_pho_TM(2,:,:))*mu_0*1e6;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TM_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TM_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TM_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TM_xy).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TM_xy).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TM_xy).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TM Fields XY Plane','Fontsize',12,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TM_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TM_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TM_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TM_yz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TM_yz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TM_yz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Y [m]','FontSize',10,'FontWeight','bold')
% ylabel('Z [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TM Fields YZ Plane','Fontsize',12,'FontWeight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure 
% subplot(2,3,1)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_phi_TM_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\phi| [V/m]')
% 
% subplot(2,3,2)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_theta_TM_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1]) 
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_\theta| [V/m]')
% 
% subplot(2,3,3)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pE_pho_TM_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|E_r| [V/m]')
% 
% subplot(2,3,4)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_phi_TM_xz).*ro_mask_int)
% colormap jet
% colorbar 
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\phi| [μT]')
% 
% 
% subplot(2,3,5)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_theta_TM_xz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_\theta| [μT]')
% 
% subplot(2,3,6)
% imagesc([-FOV/2 FOV/2],[-FOV/2 FOV/2],abs(pH_pho_TM_xz).*ro_mask_int)
% colormap jet
% colorbar
% pbaspect([1 1 1])
% hold on
% for q=1:layer_number
% plot(x_a(q,:),y_a(q,:),'k--',LineWidth=1)
% end
% hold off
% xlabel('Z [m]','FontSize',10,'FontWeight','bold')
% ylabel('X [m]','FontSize',10,'FontWeight','bold')
% title('|B_r| [μT]')
% 
% sgtitle('TM Fields ZX Plane','Fontsize',12,'FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % Salva tutte le figure aperte
% % figHandles = findall(0, 'Type', 'figure');  % Trova tutte le figure
% % for i = 1:length(figHandles)
% %     fig = figHandles(i);
% %     filename = sprintf('Figure_%d.png', i);  % Nome del file
% %     saveas(fig, filename);  % Salva la figura
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

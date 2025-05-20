%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% LEGENDRE ASSOCIATED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is an example of Legendre Function as definied in the reference
% article. Useful to check and debug.

% Reference: On the Computation af Derivatives of Legendre Functions (W.Bosh 2000)
%author: Vincenzo Miranda
%Date: 20/05/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

clc;
close all;
clear all;

% Parametri
n_tot = 5; % Numero totale di ordini n
angle = linspace(0.1, pi, 1000); % Gamma angolare
kr = 0.12; % Parametro radiale
phi = 0; % Fase
pig0=0*ones(1,length(angle));
pig1=1*ones(1,length(angle));
Pig=zeros(n_tot,n_tot+1,length(angle));

% Ciclo su n
for n = 1:n_tot
    figure; % Nuova figura per ogni n
    legend_entries = {}; % Inizializza la legenda
    
    % Ciclo su m
    for m = 0:n
        % Calcola le funzioni associate di Legendre
        Pnm = legendre(n, cos(angle));

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if m==0
            Pig(n,m+1,:)=0*ones(1,length(angle));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif m==1

        if n==1
        Pig(n,m+1,:)=pig1.*ones(1,length(angle));
        elseif n==2
        Pig(n,m+1,:)= (2*n-1)/(n-1)*cos(angle).*pig1 - n/(n-1)*pig0;
        else
        PiG11=squeeze(Pig(n-1,m+1,:))';
        PiG22=squeeze(Pig(n-2,m+1,:))';
        Pig(n,m+1,:)= (2*n-1)/(n-1)*cos(angle).*PiG11 - n/(n-1)*PiG22;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else 
         Pig(n,m+1,:) = m*Pnm(m+1,:)./ sin(angle); % Funzione settoriale   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gestione del caso speciale per m=0
        if m == 0
            Pn1 = Pnm(m+2, :); % Legendre Pn(m+1)
            Tau = -Pn1; % Tau per m=0

        % Caso limite m=n
        elseif m == n
            P_nm_minus = Pnm(m, :); % Pnm(m-1)
            Tau = n * P_nm_minus; % Relazione ricorsiva
        
        % Caso generale 0 < m < n
        else
            P_nm_plus = Pnm(m+2, :); % Pn(m+1)
            P_nm_minus = Pnm(m, :); % Pn(m-1)
            Tau = 0.5 * ((n+m)*(n-m+1)*P_nm_minus - P_nm_plus); % Ricorsione
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calcolo di N_r
        % M3_phi=   (-spherbessH(n,1,kr)).*(exp(1i.*m.*phi)).*(Tau);  
        % M4_phi=   (-spherbessH(n,2,kr)).*(exp(1i.*m.*phi)).*(Tau);
        % M3_theta= (spherbessH(n,1,kr)).*(exp(1i.*m.*phi)).*((1i.*piG));
        % M4_theta= (spherbessH(n,2,kr)).*(exp(1i.*m.*phi)).*((1i.*piG)); 


        %N_theta= (derRspherbessJ(n,kr)).*(exp(1i.*m.*phi)).*(Tau);   
        N_phi=   (derRspherbessJ(n,kr)).*((exp(1i.*m.*phi)).*(1i.*Pig));    
        %N_pho = (spherbessJ(n, kr) / kr) .* exp(1i * m * phi) .* (n * (n+1)) .* Pnm(m+1, :);


        % Traccia il modulo di N_r
        hold on;
        plot(cos(angle), squeeze(abs(N_phi(n,m+1,:))), 'LineWidth', 2);
        
        % Aggiungi la descrizione alla legenda
        legend_entries{end+1} = ['m = ' num2str(m)];

    end
    
    % Impostazioni della figura
    axis tight;
    title(['Funzioni N_\theta per n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('|N_\theta(\theta)|');
    legend(legend_entries, 'Location', 'best');
    grid on;
    hold off;
end
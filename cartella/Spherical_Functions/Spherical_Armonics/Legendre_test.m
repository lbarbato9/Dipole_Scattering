%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% LEGENDRE ASSOCIATED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is an example of Legendre Function as definied in the reference
% article. Useful to check and debug.

% Reference: On the Computation af Derivatives of Legendre Functions (W.Bosh 2000)
%author: Vincenzo Miranda
%Date: 20/05/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %LEGENDRE ASSOCIATED FUNCTION TEST

n_tot = 5;  % Numero totale di ordini
angle = linspace(0, pi, 1000);  % Intervallo di angoli da 0 a pi

% Loop per ogni n (grado del polinomio di Legendre)
for n = 1:n_tot

    figure;  % Crea una nuova figura per ogni n
    
    % Array per contenere i nomi delle legende
    legend_entries = {};  
    
    % Loop per ogni m (ordine del polinomio di Legendre)
    for m = 0:n
                
        % Calcola le funzioni associate di Legendre per ogni m
        Pnn = legendre(n, cos(angle)); 
        
        % Estrai il corrispondente Pnm per questo m (sarà in Pnn(m+1,:,:))
        Pnm = Pnn(m+1, :);  % Estrai la riga appropriata corrispondente a m
        
        % Aggiungi il valore di m alla legenda
        legend_entries{end+1} = ['m = ' num2str(m)];
        
        % Plot della funzione Pnm in funzione dell'angolo
        hold on;
        plot(cos(angle), Pnm,LineWidth=2);
    end
    
    % Impostazioni degli assi e del titolo
    axis tight;
    title(['Legendre Associated Functions Pnm per n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('Pnm(cos(\theta))');
    
    % Aggiungi la legenda
    legend(legend_entries, 'Location', 'best');
    
    hold off;
end


%Verifiche eseguite:
%Pnm(+-1)=0        per ogni m>0 nei valori theta=0 e pi deve valere 0 [OK]
%Pn0(+1)=1         [OK]
%Pn0(-1)=(-1)^n    n=0->1  n=1->-1 n=2->1   [OK]     

%NB dunque ai poli contribuisce solo m=0 per ogni n dal momento che il
%resto delle funzioni di legendre assume valore nullo.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %LEGENDRE FIRST DERIVATIVE ASSOCIATED FUNCTION


% Function for correct evaluation of Legendre associated function
% derivative fixed degree n and order m.
clc
close all
clear all

n_tot=5;
angle=linspace(0,pi,1000);


for n=1:n_tot

       figure;  % Crea una nuova figura per ogni n
    
    % Array per contenere i nomi delle legende
    legend_entries = {};  
    
    for m=0:n
        
        % Calcola le funzioni associate di Legendre per ogni m
        Pnm = legendre(n, cos(angle)); 
        
        % % Estrai il corrispondente Pnm per questo m (sarà in Pnn(m+1,:,:))
        % Pnm = Pnn(m+1, :);  % Estrai la riga appropriata corrispondente a m
        % 
        % Aggiungi il valore di m alla legenda
        legend_entries{end+1} = ['m = ' num2str(m)];
        

%Special case for m=0
if m==0
  Pn1(:,:)=Pnm(m+2,:,:);
  Tau(:,:)=(-Pn1(:,:));

%Evaluation of Derivative Pnm for m>0 
elseif m==n
   P_nm_minus(:,:)=Pnm(m,:,:);  %Pn m-1   
   Tau=n*P_nm_minus ; %Recursive relation for Der Pnm     

  elseif m<n 
  P_nm_plus(:,:)=Pnm(m+2,:,:); %Pn m+1 
  P_nm_minus(:,:)=Pnm(m,:,:);  %Pn m-1   
  Tau=0.5*( (n+m)*(n-m+1)*P_nm_minus - P_nm_plus ); %Recursive relation for Der Pnm 
 end


  % Plot della funzione Pnm in funzione dell'angolo
        hold on;
        plot(cos(angle), Tau,LineWidth=2);

end 

    % Impostazioni degli assi e del titolo
    axis tight;
    title(['Functions \tau per n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('\tau (\theta)');
    
    % Aggiungi la legenda
    legend(legend_entries, 'Location', 'best');
    hold off;

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %PIGRECO FUNCTION
                            
n_tot = 3;  % Numero totale di ordini
angle = linspace(0, pi, 1000);  % Intervallo di angoli da 0 a pi
pig0=0*ones(1,length(angle));
pig1=1*ones(1,length(angle));
Pig=zeros(n_tot,n_tot+1,length(angle));

% Loop per ogni n (grado del polinomio di Legendre)
for n = 1:n_tot

    figure;  % Crea una nuova figura per ogni n
    
    % Array per contenere i nomi delle legende
    legend_entries = {};  
    
    % Loop per ogni m (ordine del polinomio di Legendre)
    for m = 0:n
                
        % Calcola le funzioni associate di Legendre per ogni m
        Pnn = legendre(n, cos(angle)); 
        
        % Estrai il corrispondente Pnm per questo m (sarà in Pnn(m+1,:,:))
        Pnm = Pnn(m+1, :);  % Estrai la riga appropriata corrispondente a m
        
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
         Pig(n,m+1,:) = m*Pnm./ sin(angle); % Funzione settoriale   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        % Aggiungi il valore di m alla legenda
        legend_entries{end+1} = ['m = ' num2str(m)];
        
        % Plot della funzione Pnm in funzione dell'angolo
        hold on;
        plot(cos(angle), squeeze(Pig(n,m+1,:)) ,LineWidth=2);
    end
    
    % Impostazioni degli assi e del titolo
    axis tight;
    title(['Functions \pi per n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('\pi (\theta))');
    
    % Aggiungi la legenda
    legend(legend_entries, 'Location', 'best');
    hold off;
    xlim([-1 1])
end
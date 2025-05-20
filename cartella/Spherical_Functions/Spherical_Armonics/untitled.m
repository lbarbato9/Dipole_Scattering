clc
close all
clear all

theta=linspace(0,pi,2000);
n_tot=5;

P=zeros(n_tot,n_tot+1,length(theta));
Pig=zeros(n_tot,n_tot+1,length(theta));
tau=zeros(n_tot,n_tot+1,length(theta));

for n=1:n_tot

for m=0:n

P(n,m+1,:)=Pnm(n,m,theta);
Pig(n,m+1,:)=Pigreco(n,m,theta);
tau(n,m+1,:)=Tau(n,m,theta);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEGENDRE POLINOMIALS PLOT

for n=1:n_tot
    figure
    legend_entries = {};  

for m=0:n

hold on     
plot(cos(theta),squeeze(P(n,m+1,:)))
legend_entries{end+1} = ['m = ' num2str(m)];

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAU FUNCTION PLOT


for n=1:n_tot
    figure
    legend_entries = {};  

for m=0:n

hold on     
plot(cos(theta),squeeze(tau(n,m+1,:)))
legend_entries{end+1} = ['m = ' num2str(m)];

end

 % Impostazioni degli assi e del titolo
    axis tight;
    title(['\tau function for n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('\tau(\theta)');
    
    % Aggiungi la legenda
    legend(legend_entries, 'Location', 'best');
    
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIGRECO FUNCTION PLOT


for n=1:n_tot
    figure
    legend_entries = {};  

for m=0:n

hold on     
plot(cos(theta),squeeze(Pig(n,m+1,:)))
legend_entries{end+1} = ['m = ' num2str(m)];

end

 % Impostazioni degli assi e del titolo
    axis tight;
    title(['\pi function for n = ', num2str(n)]);
    xlabel('cos(\theta)');
    ylabel('\pi(\theta)');
    
    % Aggiungi la legenda
    legend(legend_entries, 'Location', 'best');
    
    hold off;
end
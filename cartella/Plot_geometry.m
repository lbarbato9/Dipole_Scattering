%% Parametri iniziali
% layer_number e a sono definiti altrove

% Creazione della sfera
[X, Y, Z] = sphere(50);  % Sfera con risoluzione di 50 punti

% Opacità crescente per i vari strati
alphas = linspace(0.8, 0.2, layer_number);  % Maggiore trasparenza per strati interni

% Definizione della scala di grigi per gli strati
grayscale = linspace(0.1, 0.5, layer_number);  % Scala di grigio dal più scuro al più chiaro

% Visualizzazione 3D interattiva
figure;
hold on;  % Mantiene il grafico attuale

% Ciclo per creare sfere concentriche
for i = 1:layer_number
    raggio_sfera = a(i); % Raggio crescente per ogni strato

    % Colore grigio per la sfera corrente
    color = [grayscale(i), grayscale(i), grayscale(i)];  % Colore in scala di grigi

    % Disegna la sfera con il raggio corrispondente e il colore
    surf(X*raggio_sfera, Y*raggio_sfera, Z*raggio_sfera, 'FaceColor', color, 'FaceAlpha', alphas(i), 'EdgeColor', 'none');  
end

% Parametri del dipolo
dipole_length = 0.1;  % Lunghezza del dipolo
P_dipole = (r_d+dipole_length/2) * dipole_moment / norm(dipole_moment);  % Punto di posizionamento del dipolo nello spazio

% Posizione delle cariche del dipolo in base al vettore dipole_moment
dipole_moment = dipole_moment / norm(dipole_moment);  % Normalizza il vettore
direction_vector = (dipole_length / 2) * dipole_moment;  % Vettore direzione del dipolo
charge1 = P_dipole - direction_vector;  % Posizione carica negativa
charge2 = P_dipole + direction_vector;  % Posizione carica positiva

% Disegna il dipolo come un segmento di linea
plot3([charge1(1), charge2(1)], [charge1(2), charge2(2)], [charge1(3), charge2(3)], 'k-', 'LineWidth', 2);

% Aggiunge punti per rappresentare le cariche
plot3(charge1(1), charge1(2), charge1(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Carica negativa
plot3(charge2(1), charge2(2), charge2(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Carica positiva

% Impostazioni del grafico
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Spherical Geometry with Dipole');
view(3);  % Imposta la vista 3D
rotate3d on;  % Abilita rotazione interattiva
grid on;
hold off;  % Rilascia il grafico attuale

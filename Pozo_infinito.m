clc;
clear;
close all;

%% --- Configuración de estilo global para las representaciones gráficas ---
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 12);

%% =========================================================================
% SECCIÓN 1: CONVERGENCIA DEL ERROR
% =========================================================================
% Evalúa el error relativo de los tres primeros autovalores (E1, E2, E3) frente a la solución analítica exacta.

disp('Calculando convergencia del error en escala Log-Log...');
N_vec = round(logspace(1, 3.3, 30)); % N desde 10 hasta 2000 en escala log
err_E1 = zeros(size(N_vec));
err_E2 = zeros(size(N_vec));
err_E3 = zeros(size(N_vec));

for idx = 1:length(N_vec)
    k = N_vec(idx);
    delta = 1 / (k + 1);
    
    % Ensamblaje vectorizado con sparse
    e = ones(k, 1);
    H = spdiags([-e, 2*e, -e], [-1, 0, 1], k, k) / (delta^2);
    
    % Extrae solo los autovalores necesarios con eigs
    num_eigs = min(10, k); 
    if k < 20
        E_num = sort(eig(full(H)));
    else
        opts.isreal = true; 
        E_num = sort(eigs(H, num_eigs, 'smallestabs', opts));
    end
    
    % Se normaliza dividiendo por pi^2 (para comparar con n^2 exacto)
    E_num = E_num / (pi^2);
    
    % Errores relativos |(E_exa - E_num) / E_exa|
    err_E1(idx) = abs(1 - E_num(1) / (1^2));
    if k >= 2
        err_E2(idx) = abs(1 - E_num(2) / (2^2));
    else
        err_E2(idx) = NaN;
    end

    if k >= 3
        err_E3(idx) = abs(1 - E_num(3) / (3^2));    
    else
        err_E3(idx) = NaN;
    end
end

figure('Name', 'Convergencia Logaritmica');
loglog(N_vec, err_E1, 'DisplayName', 'E_1');
hold on;
loglog(N_vec, err_E2, 'DisplayName', 'E_{2}');
loglog(N_vec, err_E3, 'DisplayName', 'E_{3}');
grid on;
xlabel('N');
ylabel('Error relativo');
% title('Convergencia del error de discretización espacial');
legend('Location', 'best');
% saveas(gcf, 'Convergencia_Error_Log_v2.png');


%% =========================================================================
% SECCIÓN 2: ANÁLISIS DEL ESPECTRO COMPLETO
% =========================================================================
% Esta sección evalúa qué fracción de ese espectro numérico es físicamente útil (error < 1% o < 10%)

disp('Analizando ceguera numérica para todo el espectro...');
N_max = 1000;
f_1 = zeros(1, N_max-1);
f_10 = zeros(1, N_max-1);

for k = 2:N_max
    delta = 1 / (k + 1);
    e = ones(k, 1);
    H = spdiags([-e, 2*e, -e], [-1, 0, 1], k, k) / (delta^2);
    
    E_num = eig(full(H)) / (pi^2);
    E_exa = (1:k)'.^2;
    
    dif = abs(1 - E_num ./ E_exa);
    
    f_1(k-1) = sum(dif < 0.01) / k; % Estados con error < 1%
    f_10(k-1) = sum(dif < 0.10) / k; % Estados con error < 10%
end

figure('Name', 'Fraccion Estados Validos');
plot(2:N_max, f_1, 'b-', 'DisplayName', 'Tolerancia < 1%');
hold on;
plot(2:N_max, f_10, 'r-', 'DisplayName', 'Tolerancia < 10%');
grid on;
xlabel('N');
ylabel('Fracción de espectro recuperable');
%title('Ceguera numérica: Fracción de estados precisos');
legend('Location', 'best');
ylim([0 0.5]);
%saveas(gcf, 'Fraccion_Estados_Validos_v2.png');

% --- Gráfica de por qué falla la ceguera numérica (Error vs Nivel n) ---
figure('Name', 'Error Altas Frecuencias');

plot(1:k, E_exa, 'k-', 'DisplayName', 'Exacta');
hold on;
plot(1:k, E_num, 'b--', 'DisplayName', 'Numérica');
grid on;
xlabel('n');
ylabel('E_n / \pi^2');
%title('Divergencia espectral ($N=1000$ fijo)');
legend('Location', 'northwest');

figure('Name', 'Energías n')
plot(1:k, dif, 'r-');
grid on;
xlabel('n');
ylabel('Error relativo');
%title('Disparada del error por truncamiento');

%% =========================================================================
% SECCIÓN 3: SUPERPOSICIÓN GEOMÉTRICA DE AUTOESTADOS ESPACIALES
% =========================================================================
% Grafica la amplitud de probabilidad de las 3 primeras funciones de onda, superponiendo la solución discreta numérica sobre la curva analítica exacta continua

disp('Generando autoestados espaciales...');
k = 100;
delta = 1 / (k + 1);
x = linspace(delta, 1-delta, k)'; % Nodos espaciales
e = ones(k, 1);
H = spdiags([-e, 2*e, -e], [-1, 0, 1], k, k) / (delta^2);

opts.isreal = true; 
[V, D] = eigs(H, 3, 'smallestabs', opts);

[~, sort_idx] = sort(diag(D));
V = V(:, sort_idx);

x_exact = linspace(0, 1, 500)';

figure('Name', 'Autoestados 1D');
colores = lines(3);

for n = 1:3
    
    % Solución analítica exacta
    psi_exact = sqrt(2) * sin(n * pi * x_exact);
    
    % Solución numérica
    psi_num = V(:, n);
    
    % Normalizar L2 discreto: sum(psi_num^2 * delta) = 1
    psi_num = psi_num / sqrt(sum(psi_num.^2 * delta));
    
    % Corregir signo arbitrario del autovector si está invertido
    if sign(psi_num(2)) ~= sign(psi_exact(2))
        psi_num = -psi_num;
    end
    
    % Plot exacto vs numérico
    plot(x_exact, psi_exact, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    hold on;
    plot(x, psi_num, 'LineStyle', '--', 'Color', colores(n,:), 'DisplayName', sprintf('\\Psi_%d', n));    
end

plot(NaN, NaN, 'k-', 'LineWidth', 1.0, 'DisplayName', '\Psi_i^{exacta}');
%title(sprintf('Estado $n=%d$', n));
xlabel('x [u.a.]');
ylabel('\Psi(x)');
grid on;
xlim([0 1]);
legend('Location', 'northeast');
%sgtitle('Superposición Geométrica de Autoestados ($N=100$)');
%saveas(gcf, 'Autoestados_Superposicion_v2.png');

disp('Proceso completado con éxito. Gráficas generadas.');
clc;
clear;
close all;

%% --- Configuración de estilo global para las representaciones gráficas ---
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 12);
colores = lines(4);

%% =========================================================================
% SECCIÓN 1: CONVERGENCIA Y COSTE COMPUTACIONAL
% =========================================================================
% Construye la matriz de un pozo 2D mediante el producto de Kronecker. 
% Cuantifica el error relativo de las energías y monitorea el tiempo de CPU
% P_vec representa los puntos N por cada dimensión (matriz real de N^2 x N^2)
P_vec = round(logspace(1, 3, 100)); % Puntos desde 10 hasta 1000

err_11 = zeros(size(P_vec));
err_12 = zeros(size(P_vec));
err_22 = zeros(size(P_vec));
t_cpu = zeros(size(P_vec)); % Tiempo de cómputo (CPU) 2D
err_1D = zeros(size(P_vec));
t_cpu_1D = zeros(size(P_vec)); % Tiempo de cómputo (CPU) 1D

% Energías exactas [u.a.]
E_exa_11 = 1^2 + 1^2; % = 2
E_exa_12 = 1^2 + 2^2; % = 5 (Doble degenerado con 2,1)
E_exa_22 = 2^2 + 2^2; % = 8

disp('Calculando autovalores. Observa cómo el tiempo aumenta drásticamente...');

for idx = 1:length(P_vec)
    N = P_vec(idx);
    delta = 1 / (N + 1);
    
    % Matriz de diferencias finitas 1D
    e = ones(N, 1);
    D2 = spdiags([-e, 2*e, -e], [-1, 0, 1], N, N) / delta^2;
    I = speye(N);
    
    % --- POZO 1D ---
    tic;
    E_num_1D = sort(eigs(D2, 1, 'sm')) / pi^2;
    t_cpu_1D(idx) = toc;
    err_1D(idx) = abs(1 - E_num_1D(1)/1^2); % Exacta 1D es 1^2
    
    % --- POZO 2D ---
    tic; % Iniciar cronómetro 2D
    
    % Operador Laplaciano 2D mediante producto de Kronecker
    H = kron(D2, I) + kron(I, D2); 
    
    % H es de tamaño (N^2) x (N^2)
    % Se usa eigs para buscar solo los 4 estados más bajos y evitar colapso de RAM
    E_num = sort(eigs(H, 4, 'sm')) / pi^2; 
    
    t_cpu(idx) = toc; % Detener cronómetro
    
    err_11(idx) = abs(1 - E_num(1)/E_exa_11);
    err_12(idx) = abs(1 - E_num(2)/E_exa_12); % E_num(2) y E_num(3) son degenerados (1,2 y 2,1)
    err_22(idx) = abs(1 - E_num(4)/E_exa_22); % El 4º estado es el (2,2)
end

figure('Name', 'Convergencia Pozo 2D');
loglog(P_vec, err_11, 'Color', colores(1,:), 'DisplayName', 'E_{1,1}');
hold on;
loglog(P_vec, err_12, 'Color', colores(2,:), 'DisplayName', 'E_{1,2}');
loglog(P_vec, err_22, 'Color', colores(3,:), 'DisplayName', 'E_{2,2}');
grid on;
xlabel('N');
ylabel('Error Relativo');
legend('Location', 'northeast');
% title('Convergencia en 2D: El inicio de la maldición');

% Gráfica de Coste Computacional (Eficiencia 1D vs 2D)
figure('Name', 'Eficiencia: 1D vs 2D');
loglog(err_1D, t_cpu_1D, 'b-', 'LineWidth', 2, 'DisplayName', 'Pozo 1D');
hold on;
loglog(err_11, t_cpu, 'r-', 'LineWidth', 2, 'DisplayName', 'Pozo 2D');
grid on;
xlabel('Error Relativo');
ylabel('t [s]');
set(gca, 'XDir', 'reverse'); % Invertir eje X para que "más precisión" (menor error) avance hacia la derecha
legend('Location', 'northwest');
% title('Eficiencia Computacional: El coste de la precisión');

%% =========================================================================
% SECCIÓN 2: MAPAS DE CONTORNO DE AUTOESTADOS 2D
% =========================================================================
% Calcula los autovectores espaciales y proyecta el subespacio degenerado.
% Representa los contornos de amplitud (FDM) superpuestos sobre las líneas teóricas exactas
disp('Generando superposición rigurosa de mapas de contorno...');
N_wave = 80; % Resolución visual
delta = 1 / (N_wave + 1);
e = ones(N_wave, 1);
D2 = spdiags([-e, 2*e, -e], [-1, 0, 1], N_wave, N_wave) / delta^2;
I = speye(N_wave);
H_wave = kron(D2, I) + kron(I, D2);

[V_mat, D_mat] = eigs(H_wave, 4, 'sm');
[~, sort_idx] = sort(diag(D_mat));
V_mat = V_mat(:, sort_idx);

% Remapear vectores columna 1D a matrices 2D (Espacio de configuración)
x = linspace(delta, 1-delta, N_wave);
[X, Y] = meshgrid(x, x);

Psi_11 = reshape(V_mat(:,1), N_wave, N_wave);
% Forzar fase positiva
if Psi_11(round(N_wave/2), round(N_wave/2)) < 0, Psi_11 = -Psi_11; end

% Usamos el estado excitado degenerado (1,2)
% Al ser degenerado, eigs devolverá una combinación lineal arbitraria de (1,2) y (2,1).
Psi_12 = reshape(V_mat(:,2), N_wave, N_wave);

% Normalizar ondas numéricas para que coincidan con la amplitud teórica
% eigs devuelve sum(V.^2) = 1, la física exige sum(|Psi|^2 * delta^2) = 1
Psi_11_norm = Psi_11 / delta;
Psi_12_norm = Psi_12 / delta;

% Estado fundamental analítico exacto
Psi_exa_11 = 2 * sin(pi*X) .* sin(pi*Y);

% Proyección Rigurosa del estado degenerado (E=5)
Psi_base_1 = 2 * sin(pi*X) .* sin(2*pi*Y);
Psi_base_2 = 2 * sin(2*pi*X) .* sin(pi*Y);

% Producto escalar continuo: <Psi_num | Psi_base>
alpha = sum(sum(Psi_12_norm .* Psi_base_1)) * delta^2;
beta  = sum(sum(Psi_12_norm .* Psi_base_2)) * delta^2;

% Estado mixto analítico exacto correspondiente a la convergencia
Psi_exa_12 = alpha * Psi_base_1 + beta * Psi_base_2;

figure('Name', 'Contornos FDM vs Exacta (1,1)');
hold on;
contourf(X, Y, Psi_11_norm, 10, 'LineStyle', 'none'); 
colormap('jet');
contour(X, Y, Psi_exa_11, 10, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Teórica');
axis square; grid on;
% title('Estado (1,1): Numérico vs Teórico');
xlabel('x_1 [u.a.]'); ylabel('x_2 [u.a.]');

figure('Name', 'Contornos FDM vs Exacta (1,2)');
hold on;
contourf(X, Y, Psi_12_norm, 10, 'LineStyle', 'none'); 
colormap('jet');
contour(X, Y, Psi_exa_12, 10, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Teórica');
axis square; grid on;
% title('Estado Excitado (Mezcla): Numérico vs Teórico');
xlabel('x_1 [u.a.]'); ylabel('x_2 [u.a.]');

%% =========================================================================
% SECCIÓN 3: LA MALDICIÓN DE LA DIMENSIONALIDAD
% =========================================================================
% Proyecta de forma teórica el crecimiento hiperdimensional de las matrices necesarias para extrapolar el método a D=3 y D=6
disp('Proyectando crecimiento de memoria en alta dimensionalidad...');

colores = [
    0.0, 0.0, 1.0;  % 1D: Azul puro ('b')
    1.0, 0.0, 0.0;  % 2D: Rojo puro ('r')
    0.9, 0.7, 0.0;  % 3D: Amarillo oro (Alta saturación, pero visible en fondo blanco)
    0.6, 0.0, 0.8   % 6D: Morado intenso y oscuro
];

N_hardware = logspace(0, 2, 50); % N desde 1 hasta 100
D_vec = [1, 2, 3, 6]; % 1D, 2D, 3D, 6D (Helio)
labels_D = {'1D', '2D', '3D', '6D'};

figure('Name', 'Explosion Dimensional');
for i = 1:length(D_vec)
    D = D_vec(i);
    Tamanio_Matriz = N_hardware.^(2 * D); 
    
    loglog(N_hardware, Tamanio_Matriz, 'LineWidth', 2, 'Color', colores(i,:), 'DisplayName', labels_D{i});
    hold on;
end

yline(2e9, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Límite 16 GB RAM');
grid on;
xlabel('N');
ylabel('N^{2D}'); 
legend('Location', 'northwest');
% title('Colapso del FDM: La maldición de la dimensionalidad');
xlim([1 100]);
ylim([1e1 1e25]); 

disp('Proceso completado con éxito. Gráficas bidimensionales generadas.');
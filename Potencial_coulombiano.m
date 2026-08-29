clc;
clear;
close all;

%% --- Configuración de estilo global para las representaciones gráficas ---
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 12);
colores = lines(6); 
colores(3,:) = [0.8500, 0.6200, 0.1000]; % Ajuste de contraste

d_0 = 100; % Dominio de simulación macroscopicamente grande para evitar error de frontera en los primeros estados

%% =========================================================================
% SECCIÓN 1: CONVERGENCIA PARAMÉTRICA Y RUPTURA DE DEGENERACIÓN
% =========================================================================
% Genera el análisis combinado del error de truncamiento
% El primer panel constata la dilución espacial
% El segundo panel expone la ruptura computacional de la degeneración en la capa n=3. Además, imprime por terminal la Tabla 2.1
disp('Calculando convergencia paramétrica y tabla de energías...');

N_vec_def = round(logspace(1.5, 3.5, 40));

% Matrices de error
err_1s = zeros(size(N_vec_def));
err_2s = zeros(size(N_vec_def));
err_3s = zeros(size(N_vec_def));
err_8s = zeros(size(N_vec_def));

err_3p = zeros(size(N_vec_def));
err_3d = zeros(size(N_vec_def));

% Variables para tabla de energías crudas
E_val_3s = zeros(size(N_vec_def));
E_val_3p = zeros(size(N_vec_def));
E_val_3d = zeros(size(N_vec_def));

% Energías exactas [u.a.]
E_exa_1 = -1/(2*1^2);
E_exa_2 = -1/(2*2^2);
E_exa_3 = -1/(2*3^2);
E_exa_8 = -1/(2*8^2); % Nuevo estado para demostrar error de frontera

for idx = 1:length(N_vec_def)
    N = N_vec_def(idx);
    delta = d_0 / (N + 1);
    rho = linspace(delta, d_0 - delta, N)';
    e = ones(N, 1);
    T = spdiags([-e, 2*e, -e], [-1, 0, 1], N, N) / (2*delta^2);
    
    % --- MATRIZ l=0 ---
    V_0 = spdiags(0*(1)./(2*rho.^2) - 1./rho, 0, N, N);
    H_0 = (T + V_0 + (T + V_0)') / 2;
    if N < 20
        E_0 = sort(eig(full(H_0)));
    else
        % Se buscan los 8 primeros para alcanzar el 8s
        try
            E_0 = sort(eigs(H_0, 8, -0.5)); 
        catch
            % Respaldo por si eigs falla al buscar tantos con N pequeño
            E_0 = sort(eig(full(H_0)));
        end
    end
    err_1s(idx) = abs(1 - E_0(1)/E_exa_1);
    err_2s(idx) = abs(1 - E_0(2)/E_exa_2);
    err_3s(idx) = abs(1 - E_0(3)/E_exa_3);
    if length(E_0) >= 8
        err_8s(idx) = abs(1 - E_0(8)/E_exa_8);
    else
        err_8s(idx) = NaN; % Evitar errores si eig devuelve menos de 8
    end
    E_val_3s(idx) = E_0(3);
    
    % --- MATRIZ l=1 ---
    V_1 = spdiags(1*(2)./(2*rho.^2) - 1./rho, 0, N, N);
    H_1 = (T + V_1 + (T + V_1)') / 2;
    if N < 20
        E_1 = sort(eig(full(H_1)));
    else
        E_1 = sort(eigs(H_1, 2, -0.125)); % Se buscan los 2 primeros (2p, 3p)
    end
    err_3p(idx) = abs(1 - E_1(2)/E_exa_3); % 3p es el segundo estado de l=1
    E_val_3p(idx) = E_1(2);
    
    % --- MATRIZ l=2 ---
    V_2 = spdiags(2*(3)./(2*rho.^2) - 1./rho, 0, N, N);
    H_2 = (T + V_2 + (T + V_2)') / 2;
    if N < 20
        E_2 = sort(eig(full(H_2)));
    else
        E_2 = sort(eigs(H_2, 1, -0.0555)); % Se busca el primero (3d)
    end
    err_3d(idx) = abs(1 - E_2(1)/E_exa_3); % 3d es el primer estado de l=2
    E_val_3d(idx) = E_2(1);
end

% --- IMPRESIÓN DE TABLA (N ≈ 298) ---
[~, idx_tab] = min(abs(N_vec_def - 298));
N_tab = N_vec_def(idx_tab);
disp(' ');
fprintf('--- TABLA DE ENERGÍAS PARA n=3 (N = %d) ---\n', N_tab);
fprintf('Exacta:    %10.6f u.a.\n', E_exa_3);
fprintf('3s (l=0):  %10.6f u.a. | Error: %8.4f %%\n', E_val_3s(idx_tab), err_3s(idx_tab));
fprintf('3p (l=1):  %10.6f u.a. | Error: %8.4f %%\n', E_val_3p(idx_tab), err_3p(idx_tab));
fprintf('3d (l=2):  %10.6f u.a. | Error: %8.4f %%\n', E_val_3d(idx_tab), err_3d(idx_tab));
disp('---------------------------------------------');
disp(' ');


figure('Name', 'Convergencia l = 0');

% Panel A: l=0 fijo, variando n (1s, 2s, 3s y 8s)
loglog(N_vec_def, err_1s, 'Color', colores(1,:), 'LineWidth', 2, 'DisplayName', '1s');
hold on;
loglog(N_vec_def, err_2s, 'Color', colores(2,:), 'LineWidth', 2, 'DisplayName', '2s');
loglog(N_vec_def, err_3s, 'Color', colores(3,:), 'LineWidth', 2, 'DisplayName', '3s');
loglog(N_vec_def, err_8s, 'Color', colores(4,:), 'LineWidth', 2, 'DisplayName', '8s');
grid on;
xlabel('N');
ylabel('Error Relativo');
% title('Efecto de la contracción espacial (l=0)');
legend('Location', 'southwest');

figure('Name', 'Convergencia n = 3');

% Panel B: n=3 fijo, variando l (3s, 3p, 3d)
loglog(N_vec_def, err_3s, 'Color', colores(3,:), 'LineWidth', 2, 'DisplayName', 'l=0'); % Mismo color que en Panel A
hold on;
loglog(N_vec_def, err_3p, 'Color', colores(5,:), 'LineWidth', 2, 'DisplayName', 'l=1');
loglog(N_vec_def, err_3d, 'Color', colores(6,:), 'LineWidth', 2, 'DisplayName', 'l=2');
grid on;
xlabel('N');
ylabel('Error Relativo');
% title('Ruptura de la degeneración (n=3)');
legend('Location', 'northeast');


%% =========================================================================
% SECCIÓN 2: AUTOESTADOS CON DISTINTOS MOMENTOS ANGULARES
% =========================================================================
% Representa las funciones de onda para n=3 y distintos l (3s, 3p, 3d) cerca del núcleo
disp('Calculando estructura radial de la capa n=3...');

N_wave = 1000;
delta_w = d_0 / (N_wave + 1);
rho_w = linspace(delta_w, d_0 - delta_w, N_wave)';
e_w = ones(N_wave, 1);
T_w = spdiags([-e_w, 2*e_w, -e_w], [-1, 0, 1], N_wave, N_wave) / (2*delta_w^2);

% Función exacta para la onda radial reducida u_{nl}(rho) usando Laguerre
u_exacta = @(n, l, rho) rho .* sqrt((2/n)^3 * factorial(n-l-1) / (2*n*factorial(n+l))) .* ...
    (2*rho/n).^l .* exp(-rho/n) .* double(laguerreL(n-l-1, 2*l+1, 2*rho/n));

% --- 3s (n=3, l=0) ---
V_0 = spdiags(0 - 1./rho_w, 0, N_wave, N_wave);
H_0 = (T_w + V_0 + (T_w + V_0)') / 2;
[V_mat, D] = eigs(H_0, 3, -0.5);
[~, idx] = sort(diag(D));
u_num_3s = V_mat(:,idx(3)); % 3er autovector es el 3s
u_num_3s = u_num_3s / sqrt(sum(u_num_3s.^2 * delta_w));
if u_num_3s(2) < 0, u_num_3s = -u_num_3s; end
u_exa_3s = u_exacta(3, 0, rho_w);

% --- 3p (n=3, l=1) ---
V_1 = spdiags(2./(2*rho_w.^2) - 1./rho_w, 0, N_wave, N_wave);
H_1 = (T_w + V_1 + (T_w + V_1)') / 2;
[V_mat, D] = eigs(H_1, 2, -0.125);
[~, idx] = sort(diag(D));
u_num_3p = V_mat(:,idx(2)); % 2do autovector es el 3p
u_num_3p = u_num_3p / sqrt(sum(u_num_3p.^2 * delta_w));
if u_num_3p(2) < 0, u_num_3p = -u_num_3p; end
u_exa_3p = u_exacta(3, 1, rho_w);

% --- 3d (n=3, l=2) ---
V_2 = spdiags(6./(2*rho_w.^2) - 1./rho_w, 0, N_wave, N_wave);
H_2 = (T_w + V_2 + (T_w + V_2)') / 2;
[V_mat, D] = eigs(H_2, 1, -0.0555);
u_num_3d = V_mat(:,1);
u_num_3d = u_num_3d / sqrt(sum(u_num_3d.^2 * delta_w));
if u_num_3d(2) < 0, u_num_3d = -u_num_3d; end
u_exa_3d = u_exacta(3, 2, rho_w);

figure('Name', 'Repulsion_n3_Definitiva');
hold on;
plot(rho_w, u_exa_3s, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_3s, '--', 'Color', colores(3,:), 'LineWidth', 1.5, 'DisplayName', 'u_{3s}');

plot(rho_w, u_exa_3p, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_3p, '--', 'Color', colores(5,:), 'LineWidth', 1.5, 'DisplayName', 'u_{3p}');

plot(rho_w, u_exa_3d, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_3d, '--', 'Color', colores(6,:), 'LineWidth', 1.5, 'DisplayName', 'u_{3d}');

plot(NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', 'u^{exactas}_{3l}')

xlim([0 30]); 
grid on; 
legend('Location', 'northeast'); 
xlabel('r [u.a.]');
ylabel('u(r)');
% title('Ondas de la capa n=3 cerca del núcleo');

%% =========================================================================
% SECCIÓN 3: DILUCIÓN ESPACIAL DE LOS AUTOESTADOS SFERICOS
% =========================================================================
% Visualiza las funciones de onda 1s, 2s y 3s
disp('Calculando dilución espacial de los estados con l=0...');

% Se vuelve a diagonalizar H_0
[V_mat, D] = eigs(H_0, 3, -0.5);
[~, idx] = sort(diag(D));

% Extraer 1s y 2s
u_num_1s = V_mat(:,idx(1));
u_num_1s = u_num_1s / sqrt(sum(u_num_1s.^2 * delta_w));
if u_num_1s(2) < 0, u_num_1s = -u_num_1s; end
u_exa_1s = u_exacta(1, 0, rho_w);

u_num_2s = V_mat(:,idx(2));
u_num_2s = u_num_2s / sqrt(sum(u_num_2s.^2 * delta_w));
if u_num_2s(2) < 0, u_num_2s = -u_num_2s; end
u_exa_2s = u_exacta(2, 0, rho_w);

figure('Name', 'Efecto de dilución espacial');

hold on;
plot(rho_w, u_exa_1s, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_1s, '--', 'Color', colores(1,:), 'LineWidth', 1.5, 'DisplayName', 'u_{1s}');

plot(rho_w, u_exa_2s, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_2s, '--', 'Color', colores(2,:), 'LineWidth', 1.5, 'DisplayName', 'u_{2s}');

plot(rho_w, u_exa_3s, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(rho_w, u_num_3s, '--', 'Color', colores(3,:), 'LineWidth', 1.5, 'DisplayName', 'u_{3s}');

plot(NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', 'u^{exactas}_{n0}')

xlim([0 25]);
grid on;
legend('Location', 'northeast');
xlabel('r [u.a.]');
ylabel('u(r)');
% title('(a) Expansión radial (simetría s)');

disp('Proceso completado con éxito. Gráficas del hidrógeno generadas.');
clc;
clear;
close all;

%% --- Configuración de estilo global para las representaciones gráficas ---
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 12);

color_ex = 'b';
color_err = 'r';
colores = lines(3);
% Se oscurece el amarillo para que no se pierda al hacerlo línea discontinua
colores(3,:) = [0.8500, 0.6200, 0.1000];

%% =========================================================================
% SECCIÓN 1: CONVERGENCIA DEL ERROR FRENTE A N
% =========================================================================
% Cuantifica el error relativo de energía frente a la densidad de la malla N, para una anchura de caja fija
disp('Calculando convergencia del error en escala Log-Log (L fijo)...');
L_fixed = 10; % Un L razonable para albergar los primeros estados sin ahogarlos
N_vec = round(logspace(1, 3, 30)); 
err_E0 = zeros(size(N_vec));
err_E7 = zeros(size(N_vec));
err_E13 = zeros(size(N_vec));

E_exa_0 = 0.5;
E_exa_7 = 7.5;
E_exa_13 = 13.5;

for idx = 1:length(N_vec)
    N = N_vec(idx);
    delta = L_fixed / (N + 1);
    x = linspace(-L_fixed/2 + delta, L_fixed/2 - delta, N)'; % Nodos internos
    
    % Matriz Cinética T (Diferencias finitas)
    e = ones(N, 1);
    T = spdiags([-e, 2*e, -e], [-1, 0, 1], N, N) / (2 * delta^2);
    
    % Matriz Potencial V (Armónico)
    V = spdiags(0.5 * x.^2, 0, N, N);
    
    H = T + V;
    
    num_eigs = min(15, N);
    if N < 20
        E_num = sort(eig(full(H)));
    else
        E_num = sort(eigs(H, num_eigs, 'smallestabs'));
    end
    
    err_E0(idx) = abs(1 - E_num(1) / E_exa_0);
    if N >= 8
        err_E7(idx) = abs(1 - E_num(8) / E_exa_7);
    else
        err_E7(idx) = NaN;
    end
    if N >= 14
        err_E13(idx) = abs(1 - E_num(14) / E_exa_13);
    else
        err_E13(idx) = NaN;
    end
end

figure('Name', 'Convergencia Armonico');
loglog(N_vec, err_E0, 'Color', colores(1,:), 'DisplayName', 'E_0');
hold on;
loglog(N_vec, err_E7, 'Color', colores(2,:), 'DisplayName', 'E_7');
loglog(N_vec, err_E13, 'Color', colores(3,:), 'DisplayName', 'E_{13}');
grid on;
xlabel('N');
ylabel('Error rrelativo');
% title(sprintf('Convergencia del error (Dominio truncado fijo $L=%d$)', L_fixed));
legend('Location', 'northeast');
% saveas(gcf, 'Convergencia_Armonico_Log_v2.png');

%% =========================================================================
% SECCIÓN 2: OPTIMIZACIÓN DEL TAMAÑO DEL DOMINIO
% =========================================================================
% Evalúa el error de la energía frente a la anchura de la caja (L), manteniendo la cantidad de nodos N constante
disp('Calculando optimización del dominio L (N fijo)...');
N_fixed = 1000;
L_vec = linspace(2, 100, 50);
err_L_E0 = zeros(size(L_vec));
err_L_E7 = zeros(size(L_vec));
err_L_E13 = zeros(size(L_vec));

E_exa_3 = 3.5;

for idx = 1:length(L_vec)
    L = L_vec(idx);
    delta = L / (N_fixed + 1);
    x = linspace(-L/2 + delta, L/2 - delta, N_fixed)';
    
    e = ones(N_fixed, 1);
    T = spdiags([-e, 2*e, -e], [-1, 0, 1], N_fixed, N_fixed) / (2 * delta^2);
    V = spdiags(0.5 * x.^2, 0, N_fixed, N_fixed);
    H = T + V;
    
    E_num = sort(eigs(H, 15, 'smallestabs'));
    
    err_L_E0(idx) = abs(1 - E_num(1) / E_exa_0);
    err_L_E7(idx) = abs(1 - E_num(8) / E_exa_7);
    err_L_E13(idx) = abs(1 - E_num(14) / E_exa_13);
end

figure('Name', 'Optimizacion Dominio');
semilogy(L_vec, err_L_E0, 'Color', colores(1,:), 'DisplayName', 'E_0');
hold on;
semilogy(L_vec, err_L_E7, 'Color', colores(2,:), 'DisplayName', 'E_7');
semilogy(L_vec, err_L_E13, 'Color', colores(3,:), 'DisplayName', 'E_{13}');
grid on;
xlabel('a [u.a.]');
ylabel('Error relativo');
% title(sprintf('Competición de errores: Muro artificial vs Resolución ($N=%d$)', N_fixed));
legend('Location', 'northeast');
% saveas(gcf, 'Optimizacion_Dominio_L_v2.png');

%% =========================================================================
% SECCIÓN 3: ESPECTRO NUMÉRICO FRENTE AL ANALÍTICO
% =========================================================================
% Compara el espectro de energías computado frente al comportamiento teórico esperado
disp('Calculando espectro para los parámetros unificados (L=10, N=1000)...');
L_opt_spec = 10;
N_opt_spec = 1000;
delta_spec = L_opt_spec / (N_opt_spec + 1);
x_spec = linspace(-L_opt_spec/2 + delta_spec, L_opt_spec/2 - delta_spec, N_opt_spec)';

e_spec = ones(N_opt_spec, 1);
T_spec = spdiags([-e_spec, 2*e_spec, -e_spec], [-1, 0, 1], N_opt_spec, N_opt_spec) / (2 * delta_spec^2);
V_spec = spdiags(0.5 * x_spec.^2, 0, N_opt_spec, N_opt_spec);

num_states = N_opt_spec;
E_num_spec = sort(eig(full(T_spec + V_spec)));
n_vec = 0:(num_states-1);
E_exa_spec = n_vec + 0.5;

figure('Name', 'Espectro Optimo Ampliado');

% Zoom a los primeros estados (Ajuste perfecto)
plot(n_vec, E_exa_spec, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exacta');
hold on;
plot(n_vec, E_num_spec, 'b--', 'MarkerSize', 2, 'DisplayName', 'Numérica');
xlim([0 20]); ylim([0 25]);
grid on;
xlabel('v');
ylabel('E [u.a.]');
% title('Primeros niveles (Ajuste perfecto)');
legend('Location', 'northwest');

figure('Name', 'Espectro Optimo');
plot(n_vec, E_exa_spec, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exacta');
hold on;
plot(n_vec, E_num_spec, 'b--', 'MarkerSize', 2, 'DisplayName', 'Numérica');
grid on;
xlabel('v');
ylabel('E [u.a.]');
% title('Espectro completo (Colisión y Nyquist)');
legend('Location', 'northwest');
% sgtitle(sprintf('Espectro del Oscilador Armónico (Parámetros óptimos $L=%d, N=%d$)', L_opt_spec, N_opt_spec));
% saveas(gcf, 'Espectro_Optimo_v2.png');

%% =========================================================================
% SECCIÓN 4: EVOLUCIÓN ESPACIAL DE LOS AUTOESTADOS 
% =========================================================================
% Visualiza físicamente las funciones de onda para estados de baja y alta energía
disp('Generando nueva visualización de ondas para L=10...');
L_prop = 10;
N_prop = 500;
delta_p = L_prop / (N_prop + 1);
x_p = linspace(-L_prop/2 + delta_p, L_prop/2 - delta_p, N_prop)';

e_p = ones(N_prop, 1);
T_p = spdiags([-e_p, 2*e_p, -e_p], [-1, 0, 1], N_prop, N_prop) / (2 * delta_p^2);
V_p = spdiags(0.5 * x_p.^2, 0, N_prop, N_prop);
[V_mat_p, D_p] = eigs(T_p + V_p, 25, 'smallestabs');
[E_p, idx_p] = sort(diag(D_p));

figure('Name', 'Ondas L=10');
n_plot_list = [0, 7, 13];
titles = {'Encaje holgado ($n=0$)', 'Encaje ajustado ($n=3$)', 'Asfixia por pared ($n=15$)'};

for i = 1:3
    n_state = n_plot_list(i);
    
    hold on;
    
    psi_num = V_mat_p(:, idx_p(n_state+1));
    psi_num = 2 * psi_num / sqrt(sum(psi_num.^2 * delta_p)); 
    
    % Onda analítica usando hermiteH
    H_poly = hermiteH(n_state, x_p);
    N_const = 1 / sqrt((2^n_state) * factorial(n_state) * sqrt(pi));
    psi_exa_align = N_const * H_poly .* exp(-x_p.^2 / 2);
    
    % Onda analítica extendida
    x_wide = linspace(-8, 8, 500)';
    H_poly_wide = hermiteH(n_state, x_wide);
    psi_exa_wide = N_const * H_poly_wide .* exp(-x_wide.^2 / 2);
    psi_exa_wide = 2 * psi_exa_wide + (n_state + 0.5);
    plot(x_wide, psi_exa_wide, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Alinear NUMÉRICA a la analítica
    if sum(psi_num .* psi_exa_align) < 0
        psi_num = -psi_num; 
    end
    
    psi_num = psi_num + E_p(n_state+1);
    plot(x_p, psi_num,  '--', 'Color', colores(i,:), 'LineWidth', 2, 'DisplayName', sprintf('\\Psi_{%d}', n_state));
    
    yline(n_state + 0.5, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(E_p(n_state+1), '--', 'Color', colores(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');
end

plot(NaN, NaN, 'k-', 'LineWidth', 1.0, 'DisplayName', '\Psi_v^{exacta}');

% Potencial numérico (interior de la caja)
plot(x_p, 0.5*x_p.^2, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Muros artificiales (límites de la matriz)
plot([-L_prop/2 + delta_p, -L_prop/2 + delta_p], [0.5*(-L_prop/2 + delta_p)^2, 28], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([L_prop/2 - delta_p, L_prop/2 - delta_p], [0.5*(L_prop/2 - delta_p)^2, 28], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Potencial teórico (continuación "fantasma" analítica)
x_out_L = linspace(-8, -L_prop/2, 50);
x_out_R = linspace(L_prop/2, 8, 50);
plot(x_out_L, 0.5*x_out_L.^2, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(x_out_R, 0.5*x_out_R.^2, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlim([-7 7]);
ylim([-1 16])
% title(titles{i});
xlabel('x [u.a.]');
ylabel('E [u.a.]');
legend('Location', 'best');
% sgtitle(sprintf('Evolución de los autoestados dentro de una caja fija ($L=%d$)', L_prop));
disp('Proceso completado con éxito. Gráficas del oscilador generadas.');
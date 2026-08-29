clear;
clc;
close all;

%% --- Configuración global y Parámetros del Sistema ---
V_0 = 350; % Altura del escalón 
L = 1;     % Anchura de la caja [u.a.]

N_vals = [1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 35, 40, 50];
n_states = 6;
E_num = NaN(length(N_vals), n_states);
colores_estado = lines(n_states);

%% =========================================================================
% SECCIÓN 1: ENSAMBLAJE MATRICIAL NUMÉRICO Y BÚSQUEDA ANALÍTICA 
% =========================================================================
% Construye la matriz de perturbación analítica imponiendo las reglas de selección de paridad y obtiene el espectro FDM.
% Paralelamente, localiza las raíces analíticas exactas (E_exa) solventando las ecuaciones trascendentales.

% 1.1 Cálculo de autovalores numéricos mediante diagonalización
for n_idx = 1:length(N_vals)
    n_max = N_vals(n_idx);
    H = zeros(n_max, n_max);
    for i = 1:n_max
        for j = 1:n_max
            if i == j
                H(i,j) = ((pi*i)^2 + V_0)/2;
            elseif mod(i+j,2) ~= 0 && mod(i-j,2) ~= 0
                H(i,j) = V_0*(sin(pi*(i+j)/2)/(i+j) - sin(pi*(i-j)/2)/(i-j))/pi;
            end
        end
    end
    [phi_tmp, E_tmp] = eig(H);
    diagE = diag(E_tmp);
    limit = min(n_max, n_states);
    E_num(n_idx, 1:limit) = diagE(1:limit)';
    
    % Se guarda los autovectores para N=35
    if n_max == 35
        phi_35 = phi_tmp;
    end
end

% 1.2 Búsqueda de autovalores analíticos exactos (solución trascendental)
% (Se emplea una matriz auxiliar masiva N=100 para semillar el fzero de forma robusta)
H_seed = zeros(100, 100);
for i = 1:100
    for j = 1:100
        if i == j
            H_seed(i,j) = ((pi*i)^2 + V_0)/2;
        elseif mod(i+j,2) ~= 0 && mod(i-j,2) ~= 0
            H_seed(i,j) = V_0*(sin(pi*(i+j)/2)/(i+j) - sin(pi*(i-j)/2)/(i-j))/pi;
        end
    end
end
[~, E_tmp_seed] = eig(H_seed);
semillas = sort(diag(E_tmp_seed));

g1 = @(E) sqrt(2*(V_0-E)) .* sin(sqrt(2*E)/2) .* sinh(sqrt(2*(V_0-E))/2) + ...
          sqrt(2*E) .* cos(sqrt(2*E)/2) .* cosh(sqrt(2*(V_0-E))/2);
g2 = @(E) sqrt(2*E) .* cos(sqrt(2*E)/2) .* sin(sqrt(2*(E-V_0))/2) + ...
          sqrt(2*(E-V_0)) .* cos(sqrt(2*(E-V_0))/2) .* sin(sqrt(2*E)/2);

E_exa = zeros(1,n_states);
opts_fzero = optimset('Display', 'off');
for n = 1:n_states
    if semillas(n) < V_0
        E_exa(n) = fzero(g1, semillas(n), opts_fzero);
    else
        E_exa(n) = fzero(g2, semillas(n), opts_fzero);
    end
end

% Cálculo de errores relativos
err = abs(E_num - E_exa) ./ E_exa;

%% =========================================================================
% SECCIÓN 2: CONVERGENCIA GLOBAL DEL ERROR RELATIVO
% =========================================================================
% Ilustra el error numérico frente al tamaño de la base
fig0 = figure();
hold on;
for n = 1:n_states
    plot(N_vals, err(:,n), '-', 'LineWidth', 1.5, 'Color', colores_estado(n,:), 'DisplayName', sprintf('E_%d', n));
end
set(gca, 'YScale', 'log');
xlabel('Tamaño de la base (N)');
ylabel('Error relativo');
legend('Location', 'northeast');
grid on; box on;
xlim([0 50]);
%print(fig0, 'Convergencia_error_escalon.png', '-dpng', '-r300');

%% =========================================================================
% SECCIÓN 3: DECAIMIENTO DEL PESO ESTADÍSTICO DE LA BASE
% =========================================================================
% Justifica matemáticamente el truncamiento en N=35 mostrando que la probabilidad |c_i|^2 de armónicos superiores cae a cotas despreciables (< 10^-6).
% Se divide en dos paneles: ligados (E<V0) y libres (E>V0).

% Panel (a): Estados ligados (1 a 3)
fig_coef_1 = figure('Name', 'Coeficientes Ligados');
for k = 1:3
    n = k;
    subplot(3, 1, k);
    bar(1:35, abs(phi_35(:, n)).^2, 'FaceColor', colores_estado(n,:));
    set(gca, 'YScale', 'log');
    ylim([1e-10, 1]);
    xlim([0 36]);
    ylabel('|c_i|^2');
    grid on;
    if k == 3
        xlabel('Índice de la base (i)');
    end
end
%print(fig_coef_1, 'Coeficientes_escalon_1.png', '-dpng', '-r300');

% Panel (b): Estados libres/sobre barrera (4 a 6)
fig_coef_2 = figure('Name', 'Coeficientes Libres');
for k = 1:3
    n = k + 3;
    subplot(3, 1, k);
    bar(1:35, abs(phi_35(:, n)).^2, 'FaceColor', colores_estado(n,:));
    set(gca, 'YScale', 'log');
    ylim([1e-10, 1]);
    xlim([0 36]);
    ylabel('|c_i|^2');
    grid on;
    if k == 3
        xlabel('Índice de la base (i)');
    end
end
%print(fig_coef_2, 'Coeficientes_escalon_2.png', '-dpng', '-r300');

%% =========================================================================
% SECCIÓN 4: RECONSTRUCCIÓN ESPACIAL DE LAS FUNCIONES DE ONDA
% =========================================================================
% Superpone la combinación lineal generada por los autovectores (N=35) con las funciones definidas a tramos exactas
x = linspace(0, L, 1000);
psi_num = zeros(n_states, length(x));

for state = 1:n_states
    for i = 1:35
        basis = sqrt(2/L) * sin(i * pi * x / L);
        psi_num(state, :) = psi_num(state, :) + phi_35(i, state) * basis;
    end
end

% Construcción de las funciones de onda analíticas exactas (a tramos)
psi_exa = zeros(n_states, length(x));
for n = 1:n_states
    E_val = E_exa(n);
    k = sqrt(2*E_val);
    if E_val < V_0
        kappa = sqrt(2*(V_0-E_val));
        idx_izq = (x <= L/2);
        psi_exa(n, idx_izq) = sin(k * x(idx_izq));
        idx_der = (x > L/2);
        coef_der = sin(k * L/2) / sinh(kappa * L/2);
        psi_exa(n, idx_der) = coef_der * sinh(kappa * (L - x(idx_der)));
    else
        k_prima = sqrt(2*(E_val-V_0));
        idx_izq = (x <= L/2);
        psi_exa(n, idx_izq) = sin(k * x(idx_izq));
        idx_der = (x > L/2);
        coef_der = sin(k * L/2) / sin(k_prima * L/2);
        psi_exa(n, idx_der) = coef_der * sin(k_prima * (L - x(idx_der)));
    end
    
    % Normalizar
    norm_factor = sqrt(trapz(x, abs(psi_exa(n,:)).^2));
    psi_exa(n,:) = psi_exa(n,:) / norm_factor;
    norm_factor_num = sqrt(trapz(x, abs(psi_num(n,:)).^2));
    psi_num(n,:) = psi_num(n,:) / norm_factor_num;
    
    % Corrección de fase (alinear signos)
    if dot(psi_num(n,:), psi_exa(n,:)) < 0
        psi_num(n,:) = -psi_num(n,:);
    end
end

% Panel (a): Funciones de onda ligadas bajo la barrera (E < V0)
fig1 = figure('Name', 'Ondas Ligadas');
hold on;
for i = 1:4
    plot(x, psi_exa(i,:), '-', 'Color', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(x, psi_num(i,:), '--', 'Color', colores_estado(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('\\Psi_%d', i));
end
plot(NaN, NaN, 'k-', 'LineWidth', 1.0, 'DisplayName', '\Psi_i^{exacta}');
xlabel('x [u.a.]');
ylabel('\Psi(x)');
legend('Location', 'best');
grid on; box on;
%print(fig1, 'Funcion_Onda_escalon_ligados.png', '-dpng', '-r300');

% Panel (b): Funciones de onda libres sobre la barrera (E > V0)
fig2 = figure('Name', 'Ondas Libres');
hold on;
for i = 5:6
    plot(x, psi_exa(i,:), '-', 'Color', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(x, psi_num(i,:), '--', 'Color', colores_estado(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('\\Psi_%d', i));
end
plot(NaN, NaN, 'k-', 'LineWidth', 1.0, 'DisplayName', '\Psi_i^{exacta}');
xlabel('x [u.a.]');
ylabel('\Psi(x)');
legend('Location', 'best');
grid on; box on;
%print(fig2, 'Funcion_Onda_escalon_libres.png', '-dpng', '-r300');

%% --- IMPRESIÓN DE RESULTADOS NUMÉRICOS ---
disp(' ');
disp('--- TABLA DE CONVERGENCIA ENERGÉTICA (N = 10, 35, 50) ---');
idx_10 = find(N_vals == 10);
idx_35 = find(N_vals == 35);
idx_50 = find(N_vals == 50);

fprintf('Nivel | Exacta(u.a.) | E(N=10)   | Err(%%) | E(N=35)   | Err(%%) | E(N=50)   | Err(%%)\n');
fprintf('---------------------------------------------------------------------------------------\n');
for n = 1:n_states
    fprintf('  E%d  | %10.4f   | %9.4f | %6.2f | %9.4f | %6.2f | %9.4f | %6.2f\n', ...
        n, E_exa(n), ...
        E_num(idx_10, n), err(idx_10, n)*100, ...
        E_num(idx_35, n), err(idx_35, n)*100, ...
        E_num(idx_50, n), err(idx_50, n)*100);
end
disp('---------------------------------------------------------------------------------------');
fprintf('Proceso completado con éxito. Gráficas generadas.\n');
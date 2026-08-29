clear;
clc;

%% =========================================================================
% SECCIÓN 1: CARGA DE DATOS Y ENSAMBLAJE MATRICIAL
% =========================================================================
% Carga los elementos de matriz y ensambla la Matriz Hamiltoniana por bloques.
filename = 'HelioResult_L0_N30.mat'; % <<--- CAMBIAR AQUÍ PARA PROBAR OTROS
if ~exist(filename, 'file'), error('Archivo no encontrado.'); end
load(filename);

Ze_vis = 2.0; 
l = HelioData.l; Zn = HelioData.Zn; n_vals = HelioData.n_vals; N_n = length(n_vals);
orb_chars = ['s', 'p', 'd', 'f', 'g', 'h'];
orb_char = orb_chars(l+1);

% 1. Definición de índices de los bloques
if l == 0
    start_B = 2; % Saltamos el 1s1s porque no tiene duplicado asimétrico
    N_base = 2*N_n - 1;
    idxA = 1:N_n; 
    idxB = (N_n+1):(2*N_n-1);
    
    idx_F = 1; idx_T = 2; idx_S = 3; % Fundamental, Triplete excitado, Singlete excitado
else
    start_B = 1; % Todos los estados tienen pareja, incluso el más bajo
    N_base = 2*N_n;
    idxA = 1:N_n; 
    idxB = (N_n+1):(2*N_n);
    
    idx_F = []; % No existe estado fundamental nsns
    idx_T = 1; idx_S = 2; % El primer estado ya se desdobla en Triplete y Singlete
end
N_sub = length(start_B:N_n);

H = zeros(N_base, N_base);
H_un_cuerpo = Ze_vis^2 * HelioData.T_n - Zn * Ze_vis * HelioData.Vr_n;
T_1s = 0.5 * Ze_vis^2; V_1s = -Zn * Ze_vis * 1.0; E_1s = T_1s + V_1s;

% --- BLOQUE AA (Estados |1s, ni> con |1s, nj>) ---
H(idxA, idxA) = H_un_cuerpo + E_1s * eye(N_n) + Ze_vis * HelioData.J_n;

% --- BLOQUE BB (Estados |ni, 1s> con |nj, 1s>) ---
H_sub = H_un_cuerpo(start_B:end, start_B:end) + E_1s * eye(N_sub);
H_sub = H_sub + Ze_vis * HelioData.J_n(start_B:end, start_B:end);
H(idxB, idxB) = H_sub;

% --- BLOQUE AB (Estados |1s, ni> con |nj, 1s>) ---
K_block = Ze_vis * HelioData.K_n(:, start_B:end);
H(idxA, idxB) = K_block;
H(idxB, idxA) = K_block'; % Simetría

% --- REORDENAMIENTO FÍSICO ---
P = [];
if l == 0
    P = 1; 
    for k = 2:N_n, P = [P, k, (N_n + k - 1)]; end
else
    for k = 1:N_n, P = [P, k, (N_n + k)]; end
end

H = H(P, P);

% --- FINALIZACIÓN MATRIZ ---
fprintf('\n--- MATRIZ HAMILTONIANA REORDENADA ---\n');
disp(H(1:min(6, N_base), 1:min(6, N_base)));

[V, D] = eig(H); [Energias, idx_eig] = sort(diag(D)); V = V(:, idx_eig);
fprintf('\nPrimeras Energías (u.a.):\n');
disp(Energias(1:min(6, end)));

% Reconstrucción de la base reordenada para funciones espaciales
base_original = zeros(N_base, 4);
for k = 1:N_n, base_original(k, :) = [1, 0, n_vals(k), l]; end
if l == 0
    for k = 2:N_n, base_original(N_n + k - 1, :) = [n_vals(k), l, 1, 0]; end
else
    for k = 1:N_n, base_original(N_n + k, :) = [n_vals(k), l, 1, 0]; end
end
base_reordenada = base_original(P, :);

% Función de onda hidrogenoide
R_func = @(n, L, r, Z) sqrt((2*Z/n)^3 * factorial(n-L-1)/(2*n*factorial(n+L))) .* ...
                       exp(-Z*r/n) .* (2*Z*r/n).^L .* laguerreL(n-L-1, 2*L+1, 2*Z*r/n);
R_s = @(n, L, r) R_func(n, L, r, Ze_vis);

%% =========================================================================
% SECCIÓN 2: CÁLCULO Y CACHÉ DE DENSIDADES ESPACIALES 2D
% =========================================================================
% Las densidades bidimensionales requieren un alto coste de cómputo. Esta sección las computa en el plano (r1, r2) y las guarda (.mat)
fprintf('\nPreparando mallas para representaciones gráficas 2D...\n');
r_grid = linspace(0, 25, 200);
[R1, R2] = meshgrid(r_grid, r_grid);

% Archivo de caché de Densidades
dataFile = strrep(filename, 'Result', 'Densities_Indep_y_Corr');

if exist(dataFile, 'file')
    fprintf('Cargando densidades espaciales desde %s...\n', dataFile);
    load(dataFile);
else
    fprintf('Calculando densidades espaciales 2D (esto puede tardar un poco)...\n');
    % Estado Fundamental (solo existe para l=0)
    if l == 0
        Psi_F = zeros(size(R1));
        for i = 1:N_base
            Psi_F = Psi_F + V(i, idx_F) .* (R_s(base_reordenada(i,1), base_reordenada(i,2), R1) .* R_s(base_reordenada(i,3), base_reordenada(i,4), R2));
        end
        Dens_F = abs(Psi_F).^2 .* R1.^2 .* R2.^2;
        Psi_ind_F = R_s(1, 0, R1) .* R_s(1, 0, R2);
        Dens_indep_F = abs(Psi_ind_F).^2 .* R1.^2 .* R2.^2;
    end

    % Estado Triplete y Singlete Más Bajos del Subespacio l
    Psi_T = zeros(size(R1)); Psi_S = zeros(size(R1));
    for i = 1:N_base
        Psi_T = Psi_T + V(i, idx_T) .* (R_s(base_reordenada(i,1), base_reordenada(i,2), R1) .* R_s(base_reordenada(i,3), base_reordenada(i,4), R2));
        Psi_S = Psi_S + V(i, idx_S) .* (R_s(base_reordenada(i,1), base_reordenada(i,2), R1) .* R_s(base_reordenada(i,3), base_reordenada(i,4), R2));
    end
    Dens_T = abs(Psi_T).^2 .* R1.^2 .* R2.^2;
    Dens_S = abs(Psi_S).^2 .* R1.^2 .* R2.^2;

    n_ref = n_vals(1); if l==0, n_ref = n_vals(2); end % El n de excitación más bajo
    Psi_ind_T = (R_s(1, 0, R1) .* R_s(n_ref, l, R2) - R_s(n_ref, l, R1) .* R_s(1, 0, R2))/sqrt(2);
    Psi_ind_S = (R_s(1, 0, R1) .* R_s(n_ref, l, R2) + R_s(n_ref, l, R1) .* R_s(1, 0, R2))/sqrt(2);
    Dens_indep_T = abs(Psi_ind_T).^2 .* R1.^2 .* R2.^2;
    Dens_indep_S = abs(Psi_ind_S).^2 .* R1.^2 .* R2.^2;
    
    fprintf('Guardando densidades y funciones de onda en %s...\n', dataFile);
    if l == 0
        save(dataFile, 'Dens_F', 'Dens_indep_F', 'Dens_T', 'Dens_indep_T', 'Dens_S', 'Dens_indep_S', ...
                       'Psi_F', 'Psi_ind_F', 'Psi_T', 'Psi_ind_T', 'Psi_S', 'Psi_ind_S', 'n_ref', 'r_grid', 'R1', 'R2');
    else
        save(dataFile, 'Dens_T', 'Dens_indep_T', 'Dens_S', 'Dens_indep_S', ...
                       'Psi_T', 'Psi_ind_T', 'Psi_S', 'Psi_ind_S', 'n_ref', 'r_grid', 'R1', 'R2');
    end
end

%% =========================================================================
% SECCIÓN 2.2: DIBUJO DE GRÁFICAS DE DENSIDAD 2D 
% =========================================================================
% Representa los mapas 2D para evidenciar la formación del Hueco y Montículo de Fermi
fprintf('\nGenerando representaciones gráficas 2D...\n');
if l == 0
    figure('Name', 'Estado Fundamental 1s1s');
    contourf(R1, R2, Dens_F, 10, 'LineColor', 'none'); colormap(jet); colorbar; hold on;
    clim([0, max(Dens_F(:))]);  % Fija la barra de colores a la densidad correlacionada
    contour(R1, R2, Dens_indep_F, 10, 'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
    xlabel('r_1 [u.a.]'); ylabel('r_2 [u.a.]'); axis square;
    axis([0, 2.5, 0, 2.5])
end

figure('Name', 'Triplete Más Bajo');
contourf(R1, R2, Dens_T, 10, 'LineColor', 'none'); colormap(jet); colorbar; hold on;
clim([0, max(Dens_T(:))]);  % Fija la barra de colores a la densidad correlacionada
contour(R1, R2, Dens_indep_T, 10, 'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('r_1 [u.a.]'); ylabel('r_2 [u.a.]'); axis square;
axis([0, 20, 0, 20]);

figure('Name', 'Singlete Más Bajo');
contourf(R1, R2, Dens_S, 10, 'LineColor', 'none'); colormap(jet); colorbar; hold on;
clim([0, max(Dens_S(:))]);  % Fija la barra de colores a la densidad correlacionada
contour(R1, R2, Dens_indep_S, 10, 'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('r_1 [u.a.]'); ylabel('r_2 [u.a.]'); axis square;
axis([0, 20, 0, 20]);

%% =========================================================================
% SECCIÓN 3: ESTUDIO VARIACIONAL (Ze)
% =========================================================================
% Diagonaliza la matriz iterativamente barriendo distintos valores de Ze para demostrar el efecto de apantallamiento mutuo en el mínimo de energía
fprintf('Generando curva de optimización de Ze...\n');
Ze_tests = linspace(1.5, 2.5, 50);
E_opt_1 = zeros(size(Ze_tests)); E_opt_2 = zeros(size(Ze_tests)); E_opt_3 = zeros(size(Ze_tests));

for idx_Z = 1:length(Ze_tests)
    Z_t = Ze_tests(idx_Z);
    H_temp = zeros(N_base, N_base);
    H_1body = Z_t^2 * HelioData.T_n - Zn * Z_t * HelioData.Vr_n;
    E_1s_t = 0.5 * Z_t^2 - Zn * Z_t * 1.0;
    
    H_temp(idxA, idxA) = H_1body + E_1s_t * eye(N_n) + Z_t * HelioData.J_n;
    H_temp(idxB, idxB) = H_1body(start_B:end, start_B:end) + E_1s_t * eye(N_sub) + Z_t * HelioData.J_n(start_B:end, start_B:end);
    K_block_t = Z_t * HelioData.K_n(:, start_B:end);
    H_temp(idxA, idxB) = K_block_t; H_temp(idxB, idxA) = K_block_t';
    
    H_temp = H_temp(P, P); E_temp = sort(eig(H_temp));
    E_opt_1(idx_Z) = E_temp(1);
    if length(E_temp) >= 2, E_opt_2(idx_Z) = E_temp(2); end
    if length(E_temp) >= 3, E_opt_3(idx_Z) = E_temp(3); end
end
figure('Name', 'Optimización Variacional (Ze)');
hold on;
[min_E0, idx_min0] = min(E_opt_1); % Se guarda para el Teorema del Virial

if l == 0
    plot(Ze_tests, E_opt_1, 'LineWidth', 2, 'DisplayName', sprintf('1^1S'));
    plot(Ze_tests(idx_min0), min_E0, 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
    
    if N_base >= 2
        [min_E2, idx_min2] = min(E_opt_2);
        plot(Ze_tests, E_opt_2, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', sprintf('2^3S'));
        plot(Ze_tests(idx_min2), min_E2, 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
    end
    if N_base >= 3
        [min_E3, idx_min3] = min(E_opt_3);
        plot(Ze_tests, E_opt_3, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', sprintf('2^1S'));
        plot(Ze_tests(idx_min3), min_E3, 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
    end
else
    plot(Ze_tests, E_opt_1, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3, 'DisplayName', sprintf('3^3D'));
    plot(Ze_tests(idx_min0), min_E0, 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
    
    if N_base >= 2
        [min_E2, idx_min2] = min(E_opt_2);
        plot(Ze_tests, E_opt_2, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', sprintf('3^1D'));
        plot(Ze_tests(idx_min2), min_E2, 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
    end
end

xlabel('Z_{eff}'); ylabel('E [u.a.]');
legend('Location', 'northeastoutside'); grid on;

%% =========================================================================
% SECCIÓN 4: CONVERGENCIA DE ENERGÍAS 
% =========================================================================
% Analiza cómo se estabilizan los autovalores a medida que aumenta el tamaño de la base (n_max), justificando matemáticamente su truncamiento
fprintf('Generando curva de convergencia...\n');
n_max_evals = min(n_vals):max(n_vals);
E_conv = NaN(length(n_max_evals), 3);
for idx_n = 1:length(n_max_evals)
    n_corte = n_max_evals(idx_n);
    valid_idx = find(base_reordenada(:,1) <= n_corte & base_reordenada(:,3) <= n_corte);
    if ~isempty(valid_idx)
        H_sub = H(valid_idx, valid_idx);
        E_sub = sort(eig(H_sub));
        E_conv(idx_n, 1) = E_sub(1);
        if length(E_sub) >= 2, E_conv(idx_n, 2) = E_sub(2); end
        if length(E_sub) >= 3, E_conv(idx_n, 3) = E_sub(3); end
    end
end

figure('Name', 'Convergencia de Energías');

if l == 0
    subplot(3,1,1)
    plot(n_max_evals(3:end), E_conv(3:end,1), 'LineWidth', 2, 'DisplayName', '1^1S'); hold on;
    xlabel('n_{max}'); ylabel('E [u.a.]');
    legend('Location', 'northeast'); grid on;
    subplot(3,1,2)
    plot(n_max_evals(3:end), E_conv(3:end,2), 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', '2^3S');
    xlabel('n_{max}'); ylabel('E [u.a.]');
    legend('Location', 'northeast'); grid on;
    subplot(3,1,3)
    plot(n_max_evals(3:end), E_conv(3:end,3), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '2^1S');
    xlabel('n_{max}'); ylabel('E [u.a.]');
    legend('Location', 'northeast'); grid on;
else
    subplot(2,1,1)
    plot(n_max_evals(5:end), E_conv(5:end,1), 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', sprintf('3^3D')); hold on;
    xlabel('n_{max}'); ylabel('E [u.a.]');
    legend('Location', 'northeast'); grid on;
    subplot(2,1,2)
    plot(n_max_evals(5:end), E_conv(5:end,2), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', sprintf('3^1D'));
    xlabel('n_{max}'); ylabel('E [u.a.]');
    legend('Location', 'northeast'); grid on;
end

%% =========================================================================
% SECCIÓN 5: ANÁLISIS DE PESOS (COEFICIENTES) 
% =========================================================================
% Grafica los coeficientes de los autovectores.
fprintf('Generando gráfico de pesos...\n');
figure('Name', 'Pesos de la Base');

if l == 0
    pesos_1 = abs(V(:, idx_F)).^2;
    pesos_2 = abs(V(:, idx_T)).^2;
    pesos_3 = abs(V(:, idx_S)).^2;

    subplot(3,1,1);
    bar(pesos_1);
    ylabel('|c_{ i}|^2'); set(gca, 'YScale', 'log'); grid on;
    
    subplot(3,1,2);
    bar(pesos_2);
    ylabel('|c_{ i}|^2'); set(gca, 'YScale', 'log'); grid on;
    
    subplot(3,1,3);
    bar(pesos_3);
    xlabel('Índice de la base'); ylabel('|c_{ i}|^2');
    set(gca, 'YScale', 'log'); grid on;
    
else
    pesos_1 = abs(V(:, idx_T)).^2;
    pesos_2 = abs(V(:, idx_S)).^2;
    
    subplot(2,1,1);
    bar(pesos_1);
    ylabel('|c_{ i}|^2'); set(gca, 'YScale', 'log');
    xlabel('Índice de la base'); ylabel('|c_{ i}|^2'); grid on;
    
    subplot(2,1,2);
    bar(pesos_2);
    xlabel('Índice de la base'); ylabel('|c_{ i}|^2');
    set(gca, 'YScale', 'log'); grid on;
end

%% =========================================================================
% SECCIÓN 6: DENSIDAD CONDICIONAL 
% =========================================================================
% Realiza cortes 1D fijando la posición del electrón interno (r1)
fprintf('Generando Densidades Condicionales para varios r1...\n');
r_fijos = [0.5, 1.0, 1.5]; % Valores de r1 fijos
r_var = linspace(0, 25, 200);

colores_param = [0.1882, 0.2902, 0.5843; % Azul oscuro (#304a95)
                 0.8824, 0.1020, 0.1608; % Rojo intenso (#e11a29)
                 0.3961, 0.6863, 0.1882];% Verde vivo (#65af30)

grosor = [4, 3, 2];

if l == 0
    % Array de índices para Fundamental, Triplete y Singlete
    indices_estados = [idx_F, idx_T, idx_S];
    nombres_estados = {'Fundamental 1^1S', 'Triplete 2^3S', 'Singlete 2^1S'};
    
    for k = 1:length(indices_estados)
        idx_estado = indices_estados(k);
        figure('Name', ['Densidad Condicional: ' nombres_estados{k}]);
        hold on;
        
        for j = 1:length(r_fijos)
            r_fijo = r_fijos(j);
            c = colores_param(j,:);
            g = grosor(j);
            
            Psi_cond = zeros(size(r_var));
            for i = 1:N_base
                Psi_cond = Psi_cond + V(i, idx_estado) .* (R_s(base_reordenada(i,1), base_reordenada(i,2), r_fijo) .* R_s(base_reordenada(i,3), base_reordenada(i,4), r_var));
            end
            
            % El MPI para l=0 es diferente si es fundamental o excitado
            if k == 1
                Psi_indep_cond = R_s(1, 0, r_fijo) .* R_s(1, 0, r_var);
            elseif k == 2 % Triplete MPI (antisimétrico espacial)
                Psi_indep_cond = (R_s(1, 0, r_fijo) .* R_s(n_ref, 0, r_var) - R_s(n_ref, 0, r_fijo) .* R_s(1, 0, r_var))/sqrt(2);
            elseif k == 3 % Singlete MPI (simétrico espacial)
                Psi_indep_cond = (R_s(1, 0, r_fijo) .* R_s(n_ref, 0, r_var) + R_s(n_ref, 0, r_fijo) .* R_s(1, 0, r_var))/sqrt(2);
            end
            
            Dens_cond = abs(Psi_cond).^2 .* r_var.^2; 
            Dens_indep_cond = abs(Psi_indep_cond).^2 .* r_var.^2;
            
            % Normalizamos el área
            Dens_cond = Dens_cond / trapz(r_var, Dens_cond);
            Dens_indep_cond = Dens_indep_cond / trapz(r_var, Dens_indep_cond);
            
            if k == 1
                if j == 1
                    plot(r_var, Dens_indep_cond, 'k--', 'LineWidth', 2.0, 'DisplayName', 'MPI');
                end
            else
                plot(r_var, Dens_indep_cond, 'color', c, 'linestyle', '--', 'LineWidth', 1.5, 'DisplayName', sprintf('MPI (r_1 = %.1f)', r_fijo));
            end
            
            plot(r_var, Dens_cond, 'color', c, 'LineWidth', g, 'DisplayName', sprintf('CI (r_1 = %.1f)', r_fijo));
        end
        
        xlabel('r_2 [u.a.]'); 
        ylabel('\rho(r_2|r_1)');
        xlim([0, 25])
        legend('Location', 'northeast'); grid on;
    end
else
    figure('Name', 'Densidad Condicional');
    hold on;
    for j = 1:length(r_fijos)
        r_fijo = r_fijos(j);
        c = colores_param(j,:);
        g = grosor(j);
        
        Psi_cond = zeros(size(r_var));
        for i = 1:N_base
            Psi_cond = Psi_cond + V(i, idx_T) .* (R_s(base_reordenada(i,1), base_reordenada(i,2), r_fijo) .* R_s(base_reordenada(i,3), base_reordenada(i,4), r_var));
        end
        Psi_indep_cond = (R_s(1, 0, r_fijo) .* R_s(n_ref, l, r_var) - R_s(n_ref, l, r_fijo) .* R_s(1, 0, r_var))/sqrt(2);
        
        Dens_cond = abs(Psi_cond).^2 .* r_var.^2; 
        Dens_indep_cond = abs(Psi_indep_cond).^2 .* r_var.^2;
        
        % Normalizamos
        Dens_cond = Dens_cond / trapz(r_var, Dens_cond);
        Dens_indep_cond = Dens_indep_cond / trapz(r_var, Dens_indep_cond);
        
        plot(r_var, Dens_indep_cond, 'color', c, 'linestyle', '--', 'LineWidth', 1.5, 'DisplayName', sprintf('MPI (r_1 = %.1f)', r_fijo));
        plot(r_var, Dens_cond, 'color', c, 'LineWidth', g, 'DisplayName', sprintf('Corr. (r_1 = %.1f)', r_fijo));
    end
    xlabel('r [u.a.]'); 
    ylabel('\rho(r)');
    xlim([0, 25])
    legend('Location', 'northeast'); grid on;
end

%% =========================================================================
% SECCIÓN 7: AMPLITUD DE LOS ESTADOS
% =========================================================================
% Proyecta la función total sobre el 1s interno para obtener la amplitud radial del electrón externo
fprintf('\nCalculando la amplitud radial de los estados...\n');

if l == 0
    % -----------------------------------------------------------------
    % 1. Amplitud del Estado Fundamental 1s1s (integrando directo)
    % -----------------------------------------------------------------
    Amplitud_F_Corr = trapz(r_grid, Psi_F .* R1.^2, 2); 
    Amplitud_F_Indep = trapz(r_grid, Psi_ind_F .* R1.^2, 2);
    
    Amplitud_F_Corr = Amplitud_F_Corr / sqrt(trapz(r_grid, abs(Amplitud_F_Corr).^2 .* r_grid(:).^2));
    Amplitud_F_Indep = Amplitud_F_Indep / sqrt(trapz(r_grid, abs(Amplitud_F_Indep).^2 .* r_grid(:).^2));

    figure('Name', 'Amplitud 1D del Estado Fundamental');
    plot(r_grid, Amplitud_F_Corr, 'color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'DisplayName', '1^1S'); hold on;
    plot(r_grid, Amplitud_F_Indep, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');
    
    xlabel('r [u.a.]'); ylabel('R(r)');
    xlim([0, 10])
    legend('Location', 'northeast'); grid on;

    % -----------------------------------------------------------------
    % 2. Amplitud de los Estados Excitados L=0 (Triplete y Singlete 1s2s)
    % -----------------------------------------------------------------
    % Usamos la PROYECCIÓN (el "bisturí") contra el orbital interno 1s(r1)
    Phi_1s = R_s(1, 0, R1); 
    
    Amplitud_T_out = trapz(r_grid, Psi_T .* Phi_1s .* R1.^2, 2);
    Amplitud_S_out = trapz(r_grid, Psi_S .* Phi_1s .* R1.^2, 2);
    
    % Modelo independiente (El orbital 2s puro)
    Amplitud_indep_out = R_s(2, 0, r_grid)'; 
    
    % Normalización cuántica rigurosa (Integral de |R|^2 * r^2 dr = 1)
    Amplitud_T_out = Amplitud_T_out / sqrt(trapz(r_grid, abs(Amplitud_T_out).^2 .* r_grid(:).^2));
    Amplitud_S_out = Amplitud_S_out / sqrt(trapz(r_grid, abs(Amplitud_S_out).^2 .* r_grid(:).^2));
    Amplitud_indep_out = Amplitud_indep_out / sqrt(trapz(r_grid, abs(Amplitud_indep_out).^2 .* r_grid(:).^2));
    
    % Forzamos a que el signo global sea positivo en el origen para facilitar la comparación visual
    if Amplitud_T_out(1) < 0; Amplitud_T_out = -Amplitud_T_out; end
    if Amplitud_S_out(1) < 0; Amplitud_S_out = -Amplitud_S_out; end
    if Amplitud_indep_out(1) < 0; Amplitud_indep_out = -Amplitud_indep_out; end

    figure('Name', 'Amplitud Radial de los Estados Excitados (2s)');
    plot(r_grid, Amplitud_T_out, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', '2^3S'); hold on;
    plot(r_grid, Amplitud_S_out, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '2^1S');
    plot(r_grid, Amplitud_indep_out, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');

    yline(0, 'r-', 'HandleVisibility', 'off');    
    xlabel('r [u.a.]'); ylabel('R(r)');
    xlim([0, 15])
    legend('Location', 'northeast'); grid on;

else
    % En L>0 extraemos el orbital exterior filtrando las configuraciones (1s, nl)
    Amplitud_T_out = zeros(size(r_grid));
    Amplitud_S_out = zeros(size(r_grid));
    Amplitud_indep_out = R_s(n_ref, l, r_grid)'; 
    
    for i = 1:N_base
        % Filtramos: solo sumamos si el electrón 2 es el excitado (tiene momento l)
        if base_reordenada(i,4) == l
            n_exterior = base_reordenada(i,3); 
            Amplitud_T_out = Amplitud_T_out + V(i, idx_T) * R_s(n_exterior, l, r_grid)';
            Amplitud_S_out = Amplitud_S_out + V(i, idx_S) * R_s(n_exterior, l, r_grid)';
        end
    end
    
    % Normalización cuántica rigurosa (Integral de |R|^2 * r^2 dr = 1)
    Amplitud_T_out = Amplitud_T_out / sqrt(trapz(r_grid, abs(Amplitud_T_out).^2 .* r_grid(:).^2));
    Amplitud_S_out = Amplitud_S_out / sqrt(trapz(r_grid, abs(Amplitud_S_out).^2 .* r_grid(:).^2));
    Amplitud_indep_out = Amplitud_indep_out / sqrt(trapz(r_grid, abs(Amplitud_indep_out).^2 .* r_grid(:).^2));

    
    figure('Name', sprintf('Amplitud Radial (Orbital %c)', orb_char));
    plot(r_grid, Amplitud_T_out, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', '2^3P'); hold on;
    plot(r_grid, Amplitud_S_out, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '2^1P');
    plot(r_grid, Amplitud_indep_out, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');
    
    xlabel('r [u.a.]'); 
    ylabel('R(r)');
    xlim([0, 25])
    legend('Location', 'northeast'); grid on;
end

%% =========================================================================
% SECCIÓN 8: DENSIDAD DE PROBABILIDAD MARGINAL
% =========================================================================
% Integra la probabilidad conjunta sobre una de las coordenadas para obtener el perfil de densidad radial global absoluto de los electrones
fprintf('\nCalculando la Densidad de Probabilidad Marginal...\n');

if l == 0
    r_grid_plot = r_grid(:)'; % Aseguramos vector fila
    
    % -----------------------------------------------------------------
    % 1. Densidad Marginal del Estado Fundamental 1s1s
    % -----------------------------------------------------------------
    Dens_marg_F_Corr = trapz(r_grid, Dens_F, 1);
    Dens_marg_F_Indep = trapz(r_grid, Dens_indep_F, 1);
    
    Dens_marg_F_Corr = Dens_marg_F_Corr(:)';
    Dens_marg_F_Indep = Dens_marg_F_Indep(:)';
    
    % Normalización rigurosa
    Dens_marg_F_Corr = Dens_marg_F_Corr / trapz(r_grid_plot, Dens_marg_F_Corr);
    Dens_marg_F_Indep = Dens_marg_F_Indep / trapz(r_grid_plot, Dens_marg_F_Indep);
    
    figure('Name', 'Densidad Marginal - Estado Fundamental');
    plot(r_grid_plot, Dens_marg_F_Corr, 'color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'DisplayName', '1^1S'); hold on;
    plot(r_grid_plot, Dens_marg_F_Indep, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');
    
    xlabel('r [u.a.]'); ylabel('\rho(r)');
    xlim([0, 15]);
    legend('Location', 'northeast'); grid on;

    % -----------------------------------------------------------------
    % 2. Densidad Marginal de los Estados Excitados L=0 (1s2s)
    % -----------------------------------------------------------------
    Dens_marg_T = trapz(r_grid, Dens_T, 1);
    Dens_marg_S = trapz(r_grid, Dens_S, 1);
    Dens_marg_indep = trapz(r_grid, Dens_indep_T, 1); 
    
    Dens_marg_T = Dens_marg_T(:)';
    Dens_marg_S = Dens_marg_S(:)';
    Dens_marg_indep = Dens_marg_indep(:)';
    
    % Normalización rigurosa
    Dens_marg_T = Dens_marg_T / trapz(r_grid_plot, Dens_marg_T);
    Dens_marg_S = Dens_marg_S / trapz(r_grid_plot, Dens_marg_S);
    Dens_marg_indep = Dens_marg_indep / trapz(r_grid_plot, Dens_marg_indep);
    
    figure('Name', 'Densidad Marginal - Estados Excitados (2s)');
    plot(r_grid_plot, Dens_marg_T, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', '2^3S'); hold on;
    plot(r_grid_plot, Dens_marg_S, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '2^1S');
    plot(r_grid_plot, Dens_marg_indep, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');
    
    xlabel('r [u.a.]'); ylabel('\rho(r)');
    xlim([0, 15]);
    legend('Location', 'northeast'); grid on;
    
else
    % -----------------------------------------------------------------
    % Densidad Marginal de los Estados L > 0
    % -----------------------------------------------------------------
    r_grid_plot = r_grid(:)';
    
    Dens_marg_T = trapz(r_grid, Dens_T, 1);
    Dens_marg_S = trapz(r_grid, Dens_S, 1);
    Dens_marg_indep = trapz(r_grid, Dens_indep_T, 1); 
    
    Dens_marg_T = Dens_marg_T(:)';
    Dens_marg_S = Dens_marg_S(:)';
    Dens_marg_indep = Dens_marg_indep(:)';
    
    Dens_marg_T = Dens_marg_T / trapz(r_grid_plot, Dens_marg_T);
    Dens_marg_S = Dens_marg_S / trapz(r_grid_plot, Dens_marg_S);
    Dens_marg_indep = Dens_marg_indep / trapz(r_grid_plot, Dens_marg_indep);
    
    figure('Name', sprintf('Densidad Marginal - Estados Excitados (%c)', orb_char));
    plot(r_grid_plot, Dens_marg_T, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'DisplayName', '2^3P'); hold on;
    plot(r_grid_plot, Dens_marg_S, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '2^1P');
    plot(r_grid_plot, Dens_marg_indep, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MPI');
    
    xlabel('r [u.a.]'); ylabel('\rho(r)');
    xlim([0, 25]);
    legend('Location', 'northeast'); grid on;
end

%% =========================================================================
% SECCIÓN 9: TEOREMA DEL VIRIAL (<V> / <T>) 
% =========================================================================
fprintf('\n--- TEOREMA DEL VIRIAL ---\n');
Z_opt = Ze_tests(idx_min0); 
T_mat = zeros(N_base, N_base); V_mat = zeros(N_base, N_base);
T_1body = Z_opt^2 * HelioData.T_n; V_1body = -Zn * Z_opt * HelioData.Vr_n;

T_mat(idxA, idxA) = T_1body + (0.5 * Z_opt^2) * eye(N_n);
T_mat(idxB, idxB) = T_1body(start_B:end, start_B:end) + (0.5 * Z_opt^2) * eye(N_sub);

V_mat(idxA, idxA) = V_1body - (Zn * Z_opt * 1.0) * eye(N_n) + Z_opt * HelioData.J_n;
V_mat(idxB, idxB) = V_1body(start_B:end, start_B:end) - (Zn * Z_opt * 1.0) * eye(N_sub) + Z_opt * HelioData.J_n(start_B:end, start_B:end);
V_mat(idxA, idxB) = Z_opt * HelioData.K_n(:, start_B:end); V_mat(idxB, idxA) = Z_opt * HelioData.K_n(:, start_B:end)';

T_mat = T_mat(P, P); V_mat = V_mat(P, P); H_opt = T_mat + V_mat;
[V_opt, D_opt] = eig(H_opt); [~, idx_opt] = sort(diag(D_opt)); V_opt = V_opt(:, idx_opt);

% Virial para el estado más bajo disponible
Psi_opt_ref = V_opt(:, 1);
T_exp = Psi_opt_ref' * T_mat * Psi_opt_ref;
V_exp = Psi_opt_ref' * V_mat * Psi_opt_ref;
fprintf('  -> <T> = %.4f u.a.\n  -> <V> = %.4f u.a.\n  -> Ratio <V>/<T> = %.5f\n', T_exp, V_exp, V_exp/T_exp);

%% =========================================================================
% --- SECCIÓN 10: DESDOBLAMIENTO Y RESIDUOS (Fig. 5.4) ---
% =========================================================================
% Mide la brecha energética Singlete-Triplete originada por la integral de intercambio y la contrasta contra los valores empíricos tabulados por el NIST.
fprintf('\nGenerando curva de Desdoblamiento de Intercambio y Residuos...\n');
if l == 0
    n_levels = n_vals(2:end); 
    idx_T_func = @(k) 2 * k;       % Tripletes
    idx_S_func = @(k) 2 * k + 1;   % Singletes
else
    n_levels = n_vals;
    idx_T_func = @(k) 2 * k - 1;   
    idx_S_func = @(k) 2 * k;       
end

Delta_E = zeros(size(n_levels));
for k = 1:length(n_levels)
    iT = idx_T_func(k); iS = idx_S_func(k);
    if iS <= length(Energias)
        Delta_E(k) = Energias(iS) - Energias(iT);
    end
end

% Datos experimentales (NIST Atomic Spectra Database)
Delta_E_exp = [];
if l == 0
    n_exp = 2:10;
    Delta_E_exp = [0.02925835, 0.00741788, 0.00292568, 0.00144217, ...
                   0.00081444, 0.00050419, 0.00033353, 0.00023198, 0.00016782];
elseif l == 1
    n_exp = 2:10;
    Delta_E_exp = [0.00933197, 0.00293794, 0.00125605, 0.00064588, ...
                   0.00037437, 0.00023589, 0.00015805, 0.00011099, 0.00008091];
elseif l == 2
    n_exp = 3:10;
    Delta_E_exp = [0.00001559, 0.00000900, 0.00000519, 0.00000319, ...
                   0.00000208, 0.00000142, 0.00000101, 0.00000075]; 
end

if any(Delta_E > 0)
    figure('Name', 'Desdoblamiento y Residuos');
    
    % 1. Dibujamos el Teórico
    plot(n_levels(1:9), Delta_E(1:9), 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Teórico');
    hold on;
    
    % 2. Dibujamos el Experimental y calculamos el Residuo
    if ~isempty(n_exp)
        idx_valid_num = ismember(n_levels, n_exp);
        idx_valid_exp = ismember(n_exp, n_levels);
        
        % Experimental
        plot(n_exp(idx_valid_exp), Delta_E_exp(idx_valid_exp), 'rs--', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental');
        
        % Residuo
        if any(idx_valid_num)
            n_res = n_levels(idx_valid_num);
            vec_num = Delta_E(idx_valid_num);
            vec_exp = Delta_E_exp(idx_valid_exp);
            Error_Absoluto = abs(vec_num(:) - vec_exp(:));
            
            % Pintamos el residuo en negro con triángulos punteados
            plot(n_res, Error_Absoluto, 'k^:', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Residuo');
        end
    end
    
    xlabel(sprintf('Estado n%s', orb_char)); 
    ylabel('\DeltaE [u.a.]');
    legend('Location', 'northeast');
    set(gca, 'YScale', 'log'); grid on;
end
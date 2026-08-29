clear; 
clc;

% --- Configuración del cálculo ---
l_list = 2; % Número cuántico del momento angular orbital (l)      
n_max_list = [5,8,12,15,20,25,30]; % Límite máximo del número cuántico principal para el truncamiento de la base
Zn = 2.0; % Carga nuclear efectiva            

% --- Activación del pool de workers para el cálculo paralelo de integrales ---
if isempty(gcp('nocreate'))
    fprintf('Iniciando 8 motores paralelos... ');
    parpool(8); 
end

% --- Bucle con el cálculo de los elementos de matriz ---
for l = l_list
    for n_max = n_max_list
        filename = sprintf('HelioResult_L%d_N%d.mat', l, n_max); % Nombre del archivo donde se guardan los resultados
        
        % Comprueba si el archivo de resultados actual ya existe para saltar su cálculo y ahorrar tiempo.
        if exist(filename, 'file')
            fprintf('>> Saltando L=%d, N=%d (Ya existe)\n', l, n_max);
            continue; 
        end
        
        fprintf('\n== Cálculo Incremental: L=%d, Nmax=%d ==\n', l, n_max);
        tic;
        
        n_vals = (l+1):n_max;
        N_n = length(n_vals);
        Rmax = 3 * n_max^2; % Límite superior de la integral dinámico para evitar cortes y NaN
        R = @(n, l, r, Z) sqrt((2*Z/n)^3 * factorial(n-l-1)/(2*n*factorial(n+l))) .* ...
                      exp(-Z*r/n) .* (2*Z*r/n).^l .* laguerreL(n-l-1, 2*l+1, 2*Z*r/n); % Función radial del átomo de hidrógeno

        T_n = zeros(N_n, N_n); Vr_n = zeros(N_n, N_n);
        J_n = zeros(N_n, N_n); K_n = zeros(N_n, N_n);

        % Búsqueda del archivo previo con la matriz de mayor dimensión disponible para realizar un cálculo incremental reaprovechamiento de datos ya integrados.
        d = dir(sprintf('HelioResult_L%d_N*.mat', l));
        n_old = 0; J_prev = []; K_prev = []; T_prev = []; Vr_prev = [];
        
        if ~isempty(d)
            n_existentes = cellfun(@(s) sscanf(s, sprintf('HelioResult_L%d_N%%d.mat',l)), {d.name});
            n_posibles = n_existentes(n_existentes < n_max);
            if ~isempty(n_posibles)
                best_n = max(n_posibles);
                prev_file = sprintf('HelioResult_L%d_N%d.mat', l, best_n);
                data_struct = load(prev_file); HData = data_struct.HelioData;
                n_old = size(HData.J_n, 1); 
                J_prev = HData.J_n; K_prev = HData.K_n;
                T_prev = HData.T_n; Vr_prev = HData.Vr_n;
                fprintf('  Reutilizando datos de N=%d (%d filas recuperadas)\n', best_n, n_old);
            end
        end

        % Bucle para calcular los elementos de matriz del operador 1/r, necesarios para evaluar el desajuste de carga nuclear, y las energías base.
        for i = 1:N_n
            for j = i:N_n
                if i <= n_old && j <= n_old
                    T_n(i,j) = T_prev(i,j); Vr_n(i,j) = Vr_prev(i,j);
                else
                    if n_vals(i) == n_vals(j)
                        Vr_n(i,j) = 1/n_vals(i)^2; T_n(i,j) = 1/(2*n_vals(i)^2);
                    else
                        Vr_n(i,j) = integral(@(r) R(n_vals(i),l,r,1).*R(n_vals(j),l,r,1).*r, 0, Rmax);
                        T_n(i,j) = Vr_n(i,j);
                    end
                end
                T_n(j,i) = T_n(i,j); Vr_n(j,i) = Vr_n(i,j);
            end
        end

        % Bucle para calcular los elementos de repulsión interelectrónica para la submatriz triangular superior
        calc_V_k = @(nA, lA, nC, lC, nB, lB, nD, lD, k) ...
            integral2(@(r1, r2) (R(nA,lA,r1,1).*R(nC,lC,r1,1).*r1.^2) .* ...
                                (R(nB,lB,r2,1).*R(nD,lD,r2,1).*r2.^2) .* ...
                                (min(r1,r2).^k ./ max(r1,r2).^(k+1)), ...
                       0, Rmax, 0, Rmax, 'RelTol', 1e-6); % Función que calcula el elemento de matriz de repulsión entre los estados |A,B>* y |C,D>

        J_pad = zeros(N_n, N_n); K_pad = zeros(N_n, N_n);
        if n_old > 0, J_pad(1:n_old, 1:n_old) = J_prev; K_pad(1:n_old, 1:n_old) = K_prev; end
        local_n_old = n_old;

        parfor i = 1:N_n
            ni = n_vals(i);
            temp_row_J = J_pad(i, :); temp_row_K = K_pad(i, :);
            j_start = i; if i <= local_n_old, j_start = local_n_old + 1; end
            for j = j_start:N_n
                nj = n_vals(j);
                % J: Interacción directa de Coulomb (multipolo monopolar, k=0)
                temp_row_J(j) = calc_V_k(1, 0, 1, 0, ni, l, nj, l, 0);
                % K: Interacción de intercambio (multipolo k=l dictado por simetría angular)
                temp_row_K(j) = (1/(2*l+1)) * calc_V_k(1, 0, nj, l, ni, l, 1, 0, l);
            end
            J_n(i,:) = temp_row_J; K_n(i,:) = temp_row_K;
        end
        % Simetrización explícita de las matrices (H_ij = H_ji) para reducir el tiempo de integración a la mitad y evitar asimetrías por ruido numérico.
        J_n = triu(J_n) + triu(J_n, 1)'; 
        K_n = triu(K_n) + triu(K_n, 1)';

        HelioData.T_n = T_n; HelioData.Vr_n = Vr_n;
        HelioData.J_n = J_n; HelioData.K_n = K_n;
        HelioData.l = l; HelioData.n_max = n_max; HelioData.Zn = Zn; HelioData.n_vals = n_vals;

        save(filename, 'HelioData');
        fprintf('>> %s generado en %.2f segundos.\n', filename, toc);
    end
end
fprintf('\nPROCESO COMPLETADO.\n');
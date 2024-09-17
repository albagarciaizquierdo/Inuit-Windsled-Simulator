function figures = fun_draw_results_eq(x, t, p, frame)
%% Description: 
% This function draws the results of state vector x vs time specifically
% for the equilibrium state conditions
%% Inputs: 
% x --> state vector 
% t --> timespan 
% p --> struct containing parameters
% frame --> reference frame to show kite's results 
%   frame = 1 --> SK
%   frame = 0 --> SE
%% Outputs: 
% figures 

%% Plot settings
color.blue = [0 0 1];
color.red = [1 0 0];
color.green = [0 0.6 0];
color.orange = [1 0.7 0.1];
color.purple = [0.5 0.2 0.6];

set(groot,'defaultAxesFontSize', 11);
set(groot,'defaultLineLineWidth', 1);
set(groot,'defaultAxesColorOrder', [color.blue; color.red; color.green; color.orange; color.purple]);
set(groot,'defaultAxesTickLabelInterpreter', 'latex');
set(groot,'defaultAxesGridAlpha', 0.3);
set(groot,'defaultAxesXMinorTick', 'on');
set(groot,'defaultAxesYMinorTick', 'on');
set(groot,'defaultFigureRenderer', 'painters');
set(groot,'defaultLegendBox', 'off');
set(groot,'defaultLegendInterpreter', 'latex');
set(groot,'defaultLegendLocation', 'best');
set(groot,'defaultLineMarkerSize', 1);
set(groot,'defaultTextInterpreter', 'latex');

    %% Helper function to adjust y-limits
    function adjust_ylim(ax, data)
        if max(abs(data(:))) < 10^-10
            ylim(ax, [-1, 1]);
        end
    end

    %% SYSTEM OVERVIEW IN EQUILIBRIUM
    figure;
    subplot(2,2,1)
    plot(t, x(:,1),':', t, x(:,2),'--')
    title('Position of the sled in $$S_E$$')
    xlabel('$t$ (s)')
    ylabel('Position (m)')
    legend('$x_S$','$y_S$')
    adjust_ylim(gca, x(:, 1:2));
    
    subplot(2,2,2)
    plot(t, x(:,4),':', t, x(:,5),'--')
    title('Sled velocity in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Velocity (m/s)')
    legend('$v_x^S$','$v_y^S$')
    adjust_ylim(gca, x(:, 4:5));
    
    subplot(2,2,3)
    plot(t, x(:,7),':', t, x(:,8),'--', t, x(:,9),'-.')
    title('Position of the kite in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Position (m)')
    legend({'$x_K$','$y_K$','$z_K$'},'location','southeast')
    ylim([-50 5])

    subplot(2,2,4)
    plot(t, x(:,13),':', t, x(:,14),'--', t, x(:,15),':')
    title('Kite velocity in $S_K$')
    xlabel('$t$ (s)')
    ylabel('Velocity (m/s)')
    legend({'$u$','$v$','$w$'},'location','northeast')
    adjust_ylim(gca, rad2deg(x(:, 13:15)));

%     sgtitle('System overview')
    print('mi_grafico', '-djpeg', '-r300');


    %% RESULTS
    [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
        OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m] = fun_get_results(x,p,frame);

    %% Sled forces
    figure;
    subplot(2,2,1);
    plot(t, F_S(:,1),':', t, F_S(:,2),'--', t, F_S(:,3),':')
    title('Sled force components in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Sled force (N)')
    legend({'$\textbf{F}_S \cdot\textbf{i}_E$','$\textbf{F}_S\cdot\textbf{j}_E$','$\textbf{F}_S\cdot\textbf{k}_E$'},'location','northeast')
    adjust_ylim(gca, F_S);

    subplot(2,2,2);
    plot(t, N(:,1),':', t, W_S(:,1),'--', t, F_r(:,1),'-.', t, T_ASp(:,1), t, T_ASm(:,1),'--')
    title('Sled forces in $\textbf{i}_E$')
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$'},'location','northeast')
    lgd = legend;
    lgd.NumColumns = 2;
    ylim([-50 95])

    subplot(2,2,3);
    plot(t, N(:,2),':', t, W_S(:,2),'--', t, F_r(:,2),':', t, T_ASp(:,2), t, T_ASm(:,2),'-.')
    title('Sled forces in $\textbf{j}_E$')
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$'},'Location','northeast')
    lgd = legend;
    lgd.NumColumns = 2;
    ylim([-4.5 4.5])
    
    subplot(2,2,4);
    plot(t, N(:,3),'-.', t, W_S(:,3), t, F_r(:,3),':', t, T_ASp(:,3),'--', t, T_ASm(:,3),':')
    title('Sled forces in $\textbf{k}_E$')
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$'},'location','northeast')
    lgd = legend;
    lgd.NumColumns = 2;
    ylim([-5050 5050])

%     sgtitle('Sled forces overview')

    %% Kite forces
    figure;
    subplot(2,2,1);
    plot(t, F_K(:,1),':', t, F_K(:,2),'--', t, F_K(:,3),':')
    title('Kite force components in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Sled force (N)')
    if frame == 0
        title('Kite force components in $S_E$')
        legend({'$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$'},'location','northeast')
    else
        title('Kite force components in $S_K$')
        legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
    end
    adjust_ylim(gca, F_K);

    subplot(2,2,2);
    plot(t, F_a(:,1),'-.', t, W_K(:,1), t, T_AKp(:,1),':', t, T_AKm(:,1),'--')
    if frame == 0
        title('Kite forces in $\textbf{i}_E$')
    else
        title('Kite forces in $\textbf{i}_K$')
    end
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$'},'location','east')
    lgd = legend;
    lgd.NumColumns = 2;
    adjust_ylim(gca, [F_a(:,1), W_K(:,1), T_AKp(:,1), T_AKm(:,1)]);

    subplot(2,2,3);
    plot(t, F_a(:,2),':', t, W_K(:,2),'--', t, T_AKp(:,2), t, T_AKm(:,2),'-.')
    if frame == 0
        title('Kite forces in $\textbf{j}_E$')
    else
        title('Kite forces in $\textbf{j}_K$')
    end
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$'},'Location','northeast')
    lgd = legend;
    lgd.NumColumns = 2;
    ylim([-4.5 4.5])

    subplot(2,2,4);
    plot(t, F_a(:,3),'-.', t, W_K(:,3), t, T_AKp(:,3),':', t, T_AKm(:,3),'--')
    if frame == 0
        title('Kite forces in $\textbf{k}_E$')
    else
        title('Kite forces in $\textbf{k}_K$')
    end
    xlabel('$t$ (s)')
    ylabel('Forces (N)')
    legend({'$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$'},'Location','east')
    ylim([-300 150])

%     sgtitle('Kite forces overview')
%% Traslation

    figure;
    subplot(2,1,1)
    plot(t, x(:,1), t, x(:,2))
%     title('Position of the sled in $$S_E$$')
    xlabel('$t$ (s)')
    ylabel('Sled position (m)')
    legend('$x_S$','$y_S$')
    ylim([-7 -3])
    
    subplot(2,1,2)
    plot(t, x(:,7), t, x(:,8), t, x(:,9))
%     title('Position of the kite in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Kite position (m)')
    legend({'$x_K$','$y_K$','$z_K$'},'location','southeast')
    ylim([-51 2])
    
    [F_S2, ~] = separarElementos(F_S(:,2));
    [~, F_S3] = separarElementos(F_S(:,3));
    [t2, ~] = separarElementos(t);
    [~, t3] = separarElementos(t);

    figure;
    subplot(2,1,1);
    plot(t, F_S(:,1), '-','Color',color.blue, 'LineWidth', 1); 
    hold on;  
    plot(t2, F_S2, 'o','Color',color.red, 'LineWidth', 2.5);  
    plot(t3, F_S3, 'square','Color',color.green, 'LineWidth', 1.5, 'MarkerSize', 3.5);
    hold off;
%     title('Sled force components in $S_E$')
    xlabel('$t$ (s)')
    ylabel('Sled force (N)')
    legend({'$\textbf{F}_S \cdot\textbf{i}_E$','$\textbf{F}_S\cdot\textbf{j}_E$','$\textbf{F}_S\cdot\textbf{k}_E$'},'location','northeast')
    adjust_ylim(gca, F_S);
    
    [F_K2, ~] = separarElementos(F_K(:,2));
    [~, F_K3] = separarElementos(F_K(:,3));
    [t2, ~] = separarElementos(t);
    [~, t3] = separarElementos(t);

    subplot(2,1,2);
    plot(t, F_K(:,1), '-','Color',color.blue, 'LineWidth', 1); 
    hold on;  
    plot(t2, F_K2, 'o','Color',color.red, 'LineWidth', 2.5);  
    plot(t3, F_K3, 'square','Color',color.green, 'LineWidth', 1.5, 'MarkerSize', 3.5);
    hold off;
    xlabel('$t$ (s)')
    ylabel('Kite force (N)')
    if frame == 0
%         title('Kite force components in $S_E$')
        legend({'$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$'},'location','northeast')
    else
%         title('Kite force components in $S_K$')
        legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
    end
    adjust_ylim(gca, F_K);
    
    figures = findall(0, 'Type', 'figure');

    %%
    function [pares, impares] = separarElementos(vector)
    % Separar elementos pares e impares del vector columna
    % vector: vector columna de entrada
    % pares: vector con los elementos en las posiciones pares
    % impares: vector con los elementos en las posiciones impares

    % Verificar que el vector sea columna
    if size(vector, 2) ~= 1
        error('El vector de entrada debe ser un vector columna.');
    end

    % Obtener los índices de los elementos pares e impares
    indices_pares = 2:2:length(vector);
    indices_impares = 1:2:length(vector);

    % Extraer los elementos correspondientes
    pares = vector(indices_pares);
    impares = vector(indices_impares);
end

end

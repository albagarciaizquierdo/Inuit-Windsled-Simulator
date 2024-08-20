function fun_draw_system(p,x,title)
    % Function that draws the current system given a state vector x
    %% Color settings
    color_green = [0 0.6 0];
    color_blue = [0 0.6 0.8];
    color_red = [1 0 0];

    %% State vector
    x_S = x(1);
    y_S = x(2);
    psi_S = x(3);
    x_K = x(7);
    y_K = x(8);
    z_K = x(9);
    phi = x(10);
    theta = x(11);
    psi_K = x(12);

    %% Earth reference frame S_E
    % Origin O_E
    OE = [0, 0, 0];

    figure;
    hold on;
    plot3(OE(1), OE(2), OE(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    text(OE(1), OE(2), OE(3), '$O_E$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Axes S_E
    i_E = [1, 0, 0];
    j_E = [0, 1, 0];
    k_E = [0, 0, 1];

    %% Sled reference frame
    % Origin O_S
    OE_OS = x_S * i_E + y_S * j_E;
    OS = OE_OS - OE;
    plot3(OS(1), OS(2), OS(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    if isequal(OS,OE) % Case of coincident OS and OE
        text(OS(1), OS(2)-1, OS(3), '=$O_S$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    else % General case
        text(OS(1), OS(2), OS(3), '$O_S$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end 
    
    % Axes in S_E
    i_S = [1, 0, 0];
    j_S = [0, 1, 0];
    k_S = [0, 0, 1];
    R_ES = fun_R_ES(psi_S);
    i_S = R_ES * i_S';
    j_S = R_ES * j_S';
    k_S = R_ES * k_S';
    
    if isequal(OS,OE) && isequal(i_S, i_E') && isequal(j_S, j_E') && isequal(k_S, k_E')
    quiver3(OS(1), OS(2), OS(3), i_S(1), i_S(2), i_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OS(1), OS(2), OS(3), j_S(1), j_S(2), j_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OS(1), OS(2), OS(3), k_S(1), k_S(2), k_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    text(OS(1) + i_S(1), OS(2) + i_S(2), OS(3) + i_S(3), '$\textbf{i}_S$', 'FontSize', 12, 'Color', color_blue);
    text(OS(1) + j_S(1), OS(2) + j_S(2), OS(3) + j_S(3), '$\textbf{j}_S$', 'FontSize', 12, 'Color', color_blue);
    text(OS(1) + k_S(1), OS(2) + k_S(2), OS(3) + k_S(3), '$\textbf{k}_S$', 'FontSize', 12, 'Color', color_blue);
    else
    % S_E
    quiver3(OE(1), OE(2), OE(3), i_E(1), i_E(2), i_E(3), 'Color', color_green, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OE(1), OE(2), OE(3), j_E(1), j_E(2), j_E(3), 'Color', color_green, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OE(1), OE(2), OE(3), k_E(1), k_E(2), k_E(3), 'Color', color_green, 'LineWidth', 1, 'MaxHeadSize', 1);
    text(i_E(1), i_E(2), i_E(3), '$\textbf{i}_E$', 'FontSize', 12, 'Color', color_green);
    text(j_E(1), j_E(2), j_E(3), '$\textbf{j}_E$', 'FontSize', 12, 'Color', color_green);
    text(k_E(1), k_E(2), k_E(3), '$\textbf{k}_E$', 'FontSize', 12, 'Color', color_green);
    
    % Labels S_S
    quiver3(OS(1), OS(2), OS(3), i_S(1), i_S(2), i_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OS(1), OS(2), OS(3), j_S(1), j_S(2), j_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OS(1), OS(2), OS(3), k_S(1), k_S(2), k_S(3), 'Color', color_blue, 'LineWidth', 1, 'MaxHeadSize', 1);
    text(OS(1) + i_S(1), OS(2) + i_S(2), OS(3) + i_S(3), '$\textbf{i}_S$', 'FontSize', 12, 'Color', color_blue);
    text(OS(1) + j_S(1), OS(2) + j_S(2), OS(3) + j_S(3), '$\textbf{j}_S$', 'FontSize', 12, 'Color', color_blue);
    text(OS(1) + k_S(1), OS(2) + k_S(2), OS(3) + k_S(3), '$\textbf{k}_S$', 'FontSize', 12, 'Color', color_blue);
    end
    
    %% SLED
    % Sled vertices
    A = p.l / 2 * i_S - p.w_S / 2 * j_S + OS';
    B = p.l / 2 * i_S + p.w_S / 2 * j_S + OS';
    C = -p.l / 2 * i_S + p.w_S / 2 * j_S + OS';
    D = -p.l / 2 * i_S - p.w_S / 2 * j_S + OS';

    % Sleed coordinates
    sled_x = [A(1), B(1), C(1), D(1), A(1)];
    sled_y = [A(2), B(2), C(2), D(2), A(2)];
    sled_z = [A(3), B(3), C(3), D(3), A(3)];

    fill3(sled_x, sled_y, sled_z, 'b', 'FaceAlpha', 0.1); % Shade the sled area
    plot3(sled_x, sled_y, sled_z, 'k', 'LineWidth', 1);

    %% Kite reference frame S_K
    % Origin O_K
    OE_OK = x_K * i_E + y_K * j_E + z_K * k_E;
    OK = OE_OK + OE;
    plot3(OK(1), OK(2), OK(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    text(OK(1), OK(2), OK(3), '$O_K$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Axes in S_K
    i_K = [1, 0, 0];
    j_K = [0, 1, 0];
    k_K = [0, 0, 1];
    R_EK = fun_R_EK(phi, theta, psi_K);
    i_K = R_EK * i_K';
    j_K = R_EK * j_K';
    k_K = R_EK * k_K';

    quiver3(OK(1), OK(2), OK(3), i_K(1), i_K(2), i_K(3), 'Color', color_red, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OK(1), OK(2), OK(3), j_K(1), j_K(2), j_K(3), 'Color', color_red, 'LineWidth', 1, 'MaxHeadSize', 1);
    quiver3(OK(1), OK(2), OK(3), k_K(1), k_K(2), k_K(3), 'Color', color_red, 'LineWidth', 1, 'MaxHeadSize', 1);

    % Labels
    text(OK(1) + i_K(1), OK(2) + i_K(2), OK(3) + i_K(3), '$\textbf{i}_K$', 'FontSize', 12, 'Color', color_red);
    text(OK(1) + j_K(1), OK(2) + j_K(2), OK(3) + j_K(3), '$\textbf{j}_K$', 'FontSize', 12, 'Color', color_red);
    text(OK(1) + k_K(1), OK(2) + k_K(2), OK(3) + k_K(3), '$\textbf{k}_K$', 'FontSize', 12, 'Color', color_red);

    %% Attachement points
    [ASp_AKp, ASp_AKm, up, um, ASp_OS, ASm_OS, OK_AKp, OK_AKm, T_ASp, T_ASm] = fun_tethers(p, x);
    ASp = OS - ASp_OS';
    ASm = OS - ASm_OS';
    AKp = ASp_AKp' + ASp;
    AKm = ASp_AKm' + ASm;

    plot3(ASp(1), ASp(2), ASp(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    plot3(ASm(1), ASm(2), ASm(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    plot3(AKp(1), AKp(2), AKp(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    plot3(AKm(1), AKm(2), AKm(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');

    text(ASp(1), ASp(2), ASp(3), '$A_S^+$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(ASm(1), ASm(2), ASm(3), '$A_S^-$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(AKp(1), AKp(2), AKp(3), '$A_K^+$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(AKm(1), AKm(2), AKm(3), '$A_K^-$', 'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    %% Lines
    plot3([ASp(1), AKp(1)], [ASp(2), AKp(2)], [ASp(3), AKp(3)], 'k', 'LineWidth', 1);
    plot3([ASm(1), AKm(1)], [ASm(2), AKm(2)], [ASm(3), AKm(3)], 'k', 'LineWidth', 1);

    %% Kite (represented with a rectangle contained in plane Y_K-Z_K)
    % Kite dimensions
    % Option 1: Square kite
    w_K = sqrt(p.S);
    l_K = sqrt(p.S);
    % Option 2: According to attachement points
    w_K = norm(AKp - AKm);
    l_K = p.S / w_K;

    % Kite vertices
    E = -w_K / 2 * j_K + l_K / 2 * k_K + OK';
    F = w_K / 2 * j_K + l_K / 2 * k_K + OK';
    G = +w_K / 2 * j_K - l_K / 2 * k_K + OK';
    H = -w_K / 2 * j_K - l_K / 2 * k_K + OK';

    % Kite coordinates
    kite_x = [E(1), F(1), G(1), H(1), E(1)];
    kite_y = [E(2), F(2), G(2), H(2), E(2)];
    kite_z = [E(3), F(3), G(3), H(3), E(3)];

    fill3(kite_x, kite_y, kite_z, 'r', 'FaceAlpha', 0.1); % Shade the sled area
    plot3(kite_x, kite_y, kite_z, 'k', 'LineWidth', 1);

    %% Plot3 configuration
    set(groot,'defaultAxesFontSize', 11);
    set(groot,'defaultLineLineWidth', 1);
    set(groot,'defaultAxesColorOrder', [0 0 1;1 0 0; 0 0.6 0;1 0.7 0.1; 0.5 0.2 0.6]);
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
    
    xlabel('$X_{E}$');
    ylabel('$Y_{E}$');
    zlabel('$Z_{E}$');
    grid on;
    sgtitle(title);

    % Legend with state vector
    legend_str = sprintf('State vector:\n $x_S$ = %.2f m\n $y_S$ = %.2f m\n $\\psi_S$ = %.2f deg\n $x_K$ = %.2f m\n $y_K$ = %.2f m\n $z_K$ = %.2f m\n $\\phi_K$ = %.2f deg\n $\\theta_K$ = %.2f deg\n $\\psi_K$ = %.2f deg', x_S, y_S, rad2deg(psi_S), x_K, y_K, z_K, rad2deg(phi), rad2deg(theta), rad2deg(psi_K));
    dim = [0.15 0.75 0.3 0.15]; % [x y width height] of the textbox in normalized units
    annotation('textbox', dim, 'String', legend_str, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Interpreter', 'latex', 'FontSize', 12);

    % Axes configuration
    ax = gca; % Obtain current axes object
    ax.XDir = 'reverse'; % Reverse direction for x
    ax.YDir = 'normal'; % Normal y direction
    ax.ZDir = 'reverse'; % Reverse direction for z
%     axis equal;
    view(3);
    hold off;
end

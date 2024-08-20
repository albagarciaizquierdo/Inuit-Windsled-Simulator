function figures = fun_draw_results(x, t, p, frame)
%% Frame for the kite
% frame = 1 --> SK
% frame = 0 --> SE

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
        ylim(ax, [-1e-10, 1e-10]);
    end
end

%% SLED
figure;
subplot(2,2,1)
plot(t, x(:,1), t, x(:,2))
title('Position of the sled in $$S_E$$')
xlabel('$t$ (s)')
ylabel('Position (m)')
legend('$x_S$','$y_S$')
adjust_ylim(gca, x(:, 1:2));

subplot(2,2,2)
plot(t, rad2deg(x(:,3)))
title('$\psi_S$')
xlabel('$t$ (s)')
ylabel('$\psi_S$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 3)));

subplot(2,2,3)
plot(t, x(:,4), t, x(:,5))
title('Sled velocity')
xlabel('$t$ (s)')
ylabel('$v^S$ (m/s)')
legend('$v_x^S$','$v_y^S$')
adjust_ylim(gca, x(:, 4:5));

subplot(2,2,4)
plot(t, x(:,6))
title('Sled angular velocity $r_S$')
xlabel('$t$ (s)')
ylabel('$r_S$ (rad/s)')
adjust_ylim(gca, x(:, 6));

sgtitle('Sled Overview')

figure()
subplot(2,2,1)
plot(t, x(:,4))
%     title('Velocity in i_E$')
xlabel('$t$ (s)')
ylabel('$v_x^S$ (m/s)')
ylim([-2*10^(-3),0])
%     adjust_ylim(gca, x(:, 1:2));

subplot(2,2,2)
plot(t, rad2deg(x(:,3)))
%     title('$\psi_S$')
xlabel('$t$ (s)')
ylabel('$\psi_S$ (deg)')
%     adjust_ylim(gca, rad2deg(x(:, 3)));

subplot(2,2,3)
plot(t, x(:,5))
%     title('Sled velocity')
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
%     adjust_ylim(gca, x(:, 4:5));

subplot(2,2,4)
plot(t, x(:,6))
%     title('Sled angular velocity $r_S$')
xlabel('$t$ (s)')
ylabel('$r_S$ (rad/s)')
%     adjust_ylim(gca, x(:, 6));


%% KITE
figure;
subplot(2,2,1)
plot(t, x(:,7), t, x(:,8), t, x(:,9))
title('Position $x_K,y_K,z_K$ of the kite in $S_E$')
xlabel('$t$ (s)')
ylabel('Position (m)')
legend('$x_K$','$y_K$','$z_K$')
adjust_ylim(gca, x(:, 7:9));

subplot(2,2,2)
plot(t, rad2deg(x(:,10)), t, rad2deg(x(:,11)), t, rad2deg(x(:,12)),'--')
title('Angles of the kite')
xlabel('$t$ (s)')
ylabel('Angles (deg)')
legend('$\phi_K$','$\theta_K$','$\psi_K$')
adjust_ylim(gca, rad2deg(x(:, 10:12)));

subplot(2,2,3)
plot(t, x(:,13), t, x(:,14), t, x(:,15),'--')
title('Kite velocity in $S_K$')
xlabel('$t$ (s)')
ylabel('Velocity (m/s)')
legend('$u$','$v$','$w$')
adjust_ylim(gca, rad2deg(x(:, 13:15)));

subplot(2,2,4)
plot(t, rad2deg(x(:,16)), t, rad2deg(x(:,17)), t, rad2deg(x(:,18)),'--')
title('Angular velocity of the kite in $S_K$')
xlabel('$t$ (s)')
ylabel('Angular velocity (rad/s)')
legend('$p_K$','$q_K$','$r_K$')
adjust_ylim(gca, rad2deg(x(:, 16:18)));

sgtitle('Kite Overview')
%% Positions x y
figure;
subplot(2,2,1)
plot(t, x(:,1))
xlabel('$t$ (s)')
ylabel('$x_S$ (m)')

subplot(2,2,2)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')

subplot(2,2,3)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)') 

subplot(2,2,4)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')

%% Position kite

figure;
subplot(3,1,1)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')

subplot(3,1,2)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')

subplot(3,1,3)
plot(t, x(:,9))
xlabel('$t$ (s)')
ylabel('$z_K$ (m)')

%% Kite's height

figure;
plot(t, x(:,9))
xlabel('$t$ (s)')
ylabel('$z_K$ (m)')

%% Euler angles, pqr, uvw
figure;
subplot(3,1,1)
plot(t, rad2deg(x(:,10)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$\phi$ (deg)', 'Interpreter', 'latex')

subplot(3,1,2)
plot(t, rad2deg(x(:,11)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$\theta$ (deg)', 'Interpreter', 'latex')

subplot(3,1,3)
plot(t, rad2deg(x(:,12)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$\psi_K$ (deg)', 'Interpreter', 'latex')

figure;
subplot(3,1,1)
plot(t, rad2deg(x(:,16)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$p_K$ (rad/s)', 'Interpreter', 'latex')

subplot(3,1,2)
plot(t, rad2deg(x(:,17)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$q_K$ (rad/s)', 'Interpreter', 'latex')

subplot(3,1,3)
plot(t, rad2deg(x(:,18)))
xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$r_K$ (rad/s)', 'Interpreter', 'latex')

figure;
subplot(3,1,1)
plot(t, x(:,13))
xlabel('$t$ (s)')
ylabel('$u$ (m/s)')

subplot(3,1,2)
plot(t, x(:,14))
xlabel('$t$ (s)')
ylabel('$v$ (m/s)')

subplot(3,1,3)
plot(t, x(:,15))
xlabel('$t$ (s)')
ylabel('$w$ (m/s)')

%%
%     %% 3D SLED
%     figure()
%     plot3(x(:,1), x(:,2), t, 'LineWidth', 2);
%     grid on;
%     xlabel('$x_S$ (m)');
%     ylabel('$y_S$ (m)');
%     zlabel('$t$ (s)');
%     title('3D position of the sled in $S_E$');

%     %% 3D KITE
%     figure()
%     plot3(x(:,7), x(:,8), x(:,9), 'LineWidth', 2);
%     grid on;
%     xlabel('$x_K$ (m)');
%     ylabel('$y_K$ (m)');
%     zlabel('$z_K$ (m)');
%     title('3D position of the kite in $S_K$');

%% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m] = fun_get_results(x,p,frame);

% Aerodynamic
figure;
subplot(1,2,1)
plot(t, rad2deg(alpha), t, rad2deg(beta))
title('AoA and sideslip angles')
xlabel('$t$ (s)')
ylabel('Angles (deg)')
legend('$\alpha$','$\beta$')
adjust_ylim(gca, rad2deg([alpha, beta]));

subplot(1,2,2)
plot(t, v_a(:,1), t, v_a(:,2), t, v_a(:,3))
if frame == 0
    title('Aerodynamic velocity with time in $S_E$')
    legend('$\textbf{v}_a \cdot\textbf{i}_E$','$\textbf{v}_a\cdot\textbf{j}_E$','$\textbf{v}_a\cdot\textbf{k}_E$')
else
    title('Aerodynamic velocity with time in $S_K$')
    legend('$\textbf{v}_a \cdot\textbf{i}_K$','$\textbf{v}_a\cdot\textbf{j}_K$','$\textbf{v}_a\cdot\textbf{k}_K$')
end
xlabel('$t$ (s)')
ylabel('$\textbf{v}_a$ (m/s)')
adjust_ylim(gca, v_a);

% Tethers length
figure()
plot(t, l_p)
title('Tethers length with time')
xlabel('$t$ (s)')
ylabel('Tethers length (m)')
adjust_ylim(gca, l_p);

% Sled forces
figure;
subplot(2,2,1);
plot(t, F_S(:,1), t, F_S(:,2), t, F_S(:,3),'--')
title('Sled force with time')
xlabel('$t$ (s)')
ylabel('Sled force (N)')
legend('$\textbf{F}_S \cdot\textbf{i}_E$','$\textbf{F}_S\cdot\textbf{j}_E$','$\textbf{F}_S\cdot\textbf{k}_E$')
adjust_ylim(gca, F_S);

subplot(2,2,2);
plot(t, N(:,1), t, W_S(:,1),'--', t, F_r(:,1), t, T_ASp(:,1), t, T_ASm(:,1),'--')
title('Sled forces in $\textbf{i}_E$')
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [N(:,1), W_S(:,1), F_r(:,1), T_ASp(:,1), T_ASm(:,1)]);

subplot(2,2,3);
plot(t, N(:,2), t, W_S(:,2),'-.', t, F_r(:,2),':', t, T_ASp(:,2), t, T_ASm(:,2),'--')
title('Sled forces in $\textbf{j}_E$')
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [N(:,2), W_S(:,2), F_r(:,2), T_ASp(:,2), T_ASm(:,2)]);

subplot(2,2,4);
plot(t, N(:,3), t, W_S(:,3), t, F_r(:,3), t, T_ASp(:,3), t, T_ASm(:,3),'--')
title('Sled forces in $\textbf{k}_E$')
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$N$','$W_S$','$F_r$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [N(:,3), W_S(:,3), F_r(:,3), T_ASp(:,3), T_ASm(:,3)]);

sgtitle('Sled Forces Overview')

% Kite forces
figure;
subplot(2,2,1);
plot(t, F_K(:,1), t, F_K(:,2), t, F_K(:,3))
title('Kite force with time')
xlabel('$t$ (s)')
ylabel('Sled force (N)')
if frame == 0
    legend('$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$')
else
    legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
end
adjust_ylim(gca, F_K);

subplot(2,2,2);
plot(t, F_a(:,1), t, W_K(:,1), t, T_AKp(:,1), t, T_AKm(:,1),'--')
if frame == 0
    title('Kite forces in $\textbf{i}_E$')
else
    title('Kite forces in $\textbf{i}_K$')
end
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [F_a(:,1), W_K(:,1), T_AKp(:,1), T_AKm(:,1)]);

subplot(2,2,3);
plot(t, F_a(:,2), t, W_K(:,2), t, T_AKp(:,2), t, T_AKm(:,2),'--')
if frame == 0
    title('Kite forces in $\textbf{j}_E$')
else
    title('Kite forces in $\textbf{j}_K$')
end
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [F_a(:,2), W_K(:,2), T_AKp(:,2), T_AKm(:,2)]);

subplot(2,2,4);
plot(t, F_a(:,3), t, W_K(:,3), t, T_AKp(:,3), t, T_AKm(:,3),'--')
if frame == 0
    title('Kite forces in $\textbf{j}_E$')
else
    title('Kite forces in $\textbf{j}_K$')
end
xlabel('$t$ (s)')
ylabel('Forces (N)')
legend('$F_a$','$W_K$','$T_{A_S}^+$','$T_{A_S}^-$')
adjust_ylim(gca, [F_a(:,3), W_K(:,3), T_AKp(:,3), T_AKm(:,3)]);

sgtitle('Kite Forces Overview')

%% VARYING XS
figure;
subplot(2,1,1)
plot(t, x(:,1))
xlabel('$t$ (s)')
ylabel('$x_S$ (m)')
%     ylim([-2*10^(-3) 0])

subplot(2,1,2)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')
adjust_ylim(gca, x(:, 1:2));

figure;
subplot(2,1,1)
plot(t, x(:,13), t, x(:,14), t, x(:,15))
xlabel('$t$ (s)')
ylabel('Kite velocity (m/s)')
legend('$u$','$v$','$w$')
adjust_ylim(gca, rad2deg(x(:, 13:15)));

subplot(2,1,2);
plot(t, F_K(:,1),'Color',color.blue, 'LineWidth', 1.5)
hold on
plot(t(1:10:end), F_K((1:10:end),2),'o','Color',color.red, 'LineWidth', 2.5)
plot(t, F_K(:,3))
xlabel('$t$ (s)')
ylabel('Kite force (N)')
if frame == 0
    legend('$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$')
else
    legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
end
adjust_ylim(gca, F_K);

%% VARYING YS
figure;
subplot(2,1,1)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)')

subplot(2,1,2)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')

figure;
subplot(2,1,1)
plot(t, x(:,4))
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
ylim([-4*10^(-4) 4*10^(-4)])

subplot(2,1,2)
plot(t, x(:,14))
xlabel('$t$ (s)')
ylabel('$v$ (m/s)')
adjust_ylim(gca, rad2deg(x(:, 14)));

figure;
subplot(2,1,1);
plot(t(1:30:end), F_S((1:30:end),1),'o','Color',color.blue, 'LineWidth', 2.5)
hold on
plot(t, F_S(:,2),'Color',color.red)
plot(t, F_S(:,3),'--')    
xlabel('$t$ (s)')
ylabel('Sled force (N)')
legend('$\textbf{F}_S \cdot\textbf{i}_E$','$\textbf{F}_S\cdot\textbf{j}_E$','$\textbf{F}_S\cdot\textbf{k}_E$')

subplot(2,1,2);
plot(t(1:30:end), F_K((1:30:end),1),'o','Color',color.blue, 'LineWidth', 2.5)
hold on
plot(t, F_K(:,2),'Color',color.red)
plot(t, F_K(:,3),'--')
xlabel('$t$ (s)')
ylabel('Kite force (N)')
if frame == 0
    legend('$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$')
else
    legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
end
adjust_ylim(gca, F_K);

%%
figures = findall(0, 'Type', 'figure');
end

function figures = fun_draw_results_stab(x, t, p, frame)
%% Description: 
% This function draws the results of state vector x vs time specifically
% for the stability analysis with numerical integration
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
    if max(abs(data(:))) < 10^-5
        ylim(ax, [-1e-4, 1e-4]);
    end
end

%% RESULTS
[F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,...
    OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m] = fun_get_results(x,p,frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER
%% 28
figure;
subplot(2,1,1)
plot(t, rad2deg(alpha))
xlabel('$t$ (s)')
ylabel('$\alpha$ (deg)')
adjust_ylim(gca, rad2deg([alpha]));

subplot(2,1,2)
plot(t, rad2deg(beta))
xlabel('$t$ (s)')
ylabel('$\beta$ (deg)')
adjust_ylim(gca, rad2deg([beta]));

%% 27
figure;
subplot(2,1,1)
plot(t, l_p)
xlabel('$t$ (s)')
ylabel('$l^+$ (m)')
adjust_ylim(gca, l_p);

subplot(2,1,2)
plot(t, l_p)
xlabel('$t$ (s)')
ylabel('$l^-$ (m)')
adjust_ylim(gca, l_m);

%% COMBINATION
%% 26
figure;
subplot(2,1,1);
% plot(t(1:30:end), F_S((1:30:end),1),'o','Color',color.blue, 'LineWidth', 2.5)
plot(t, F_S(:,1))
hold on
plot(t, F_S(:,2))
% plot(t, F_S(:,3),'--')    
plot(t, F_S(:,3))
xlabel('$t$ (s)')
ylabel('Sled force (N)')
legend('$\textbf{F}_S \cdot\textbf{i}_E$','$\textbf{F}_S\cdot\textbf{j}_E$','$\textbf{F}_S\cdot\textbf{k}_E$')
adjust_ylim(gca, F_S);


subplot(2,1,2);
plot(t, F_K(:,1),'Color',color.blue, 'LineWidth', 1.5)
hold on
% plot(t(1:10:end), F_K((1:10:end),2),'o','Color',color.red, 'LineWidth', 2.5)
plot(t, F_K(:,2))
plot(t, F_K(:,3))
xlabel('$t$ (s)')
ylabel('Kite force (N)')
if frame == 0
    legend('$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$')
else
    legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
end
adjust_ylim(gca, F_K);

%% 25
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,6)))
xlabel('$t$ (s)')
ylabel('$r_S$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:,6)))

subplot(2,1,2)
plot(t, rad2deg(x(:,18)))
xlabel('$t$ (s)')
ylabel('$r_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 18)))

%% 24
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,3)))
xlabel('$t$ (s)')
ylabel('$\psi_S$ (deg)')
adjust_ylim(gca, rad2deg(x(:,3)));
 
subplot(2,1,2)
plot(t, rad2deg(x(:,12)))
xlabel('$t$ (s)')
ylabel('$\psi_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:,12)));

%% 23
figure;
subplot(2,1,1)
plot(t, x(:,4), t, x(:,5))
xlabel('$t$ (s)')
ylabel('Sled velocity (m/s)')
legend('$v_x^S$','$v_y^S$')
adjust_ylim(gca, (x(:, 4:5)));
 
subplot(2,1,2)
plot(t, x(:,13), t, x(:,14), t, x(:,15))
xlabel('$t$ (s)')
ylabel('Kite velocity (m/s)')
legend('$u$','$v$','$w$')
adjust_ylim(gca, (x(:, 13:15)));

%% 22
figure;
subplot(2,1,1)
plot(t, x(:,5))
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
adjust_ylim(gca,(x(:,5)));

subplot(2,1,2)
plot(t, x(:,14))
xlabel('$t$ (s)')
ylabel('$v$ (m/s)')
adjust_ylim(gca,(x(:,14)));

%% 21
figure;
subplot(2,1,1)
plot(t, x(:,4))
xlabel('$t$ (s)')
ylabel('$v_x^S$ (m/s)')
adjust_ylim(gca,(x(:,4)));

subplot(2,1,2)
plot(t, x(:,13))
xlabel('$t$ (s)')
ylabel('$u$ (m/s)')
adjust_ylim(gca,(x(:,13)));

%% 20
figure;
subplot(2,1,1)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)')
adjust_ylim(gca,(x(:,2)));
% ylim([-2.6*10^(-4)  -2.5*10^(-4)])

subplot(2,1,2)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')
adjust_ylim(gca,(x(:,8)));
% ylim([-2.6*10^(-4)  -2.5*10^(-4)])

%% 19
figure;
subplot(2,1,1)
plot(t, x(:,1))
xlabel('$t$ (s)')
ylabel('$x_S$ (m)')
adjust_ylim(gca, x(:,1));

subplot(2,1,2)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')
adjust_ylim(gca, x(:, 7));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KITE
%% 18
figure;
subplot(2,1,1)
plot(t, x(:,13), t, x(:,14), t, x(:,15))
xlabel('$t$ (s)')
ylabel('Kite velocity (m/s)')
legend('$u$','$v$','$w$')
adjust_ylim(gca, (x(:, 13:15)));

subplot(2,1,2);
plot(t, F_K(:,1))
hold on
% plot(t(1:10:end), F_K((1:10:end),2),'o','Color',color.red, 'LineWidth', 2.5)
plot(t, F_K(:,2))
plot(t, F_K(:,3))
xlabel('$t$ (s)')
ylabel('Kite force (N)')
if frame == 0
    legend('$\textbf{F}_K \cdot\textbf{i}_E$','$\textbf{F}_K\cdot\textbf{j}_E$','$\textbf{F}_K\cdot\textbf{k}_E$')
else
    legend('$\textbf{F}_K \cdot\textbf{i}_K$','$\textbf{F}_K\cdot\textbf{j}_K$','$\textbf{F}_K\cdot\textbf{k}_K$')
end
adjust_ylim(gca, F_K);

%% 17
figure;
subplot(3,1,1)
plot(t, F_K(:,1))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_K\cdot\textbf{i}_E$ (N)')
adjust_ylim(gca, F_K(:,1));

subplot(3,1,2)
plot(t, F_K(:,2))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_K\cdot\textbf{j}_E$ (N)')
adjust_ylim(gca, F_K(:,2));

subplot(3,1,3)
plot(t, F_K(:,3))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_K\cdot\textbf{k}_E$ (N)')
adjust_ylim(gca, F_K(:,3));

%% 16
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,12)))
xlabel('$t$ (s)')
ylabel('$\psi_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 12)));

subplot(2,1,2)
plot(t, rad2deg(x(:,18)))
xlabel('$t$ (s)')
ylabel('$r_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 18)));
%% 15
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,11)))
xlabel('$t$ (s)')
ylabel('$\theta_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 11)));

subplot(2,1,2)
plot(t, rad2deg(x(:,17)))
xlabel('$t$ (s)')
ylabel('$q_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 17)));

%% 14
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,10)))
xlabel('$t$ (s)')
ylabel('$\phi_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 10)));

subplot(2,1,2)
plot(t, rad2deg(x(:,16)))
xlabel('$t$ (s)')
ylabel('$p_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 16)));

%% 13
figure;
subplot(2,1,1)
plot(t, x(:,9))
xlabel('$t$ (s)')
ylabel('$z_K$ (m)')
adjust_ylim(gca, x(:, 9));

subplot(2,1,2)
plot(t, x(:,15))
xlabel('$t$ (s)')
ylabel('$w$ (m/s)')
adjust_ylim(gca, x(:, 15));
%% 12
figure;
subplot(2,1,1)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')
adjust_ylim(gca, x(:, 8));

subplot(2,1,2)
plot(t, x(:,14))
xlabel('$t$ (s)')
ylabel('$v$ (m/s)')
adjust_ylim(gca, x(:, 14));
%% 11
figure;
subplot(2,1,1)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')
adjust_ylim(gca, x(:, 7));

subplot(2,1,2)
plot(t, x(:,13))
xlabel('$t$ (s)')
ylabel('$u$ (m/s)')
adjust_ylim(gca, x(:, 13));
%% 10
subplot(3,1,1)
plot(t, rad2deg(x(:,16)))
xlabel('$t$ (s)')
ylabel('$p_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 16)));

subplot(3,1,2)
plot(t, rad2deg(x(:,17)))
xlabel('$t$ (s)')
ylabel('$q_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 17)));

subplot(3,1,3)
plot(t, rad2deg(x(:,18)))
xlabel('$t$ (s)')
ylabel('$r_K$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 18)));
%% 9.2
figure;
plot(t, rad2deg(x(:,11)))
xlabel('$t$ (s)')
ylabel('$\theta_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 11)));
%% 9
figure;
subplot(3,1,1)
plot(t, rad2deg(x(:,10)))
xlabel('$t$ (s)')
ylabel('$\phi_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 10)));

subplot(3,1,2)
plot(t, rad2deg(x(:,11)))
xlabel('$t$ (s)')
ylabel('$\theta_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 11)));

subplot(3,1,3)
plot(t, rad2deg(x(:,12)))
xlabel('$t$ (s)')
ylabel('$\psi_K$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 12)));
%% 8.2
figure;
subplot(2,1,1)
plot(t, x(:,13))
xlabel('$t$ (s)')
ylabel('$u$ (m/s)')
adjust_ylim(gca, x(:, 13));

subplot(2,1,2)
plot(t, x(:,15))
xlabel('$t$ (s)')
ylabel('$w$ (m/s)')
adjust_ylim(gca, x(:, 15));
%% 8
figure;
subplot(3,1,1)
plot(t, x(:,13))
xlabel('$t$ (s)')
ylabel('$u$ (m/s)')
adjust_ylim(gca, x(:, 13));

subplot(3,1,2)
plot(t, x(:,14))
xlabel('$t$ (s)')
ylabel('$v$ (m/s)')
adjust_ylim(gca, x(:, 14));

subplot(3,1,3)
plot(t, x(:,15))
xlabel('$t$ (s)')
ylabel('$w$ (m/s)')
adjust_ylim(gca, x(:, 15));
%% 7.2
figure;
subplot(2,1,1)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')
adjust_ylim(gca, x(:,7));

subplot(2,1,2)
plot(t, x(:,9))
xlabel('$t$ (s)')
ylabel('$z_K$ (m)')
adjust_ylim(gca, x(:,9));
%% 7
figure;
subplot(3,1,1)
plot(t, x(:,7))
xlabel('$t$ (s)')
ylabel('$x_K$ (m)')
adjust_ylim(gca, x(:,7));

subplot(3,1,2)
plot(t, x(:,8))
xlabel('$t$ (s)')
ylabel('$y_K$ (m)')
adjust_ylim(gca, x(:,8));

subplot(3,1,3)
plot(t, x(:,9))
xlabel('$t$ (s)')
ylabel('$z_K$ (m)')
adjust_ylim(gca, x(:,9));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLED
%% 6 
figure;
subplot(3,1,1)
plot(t, F_S(:,1))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_S\cdot\textbf{i}_E$ (N)')
adjust_ylim(gca, F_S(:,1));

subplot(3,1,2)
plot(t, F_S(:,2))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_S\cdot\textbf{j}_E$ (N)')
adjust_ylim(gca, F_S(:,2));

subplot(3,1,3)
plot(t, F_S(:,3))
xlabel('$t$ (s)')
ylabel('$\textbf{F}_S\cdot\textbf{k}_E$ (N)')
adjust_ylim(gca, F_S(:,3));
%% 5.5
figure;
subplot(3,1,1)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)')
adjust_ylim(gca, x(:, 2));

subplot(3,1,2)
plot(t, x(:,5))
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
adjust_ylim(gca, x(:, 5));

subplot(3,1,3)
plot(t, rad2deg(x(:,3)))
xlabel('$t$ (s)')
ylabel('$\psi_S$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 3)));
%% 5
figure;
subplot(2,1,1)
plot(t, rad2deg(x(:,3)))
xlabel('$t$ (s)')
ylabel('$\psi_S$ (deg)')
adjust_ylim(gca, rad2deg(x(:, 3)));

subplot(2,1,2)
plot(t, rad2deg(x(:,6)))
xlabel('$t$ (s)')
ylabel('$r_S$ (deg/s)')
adjust_ylim(gca, rad2deg(x(:, 6)));
%% 4
figure;
subplot(2,1,1)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)')
adjust_ylim(gca, x(:, 2));

subplot(2,1,2)
plot(t, x(:,5))
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
adjust_ylim(gca, x(:, 5));
%% 3
figure;
subplot(2,1,1)
plot(t, x(:,1))
xlabel('$t$ (s)')
ylabel('$x_S$ (m)')
adjust_ylim(gca, x(:, 1));

subplot(2,1,2)
plot(t, x(:,4))
xlabel('$t$ (s)')
ylabel('$v_x^S$ (m/s)')
adjust_ylim(gca, x(:, 4));

%% 2
figure;
subplot(2,1,1)
plot(t, x(:,4))
xlabel('$t$ (s)')
ylabel('$v_x^S$ (m/s)')
adjust_ylim(gca, x(:, 4));

subplot(2,1,2)
plot(t, x(:,5))
xlabel('$t$ (s)')
ylabel('$v_y^S$ (m/s)')
adjust_ylim(gca, x(:, 5));
%% 1
figure;
subplot(2,1,1)
plot(t, x(:,1))
xlabel('$t$ (s)')
ylabel('$x_S$ (m)')
adjust_ylim(gca, x(:, 1));

subplot(2,1,2)
plot(t, x(:,2))
xlabel('$t$ (s)')
ylabel('$y_S$ (m)')
adjust_ylim(gca, x(:, 2));

adjust_ylim(gca, x(:, 1:2));
%%
figures = findall(0, 'Type', 'figure');
end

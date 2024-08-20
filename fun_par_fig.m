function figures = fun_par_fig(flag,p,points,vw_max,frame,x_red0_design,color)
%% Function description:
% Generate several figures for the parametric study based on the slection
% of the "flag" variable
%% Inputs:
% flag: Numeric value between 1 to 13 selecting the variable to be varied
% in the parametric study as follows:
    % Flag = 1 --> Kite position vs wind velocity with vw_min 
    % Flag = 2 --> Length and angles for filtered velocity vs vw
    % Flag = 3 --> All vs vw for different l0
    % Flag = 4 --> z_K vs vw for different m_S
    % Flag = 5 --> All vs vw for different S
    % Flag = 6 --> All vs vw for different yA_K
    % Flag = 7 --> All vs vw for different zA_K
    % Flag = 8 --> All vs vw for different xA_K
    % Flag = 9 --> All vs vw for different w_S
    % Flag = 10 --> All vs vw for different l_S
    % Flag = 11 --> z_K vs vw for different xA_S
    % Flag = 12 --> z_K vs vw for different yA_S
    % Flag = 13 --> z_K vs vw for different mu_s
% p --> parameters struct
% points --> number of points to be calculated
% vw_max --> maximum speed of wind desired to simulate
% frame --> selection of kite forces and moments reference frame
    % frame = 0 --> SE
    % frame = 1 --> SK
% x_red0_design --> initial guess (reduces state vector) for vw=10m/s
% color --> struct with colors definitions for graphs
%% Output
% Figures --> variable containing all the figures

%%
global p
switch flag
    case 1
        [x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_vw(p,points,x_red0_design,frame,vw_max);

        figure()
        hold on
        plot(v_w, x_eq(7,:), 'DisplayName', '$x_K$')
        plot(v_w, x_eq(8,:), 'DisplayName', '$y_K$')
        plot(v_w, x_eq(9,:), 'DisplayName', '$z_K$')
        plot(v_w_min, 0, 'k.', 'MarkerSize', 10, 'DisplayName', '$z_K=0$')
        y_limits = ylim;
        fill_area = fill([min(v_w), v_w_min, v_w_min, min(v_w)], ...
             [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
             [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        set(get(get(fill_area, 'Annotation'), 'LegendInformation'), ...
            'IconDisplayStyle', 'on');
        set(fill_area, 'DisplayName', 'Non-physical area');
        xline(v_w_min, 'k--', 'LineWidth', 1, 'DisplayName', ['$v_w|_{min} = ', num2str(v_w_min, '%.2f'), '$ m/s']);
        xlim([0 14])
        hold off
%         title('Kite position vs wind velocity')
        xlabel('$v_w$ (m/s)', 'Interpreter', 'latex')
        ylabel('Kite position (m)', 'Interpreter', 'latex')
        leg = legend('Location', 'best', 'Interpreter', 'latex');
        leg.Color = 'white';
        leg.Box = 'on';
        leg.EdgeColor = 'black';
        ax = gca;
        ax.Box = 'on';
        grid on
        
    case 2
        [x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
        % Filter above min velocity
        [v_w_filt,gamma_filt,alpha_filt,l_filt,x_eq_filt] = fun_par_filter(v_w_min,v_w,gamma,alpha,l,x_eq);

        figure()
        plot(v_w_filt,l_filt)
%         title('Tethers length vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$l$ (m)')

        figure()
        plot(v_w_filt,gamma_filt)
%         title('Elevation angle vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')

        figure()
        plot(v_w_filt,rad2deg(alpha_filt))
%         title('Angle of attack vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        ylim([6.1 9.7])

        figure()
        plot(gamma_filt,rad2deg(alpha_filt))
%         title('Angle of attack vs elevation angle')
        xlabel('$\gamma$ (deg)')
        ylabel('$\alpha$ (deg)')
        ylim([6.25 9.7])
       
    case 3
        l0= zeros(1,4);
        for j = 1:4
            l0(j) = p.l0;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w_1, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w_2, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,0.6*points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w_3, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,0.45*points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w_4, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,0.35*points,x_red0_design,frame,vw_max);
            end
            p.l0 = p.l0 - 10;
        end
        % z_k and min vel
        figure()
        plot(v_w_1,x_eq_1(9,:),v_w_2,x_eq_2(9,:),v_w_3,x_eq_3(9,:),v_w_4,x_eq_4(9,:),v_w_1,zeros(1,length(v_w_1)),'k--')
%         title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'pentagram', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 12);
        plot(v_w_min_2, 0, 'pentagram', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 9);
        plot(v_w_min_3, 0, 'pentagram', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 6);
        plot(v_w_min_4, 0, 'pentagram', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        leg = legend(['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']);
        leg.Color = 'white';
        leg.Box = 'on';
        leg.EdgeColor = 'black';
        xlim([0,14])
        hold off
        
        [v_w_1,gamma_1,alpha_1,l_1,x_eq_1] = fun_par_filter(v_w_min_1,v_w_1,gamma_1,alpha_1,l_1,x_eq_1);
        [v_w_2,gamma_2,alpha_2,l_2,x_eq_2] = fun_par_filter(v_w_min_2,v_w_2,gamma_2,alpha_2,l_2,x_eq_2);
        [v_w_3,gamma_3,alpha_3,l_3,x_eq_3] = fun_par_filter(v_w_min_3,v_w_3,gamma_3,alpha_3,l_3,x_eq_3);
        [v_w_4,gamma_4,alpha_4,l_4,x_eq_4] = fun_par_filter(v_w_min_4,v_w_4,gamma_4,alpha_4,l_4,x_eq_4);
        
        figure()
        plot(v_w_1,x_eq_1(7,:),v_w_2,x_eq_2(7,:),v_w_3,x_eq_3(7,:),v_w_4,x_eq_4(7,:))
%         title('$x_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$x_K$ (m)')
        legend(['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m'])
        xlim([3.2,14])
        
        figure()
        hold on
        plot(v_w_1, gamma_1, '-', 'MarkerSize', 3);
        plot(v_w_2, gamma_2, '*', 'MarkerSize', 3);
        plot(v_w_3, gamma_3, 'o', 'MarkerSize', 5);
        plot(v_w_4, gamma_4, 'square', 'MarkerSize', 7);
        hold off
%         title('Elevation angle vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        xlim([3.2,14])
        legend(['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m'])
        ax = gca;
        ax.Box = 'on';
        
        figure()
        hold on
        plot(v_w_1, rad2deg(alpha_1), '-', 'MarkerSize', 3);
        plot(v_w_2, rad2deg(alpha_2), '*', 'MarkerSize', 3);
        plot(v_w_3, rad2deg(alpha_3), 'o', 'MarkerSize', 5);
        plot(v_w_4, rad2deg(alpha_4), 'square', 'MarkerSize', 7);
        hold off
        %         title('Angle of attack vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        legend(['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m'])
        ax = gca;
        ax.Box = 'on';
        xlim([3.2,14])
        
        figure()
        plot(v_w_1,l_1,v_w_2,l_2,v_w_3,l_3,v_w_4,l_4)
%         title('Tethers length vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$l$ (m)')
        ylim([19 52])
        legend({['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m']},'location','northeast')
        xlim([3.2,14])
        
        figure()
        hold on
        plot(v_w_1, l_1-l0(1), '-', 'MarkerSize', 3);
        plot(v_w_2, l_2-l0(2), '*', 'MarkerSize', 3);
        plot(v_w_3, l_3-l0(3), 'o', 'MarkerSize', 5);
        plot(v_w_4, l_4-l0(4), 'square', 'MarkerSize', 7);
        hold off
        xlabel('$v_w$ (m/s)')
        ylabel('$l-l_0$ (m)')
        legend({['$l_0 = ',num2str(l0(1),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(2),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(3),'%.2f'),'$ m'],['$l_0 = ',num2str(l0(4),'%.2f'),'$ m']},'location','best')
        xlim([3.2,14])
        
    case 4
        mS= zeros(1,4);

        for j = 1:4
            mS(j) = p.m_S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.m_S = p.m_S + 500;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:))
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        legend(['$m_S = ',num2str(mS(1),'%.2f'),'$ kg'],['$m_S = ',num2str(mS(2),'%.2f'),'$ kg'],['$m_S = ',num2str(mS(3),'%.2f'),'$ kg'],['$m_S = ',num2str(mS(4),'%.2f'),'$ kg'])
    
    case 5
        S= zeros(1,4);
        p.S = 20;
        for j = 1:4
            S(j) = p.S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w_1, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w_2, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w_3, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w_4, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.S = p.S + 20;
        end
        
        figure()
        hold on
        plot(v_w_1, x_eq_1(9,:));
        plot(v_w_2, x_eq_2(9,:));
        plot(v_w_3, x_eq_3(9,:));
        plot(v_w_4, x_eq_4(9,:));
        plot(v_w_1,zeros(1,length(v_w_1)),'k--')        
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        xlim([0 14])
        ax = gca;
        ax.Box = 'on';
        leg = legend(['$S = ',num2str(S(1),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(2),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(3),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(4),'%.2f'),'$ m$^2$'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']);
        leg.Color = 'white';
        leg.Box = 'on';
        leg.EdgeColor = 'black';
        hold off

        [v_w_1,gamma_1,alpha_1,l_1,x_eq_1] = fun_par_filter(v_w_min_1,v_w_1,gamma_1,alpha_1,l_1,x_eq_1);
        [v_w_2,gamma_2,alpha_2,l_2,x_eq_2] = fun_par_filter(v_w_min_2,v_w_2,gamma_2,alpha_2,l_2,x_eq_2);
        [v_w_3,gamma_3,alpha_3,l_3,x_eq_3] = fun_par_filter(v_w_min_3,v_w_3,gamma_3,alpha_3,l_3,x_eq_3);
        [v_w_4,gamma_4,alpha_4,l_4,x_eq_4] = fun_par_filter(v_w_min_4,v_w_4,gamma_4,alpha_4,l_4,x_eq_4);
        
        figure()
        hold on
        plot(v_w_1, x_eq_1(7,:));
        plot(v_w_2, x_eq_2(7,:));
        plot(v_w_3, x_eq_3(7,:));
        plot(v_w_4, x_eq_4(7,:));
        xlabel('$v_w$ (m/s)')
        ylabel('$x_K$ (m)')
        xlim([1.6 14])
        ax = gca;
        ax.Box = 'on';
        legend(['$S = ',num2str(S(1),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(2),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(3),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(4),'%.2f'),'$ m$^2$'])
        
        figure()
        plot(v_w_1,l_1,v_w_2,l_2,v_w_3,l_3,v_w_4,l_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$l$ (m)')
        xlim([1.6 14])
        ax = gca;
        ax.Box = 'on';
        legend(['$S = ',num2str(S(1),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(2),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(3),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(4),'%.2f'),'$ m$^2$'])

        
        figure()
        hold on
        plot(v_w_1, rad2deg(alpha_1));
        plot(v_w_2, rad2deg(alpha_2));
        plot(v_w_3, rad2deg(alpha_3));
        plot(v_w_4, rad2deg(alpha_4));
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        xlim([1.6 14])
        ax = gca;
        ax.Box = 'on';
        legend(['$S = ',num2str(S(1),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(2),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(3),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(4),'%.2f'),'$ m$^2$'])

        figure()
        plot(v_w_1,gamma_1,v_w_2,gamma_2,v_w_3,gamma_3,v_w_4,gamma_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        xlim([1.6 14])
        ax = gca;
        ax.Box = 'on';
        legend(['$S = ',num2str(S(1),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(2),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(3),'%.2f'),'$ m$^2$'],['$S = ',num2str(S(4),'%.2f'),'$ m$^2$'])
    case 6
        yA_K= zeros(1,4);
        % p.yA_K = p.b/5;
        for j = 1:4
            yA_K(j) = p.yA_K;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w_1, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w_2, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,0.7*points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w_3, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,0.5*points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w_4, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,0.3*points,x_red0_design,frame,vw_max);
            end
            p.yA_K = p.yA_K + 1;
        end

        figure()
        hold on
        plot(v_w_1, x_eq_1(9,:), '-', 'MarkerSize', 3);
        plot(v_w_2, x_eq_2(9,:), '*', 'MarkerSize', 3);
        plot(v_w_3, x_eq_3(9,:), 'o', 'MarkerSize', 5);
        plot(v_w_4, x_eq_4(9,:), 'square', 'MarkerSize', 7);
        plot(v_w_1,zeros(1,length(v_w_1)),'k--')
        plot(v_w_min_1, 0, 'pentagram', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 12);
        plot(v_w_min_2, 0, 'pentagram', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 9);
        plot(v_w_min_3, 0, 'pentagram', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 6);
        plot(v_w_min_4, 0, 'pentagram', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        xlim([0 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$y_A^K = ',num2str(yA_K(1),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(2),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(3),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off
        
        [v_w_1,gamma_1,alpha_1,l_1,x_eq_1] = fun_par_filter(v_w_min_1,v_w_1,gamma_1,alpha_1,l_1,x_eq_1);
        [v_w_2,gamma_2,alpha_2,l_2,x_eq_2] = fun_par_filter(v_w_min_2,v_w_2,gamma_2,alpha_2,l_2,x_eq_2);
        [v_w_3,gamma_3,alpha_3,l_3,x_eq_3] = fun_par_filter(v_w_min_3,v_w_3,gamma_3,alpha_3,l_3,x_eq_3);
        [v_w_4,gamma_4,alpha_4,l_4,x_eq_4] = fun_par_filter(v_w_min_4,v_w_4,gamma_4,alpha_4,l_4,x_eq_4);
        
        
        figure()
        hold on
        plot(v_w_1, x_eq_1(7,:), '-', 'MarkerSize', 3);
        plot(v_w_2, x_eq_2(7,:), '*', 'MarkerSize', 3);
        plot(v_w_3, x_eq_3(7,:), 'o', 'MarkerSize', 5);
        plot(v_w_4, x_eq_4(7,:), 'square', 'MarkerSize', 7);
        xlabel('$v_w$ (m/s)')
        ylabel('$x_K$ (m)')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$y_A^K = ',num2str(yA_K(1),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(2),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(3),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        
        figure()
        hold on
        plot(v_w_1, gamma_1, '-', 'MarkerSize', 3);
        plot(v_w_2, gamma_2, '*', 'MarkerSize', 3);
        plot(v_w_3, gamma_3, 'o', 'MarkerSize', 5);
        plot(v_w_4, gamma_4, 'square', 'MarkerSize', 7);
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$y_A^K = ',num2str(yA_K(1),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(2),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(3),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(4),'%.2f'),'$ m']},'LOCATION','best')
            
        figure()
        hold on
        plot(v_w_1, rad2deg(alpha_1), '-', 'MarkerSize', 3);
        plot(v_w_2, rad2deg(alpha_2), '*', 'MarkerSize', 3);
        plot(v_w_3, rad2deg(alpha_3), 'o', 'MarkerSize', 5);
        plot(v_w_4, rad2deg(alpha_4), 'square', 'MarkerSize', 7);
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$y_A^K = ',num2str(yA_K(1),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(2),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(3),'%.2f'),'$ m'],['$y_A^K = ',num2str(yA_K(4),'%.2f'),'$ m']},'LOCATION','best')

        
    case 7
        zA_K= zeros(1,4);
        for j = 1:4
            zA_K(j) = p.zA_K;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w_1, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w_2, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w_3, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w_4, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.zA_K = p.zA_K - 0.5;
        end

        figure()
        hold on
        plot(v_w_1, x_eq_1(9,:));
        plot(v_w_2, x_eq_2(9,:));
        plot(v_w_3, x_eq_3(9,:));
        plot(v_w_4, x_eq_4(9,:));
        plot(v_w_1,zeros(1,length(v_w_1)),'k--')        
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        leg = legend({['$z_A^K = ',num2str(zA_K(1),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(2),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(3),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','best');
        leg.Color = 'white';
        leg.Box = 'on';
        leg.EdgeColor = 'black';
        ax = gca;
        ax.Box = 'on';
        hold off
        
        [v_w_1,gamma_1,alpha_1,l_1,x_eq_1] = fun_par_filter(v_w_min_1,v_w_1,gamma_1,alpha_1,l_1,x_eq_1);
        [v_w_2,gamma_2,alpha_2,l_2,x_eq_2] = fun_par_filter(v_w_min_2,v_w_2,gamma_2,alpha_2,l_2,x_eq_2);
        [v_w_3,gamma_3,alpha_3,l_3,x_eq_3] = fun_par_filter(v_w_min_3,v_w_3,gamma_3,alpha_3,l_3,x_eq_3);
        [v_w_4,gamma_4,alpha_4,l_4,x_eq_4] = fun_par_filter(v_w_min_4,v_w_4,gamma_4,alpha_4,l_4,x_eq_4);
        
        figure()
        hold on
        plot(v_w_1, x_eq_1(7,:));
        plot(v_w_2, x_eq_2(7,:));
        plot(v_w_3, x_eq_3(7,:));
        plot(v_w_4, x_eq_4(7,:));
        xlabel('$v_w$ (m/s)')
        ylabel('$x_K$ (m)')
        legend({['$z_A^K = ',num2str(zA_K(1),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(2),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(3),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        xlim([3.2 14])
        ylim([-18 -53])
        ax = gca;
        ax.Box = 'on';
        hold off

        figure()
        plot(v_w_1,gamma_1,v_w_2,gamma_2,v_w_3,gamma_3,v_w_4,gamma_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        legend({['$z_A^K = ',num2str(zA_K(1),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(2),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(3),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        
        figure()
        plot(v_w_1,rad2deg(alpha_1),v_w_2,rad2deg(alpha_2),v_w_3,rad2deg(alpha_3),v_w_4,rad2deg(alpha_4))
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        legend({['$z_A^K = ',num2str(zA_K(1),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(2),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(3),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        
        figure()
        plot(v_w_1,l_1,v_w_2,l_2,v_w_3,l_3,v_w_4,l_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$l$ (m)')
        legend({['$z_A^K = ',num2str(zA_K(1),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(2),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(3),'%.2f'),'$ m'],['$z_A^K = ',num2str(zA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        
    case 8 
        xA_K= zeros(1,4);
        for j = 1:4
            xA_K(j) = p.xA_K;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w_1, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w_2, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w_3, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w_4, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.xA_K = p.xA_K - 0.15;
        end
        figure()
        hold on
        plot(v_w_1, x_eq_1(9,:));
        plot(v_w_2, x_eq_2(9,:));
        plot(v_w_3, x_eq_3(9,:));
        plot(v_w_4, x_eq_4(9,:));
        plot(v_w_1,zeros(1,length(v_w_1)),'k--')        
        plot(v_w_min_1, 0, 'pentagram', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 12);
        plot(v_w_min_2, 0, 'pentagram', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 9);
        plot(v_w_min_3, 0, 'pentagram', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 6);
        plot(v_w_min_4, 0, 'pentagram', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        leg = legend({['$x_A^K = ',num2str(xA_K(1),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(2),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(3),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast');
        leg.Color = 'white';
        leg.Box = 'on';
        leg.EdgeColor = 'black';
        ax = gca;
        ax.Box = 'on';
        hold off
        
        [v_w_1,gamma_1,alpha_1,l_1,x_eq_1] = fun_par_filter(v_w_min_1,v_w_1,gamma_1,alpha_1,l_1,x_eq_1);
        [v_w_2,gamma_2,alpha_2,l_2,x_eq_2] = fun_par_filter(v_w_min_2,v_w_2,gamma_2,alpha_2,l_2,x_eq_2);
        [v_w_3,gamma_3,alpha_3,l_3,x_eq_3] = fun_par_filter(v_w_min_3,v_w_3,gamma_3,alpha_3,l_3,x_eq_3);
        [v_w_4,gamma_4,alpha_4,l_4,x_eq_4] = fun_par_filter(v_w_min_4,v_w_4,gamma_4,alpha_4,l_4,x_eq_4);
        
        figure()
        hold on
        plot(v_w_1, x_eq_1(7,:));
        plot(v_w_2, x_eq_2(7,:));
        plot(v_w_3, x_eq_3(7,:));
        plot(v_w_4, x_eq_4(7,:));
        xlabel('$v_w$ (m/s)')
        ylabel('$x_K$ (m)')
        legend({['$x_A^K = ',num2str(xA_K(1),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(2),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(3),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(4),'%.2f'),'$ m']},'LOCATION','best')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        hold off

        figure()
        plot(v_w_1,gamma_1,v_w_2,gamma_2,v_w_3,gamma_3,v_w_4,gamma_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$x_A^K = ',num2str(xA_K(1),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(2),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(3),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(4),'%.2f'),'$ m']},'LOCATION','best')

        
        figure()
        plot(v_w_1,rad2deg(alpha_1),v_w_2,rad2deg(alpha_2),v_w_3,rad2deg(alpha_3),v_w_4,rad2deg(alpha_4))
        xlabel('$v_w$ (m/s)')
        ylabel('$\alpha$ (deg)')
        xlim([3.2 14])
        ylim([5 10])
        ax = gca;
        ax.Box = 'on';
        legend({['$x_A^K = ',num2str(xA_K(1),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(2),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(3),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(4),'%.2f'),'$ m']},'LOCATION','best')

        
        figure()
        plot(v_w_1,l_1,v_w_2,l_2,v_w_3,l_3,v_w_4,l_4)
        xlabel('$v_w$ (m/s)')
        ylabel('$l$ (m)')
        xlim([3.2 14])
        ax = gca;
        ax.Box = 'on';
        legend({['$x_A^K = ',num2str(xA_K(1),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(2),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(3),'%.2f'),'$ m'],['$x_A^K = ',num2str(xA_K(4),'%.2f'),'$ m']},'LOCATION','best')
    
    case 9
        w_S= zeros(1,4);
        for j = 1:4
            w_S(j) = p.w_S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.w_S = p.w_S + 1;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:),v_w,zeros(1,length(v_w)),'k--')
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        legend({['$w_S = ',num2str(w_S(1),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(2),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(3),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off

        figure()
        plot(v_w,gamma_1,v_w,gamma_2,v_w,gamma_3,v_w,gamma_4)
        title('Elevation angle vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$\gamma$ (deg)')
        legend({['$w_S = ',num2str(w_S(1),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(2),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(3),'%.2f'),'$ m'],['$w_S = ',num2str(w_S(4),'%.2f'),'$ m']},'LOCATION','northeast')
    case 10
        l_S= zeros(1,4);
        for j = 1:4
            l_S(j) = p.l_S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.l_S = p.l_S + 10;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:),v_w,zeros(1,length(v_w)),'k--')
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        legend({['$l_S = ',num2str(l_S(1),'%.2f'),'$ m'],['$l_S = ',num2str(l_S(2),'%.2f'),'$ m'],['$l_S = ',num2str(l_S(3),'%.2f'),'$ m'],['$l_S = ',num2str(l_S(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off
    case 11
        xA_S= zeros(1,4);
        p.xA_S = 0;
        for j = 1:4
            xA_S(j) = p.xA_S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.xA_S = p.xA_S - p.l_S/6;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:),v_w,zeros(1,length(v_w)),'k--')
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        legend({['$x_A^S = ',num2str(xA_S(1),'%.2f'),'$ m'],['$x_A^S = ',num2str(xA_S(2),'%.2f'),'$ m'],['$x_A^S = ',num2str(xA_S(3),'%.2f'),'$ m'],['$x_A^S = ',num2str(xA_S(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off
        
        case 12
        yA_S= zeros(1,4);
        p.yA_S = 0;
        for j = 1:4
            yA_S(j) = p.yA_S;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.yA_S = p.yA_S + p.w_S/6;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:),v_w,zeros(1,length(v_w)),'k--')
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        legend({['$y_A^S = ',num2str(yA_S(1),'%.2f'),'$ m'],['$y_A^S = ',num2str(yA_S(2),'%.2f'),'$ m'],['$y_A^S = ',num2str(yA_S(3),'%.2f'),'$ m'],['$y_A^S = ',num2str(yA_S(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off
        
    case 13
        mu_s= zeros(1,4);
        p.mu_s = 0.01;
        for j = 1:4
            mu_s(j) = p.mu_s;
            p.v_w = 10;
            switch j
                case 1
                    [x_eq_1, v_w, gamma_1, alpha_1, l_1, v_w_min_1] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 2
                    [x_eq_2, v_w, gamma_2, alpha_2, l_2, v_w_min_2] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 3
                    [x_eq_3, v_w, gamma_3, alpha_3, l_3, v_w_min_3] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
                case 4
                    [x_eq_4, v_w, gamma_4, alpha_4, l_4, v_w_min_4] = fun_par_vw(p,points,x_red0_design,frame,vw_max);
            end
            p.mu_s = p.mu_s + 0.2;
        end
        figure()
        plot(v_w,x_eq_1(9,:),v_w,x_eq_2(9,:),v_w,x_eq_3(9,:),v_w,x_eq_4(9,:),v_w,zeros(1,length(v_w)),'k--')
        title('$z_K$ vs wind velocity')
        xlabel('$v_w$ (m/s)')
        ylabel('$z_K$ (m)')
        hold on
        plot(v_w_min_1, 0, 'o', 'Color', color.blue, 'MarkerFaceColor', color.blue, 'MarkerSize', 3);
        plot(v_w_min_2, 0, 'o', 'Color', color.red, 'MarkerFaceColor', color.red, 'MarkerSize', 3);
        plot(v_w_min_3, 0, 'o', 'Color', color.green, 'MarkerFaceColor', color.green, 'MarkerSize', 3);
        plot(v_w_min_4, 0, 'o', 'Color', color.orange, 'MarkerFaceColor', color.orange, 'MarkerSize', 3);
        legend({['$\mu_s = ',num2str(mu_s(1),'%.2f'),'$ m'],['$\mu_s = ',num2str(mu_s(2),'%.2f'),'$ m'],['$\mu_s = ',num2str(mu_s(3),'%.2f'),'$ m'],['$\mu_s = ',num2str(mu_s(4),'%.2f'),'$ m'],['$z_K = 0$ m'],['$v_w|_{min} =',num2str(v_w_min_1,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_2,'%.2f'),'$ m/s'],['$v_w|_{min} = ',num2str(v_w_min_3,'%.2f'),'$ m/s'],['$v_w|_{min} =',num2str(v_w_min_4,'%.2f'),'$ m/s']},'LOCATION','northeast')
        hold off
end
p = parameters;
figures = findall(0, 'Type', 'figure');

end
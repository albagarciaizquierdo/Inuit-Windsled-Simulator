function p = parameters
% Sled variables
p.m_S   = 500;                      % Sled mass (kg)
p.l_S     = 3.5;                      % Sled length (m)
p.w_S   = 3;                        % Sled width (m)
p.Iz_S  = p.m_S*(p.l_S^2+p.w_S^2)/12; % Sled inertia in z-axis, (kg*m^2)

% Tethers (Material: Dacron)
p.l0    = 50;               % Tethers natural length (m)
p.d     = 2e-3;             % Tethers diameter (m)
p.A_t   = pi*p.d^2/4;       % Tethers cross sectional area (m^2)
p.E     = 3e9;              % Young's modulus (Dacron) (N/m^2)
p.k     = (p.E*p.A_t)/p.l0; % Stiffness constant (N/m)
p.V_t   = p.A_t*p.l0;       % Tethers volume (m^3)
p.rho_t = 1.39e3;           % Tethers density (Dacron)(kg/m^3)
p.m_t   = p.V_t*p.rho_t;    % Tethers mass (kg)

% Kite variables (Fabric material: Nylon)
p.S         = 20;                    % Kite surface (m^2)
p.c         = 2;                     % Kite chord (m)
p.b         = 5;                     % Kite span (m)
p.rho_f     = 1.14e3;                % Kite fabric density (nylon) (kg/m^3)
p.t         = 0.02e-2;               % Kite fabric thickness (m)
p.V_f       = p.S*p.t;               % Kite fabric volume (m^3)
p.m_f       = p.V_f*p.rho_f;         % Kite fabric mass (kg)
p.m_K       = p.m_f + p.m_t;         % Kite mass (kg)
p.eps       = p.b/p.c;               % Kite span to chord ratio (-)
p.Ix_K      = 1.32e-4*p.m_K*p.l0^2;  % Kite inertia in x-axis (kg·m^2)
p.Iy_K      = 2.92e-5*p.m_K*p.l0^2;  % Kite inertia in y-axis (kg·m^2)
p.Iz_K      = 1.12e-4*p.m_K*p.l0^2;  % Kite inertia in z-axis (kg·m^2)
p.Cx0       = -0.065;                % (-)
p.Cxalpha   = 0.176;                 % (-)
p.Cybeta    = -1.57;                 % (-)
p.Cz0       = 0.116;                 % (-)
p.Czalpha   = -2.97;                 % (-)
p.Clbeta    = -0.037;                % (-)
p.Clp       = -0.15;                 % (-)
p.Cm0       = 0.1332;                % (-)
p.Cmalpha   = -0.7633;               % (-)
p.Cmq       = -0.165;                % (-)
p.Cnbeta    = -0.027;                % (-)
p.Cnr       = -0.002;                % (-)

% Attachment points
p.xA_S = -p.l_S/2;    % Sled attachment points distance to OS in x_S direction (m)
p.yA_S = p.w_S/4;   % Sled attachment points distance to OS in y_S direction (m)
p.xA_K = p.b/10;    % Kite attachment points distance to OK in x_K direction (m)
p.yA_K = p.b/2;     % Kite attachment points distance to OK in y_K direction (m)
p.zA_K = 0;         % Kite attachment points distance to OK in z_K direction (m)

%
p.V_T = 7;    % Reference velocity (m/s)

% Environment
p.rho   = 1.225;    % Air density (kg/m^3)
p.v_w   = 10;       % Wind velocity (m/s)
p.mu_s  = 0.05;     % Coefficient of static friction between ice and teflon
p.mu_k  = 0.04;     % Coefficient of kinetic friction between ice and teflon
p.g     = 9.81;     % Gravity acceleration constant (m/s^2)

end
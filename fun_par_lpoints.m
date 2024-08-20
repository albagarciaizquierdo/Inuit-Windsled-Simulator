function [x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_lpoints(p,points,x_red0_design,frame)

global p
gamma = zeros(1,points);
v_w = zeros(1,points);
v_a = zeros(3,points);
alpha = zeros(1,points);
beta = zeros(1,points);
l = zeros(1,points);
x_red0 = zeros(3,points);
X_red_eq = zeros(3,points);
x_eq = zeros(18,points);

p.v_w = 10;

% Varying wind velocity
for i=1:points
    if i == 1
        x_red0(:,i) = x_red0_design;
        v_w(i) = p.v_w;
    else
        x_red0(:,i) = X_red_eq(:,i-1);
        p.v_w = p.v_w - 10/points;
        v_w(i) = p.v_w;
    end
    [X_red_eq(:,i),~,~] = my_fzero("fun_equilibrio_red",x_red0(:,i),1e-8,30,1e-6);
    x_eq(:,i) = equilibrium_conditions(X_red_eq(:,i));
    [gamma(i),~] = fun_gamma(x_eq(:,i),p);
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,alpha(i),beta(i),v_a(:,i),l(:,i),~,~] = fun_get_results(x_eq(:,i),p,frame);
end
v_w = flip(v_w);
x_eq = fliplr(x_eq);
gamma = flip(gamma);
alpha = flip(alpha);
l = flip(l);

% Finding the exact wind velocity where z_K is 0
v_w_min = 0;
z_K = x_eq(9, :);
for i = 1:(points - 1)
    if z_K(i) * z_K(i + 1) < 0  % Check if there is a sign change
        % Interpolate to find the exact v_w where z_K is min
        v_w_min = interp1(z_K([i, i + 1]), v_w([i, i + 1]), 0);
    end
end

end
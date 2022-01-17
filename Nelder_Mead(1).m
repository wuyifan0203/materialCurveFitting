% RFDM
% Nelder Mead
% Lumerical Material Database
% Data Format
% Column 1   Column 2   Column 3    Column 4              Column 5              Column 6
% epsIm      epsRe      frequency   refractive_index_Im   refractive_index_Re   wavelength  


clear all
close all
clc

data = load('Ag_Lum_Palik(0-2um).txt');

%% Read Data
lambda = transpose(data(:,6) * (1e6));  % in microns
n_real = transpose(data(:,5));
frequency = transpose(data(:,3));
eps_mea = transpose(data(:,2) - 1j * data(:,1));

lambda = lambda(end:-1:1);
eps_mea = eps_mea(end:-1:1);
n_real = n_real(end:-1:1);
frequency = frequency(end:-1:1);
%% Select The Range For Curve Fitting
i_min = min(find(lambda >= 0.2));
i_max = max(find(lambda <= 2.0));

lambda = lambda(i_min:i_max);
eps_mea = eps_mea(i_min:i_max);
n_real = n_real(i_min:i_max);
frequency = frequency(i_min:i_max);

eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*mu0);
omega_rad = 2 * pi * c0 ./ (lambda * 1e-6);
h = 6.6260755e-34; % J * s
eV_J = 1.6021766208e-19; % 1 eV =1.6021766208e-19 J
%omega_eV = omega_rad / (2*pi) * h /eV_J;
omega = omega_rad / (2*pi) * h /eV_J;

% omega = 1 ./ lambda;
eps_real_mea = real(eps_mea);
eps_imag_mea = imag(eps_mea);
%% Nelder-Mead Simplex Algorithm
Tol = 0.1;
MaxIter = 10000;

lambda_rr1 = 1.1e1;
lambda_ri1 = 1.2e1;
lambda_pr1 = 1.3e1;
lambda_pi1 = 1.4e1;

lambda_rr2 = 1.5e1;
lambda_ri2 = 1.6e1;
lambda_pr2 = 1.7e1;
lambda_pi2 = 1.8e1;

lambda_rr3 = 1.9e1;
lambda_ri3 = 2.0e1;
lambda_pr3 = 2.1e1;
lambda_pi3 = 2.2e1;

% Define a simplex of N+1 points and N parameters for optimization
% The order of the parameters: rr1, ri1, pr1, pi1, rr2, ri2, pr2, pi2, rr3, ri3, pr3, pi3
p = zeros(13, 12);
p_temp = zeros(13, 12);
I = diag([lambda_rr1, lambda_ri1, lambda_pr1, lambda_pi1, lambda_rr2, lambda_ri2,...
          lambda_pr2, lambda_pi2, lambda_rr3, lambda_ri3, lambda_pr3, lambda_pi3]);
    
p0 = rand(1, 12);
% p0 = [-4.266445917421868,...
%       -4.820807690100082,...
%       0.565972577620397,...
%       4.326478532083442,...
%       4.217529894360474,...
%       1.883722191247231,...
%       0.512865382759353,...
%       -2.079310847097983 ,...
%       -6.791269230123738 ,...
%       -55.261295529796950 ,...
%       1.958510107371007,...
%       3.728016339168878];

p(1,:) = p0;
p(2,:) = I(1,:) + p(1,:); p(3,:) = I(2,:) + p(1,:); p(4,:) = I(3,:) + p(1,:);
p(5,:) = I(4,:) + p(1,:); p(6,:) = I(5,:) + p(1,:); p(7,:) = I(6,:) + p(1,:);
p(8,:) = I(7,:) + p(1,:); p(9,:) = I(8,:) + p(1,:); p(10,:) = I(9,:) + p(1,:);
p(11,:) = I(10,:) + p(1,:);
%% Output the coefficients in rad/s
% rr1 = 5.551383614511336;
% ri1 = -4.674829278492341;
% pr1 = 0.480299752849433;
% pi1 = 3.372749338864088;
% 
% rr2 = -1.251536752771733;
% ri2 = 9.134308097402470;
% pr2 = 0.522692455226915;
% pi2 = -4.240471807705853;
% 
% rr3 = -28.834770735842053;
% ri3 = -120.4328257536638;
% pr3 = 2.915192613409641;
% pi3 = 3.518883159679624;
% 
% r1_eV = rr1 + 1j * ri1;
% p1_eV = -pr1^2 + 1j * pi1;
% 
% r2_eV = rr2 + 1j * ri2;
% p2_eV = -pr2^2 + 1j * pi2;
% 
% r3_eV = rr3 + 1j * ri3;
% p3_eV = -pr3^2 + 1j * pi3;
% 
% h = 6.6260755e-34; % J * s
% eV_J = 1.6021766208e-19; % 1 eV =1.6021766208e-19 J
% 
% r1 = 2*pi / h * eV_J * r1_eV
% r2 = 2*pi / h * eV_J * r2_eV
% r3 = 2*pi / h * eV_J * r3_eV
% 
% p1 = 2*pi / h * eV_J * p1_eV
% p2 = 2*pi / h * eV_J * p2_eV
% p3 = 2*pi / h * eV_J * p3_eV
%% Nelder Mead Initialization
n = 13; % The number of points of the simplex

nf = length(omega);
err_cf = zeros(n, 1);

counter = 1;
err = 100;
RMSE(1) = err;

p00 = p(1, :);  % The final p
while(counter <= MaxIter && err > Tol)
    %% Stage 1
    for i = 1:n
        err_cf(i) = f(p(i, :), omega, eps_mea(1, :));
    end

    [f_p, index] = sort(err_cf); % Ascending Order

    for i = 1:n
        p_temp(i, :) = p(index(i), :);
    end

    p = p_temp;

    p_s = p(1, :);
    p_l = p(end, :);
    p_nl = p(end - 1, :);

    f_s = f_p(1);
    f_l = f_p(end);
    f_nl = f_p(end - 1);

%    RMSE(counter) = f(p_s, omega, eps_mea);
    %% Stage 2
    rho = 1;
    x = 2; % x > 1
    r = 0.5; % 0 < r < 1
    sigma = 0.5;

    p_g = 1/(n-1) * (sum(p) - p_l);
    p_r = p_g + rho * (p_g - p_l);% reflection
    % Stage 2.1
    f_p_r = f(p_r, omega, eps_mea(1, :));

    if ((f_p_r >= f_s) && (f_p_r<f_nl))
        p(end, :) = p_r;
    elseif (f_p_r < f_s)
        p_e = p_g + x * (p_r - p_g); % expansion
        f_p_e = f(p_e, omega, eps_mea(1, :));
        if (f_p_e < f_p_r) % Stage 2.i.1
            p(end, :) = p_e;
%            RMSE(counter) = f(p_e, omega, eps_mea);
        elseif (f_p_e >= f_p_r)
            p(end, :) = p_r;
%            RMSE(counter) = f(p_r, omega, eps_mea);
        end
    elseif (f_p_r >= f_nl) % contraction
        if (f_p_r < f_l) % outside contraction
            p_c = p_g + r * (p_r - p_g);
        elseif (f_p_r > f_l) % inside contraction
            p_c = p_g + r * (p_l - p_g);
        end
        f_c = f(p_c, omega, eps_mea(1, :));
        if (f_c <= f_l)
            p(end, :) = p_c;
        elseif (f_c > f_l)
            p(2:end, :) = p_s + sigma * (p(2:end, :) - p_s); % shrinkage
        end
    end
    
    if (min(-imag(RFDM(p(1, :), omega))) >= 0)
        p00 = p(1, :);
        err = f(p00, omega, eps_mea);
        RMSE(counter) = err;
    else
        if (counter > 1)
            RMSE(counter) = RMSE(counter - 1);
        else
            RMSE(counter) = err;
        end
    end
    
    counter = counter + 1;
end

rmse_best = err
%% Curve Fitting Result
eps_curveFitting = RFDM(p00, omega);
%% Plot
figure(1)
hold on
plot(lambda, eps_real_mea, 'ko');
plot(lambda, real(eps_curveFitting), 'r-');
hold off
xlabel('\lambda (microns)');
ylabel('\Re[\epsilon]');
legend('Measurement','RFDM');

figure(2)
hold on
plot(lambda, -eps_imag_mea, 'ko');
plot(lambda, -imag(eps_curveFitting), 'r-');
hold off
xlabel('\lambda (microns)');
ylabel('\Im[\epsilon]');
legend('Measurement','RFDM');

figure(3)
plot(RMSE);
xlabel('The number of iterations');
ylabel('RMSE');
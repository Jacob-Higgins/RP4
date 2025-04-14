clc; clear all; close all;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLegendFontName', 'Times New Roman') % Set font for legends
set(0, 'DefaultAxesFontSize', 8)
%% Problem Parameters

throat_radius = 0.001;
exit_radius = 0.011;
nozzle_length = 0.028;

initial_exit_mach_guess = 2;

p0 = 11e5;%11
t0 = 1200;%1200, if using cea reactant temp should be 579.4k

%% Glass Material Properties

GlassProps.Type = {'Borosilicate Glass', 'Fused Silica', 'Quartz'};
GlassProps.Conductivity = [1.12, 1.38, 1.35];% W/mK
GlassProps.SpecHeatCap = [830, 750, 750];% J/kgK
GlassProps.ThermExpan = [3.25e-6, 0.55e-6, 0.45e-6];% microstrain/K
GlassProps.FracTough = [0.6e6, 0.65e6, 0.67e6];% MPa m^0.5
GlassProps.YoungsMod = [64e9,72e9,72e9]; % GPa

%% Solution Choices

Boundary_Layer_Model = "legge";


%% CEA Outputs

gas_mix_mol_mass = 1.3026e-2; %beofficefore 0.0105
%gas_mix_viscosity = 1.5075e-5; % before 4.47e-6
gas_mix_viscosity = 4.47e-6;

Universal_R = 8.314;
R_gas_mix = Universal_R/gas_mix_mol_mass;

rho0 = p0 / (R_gas_mix*t0);

exit_temperature = 317.79;
exit_pressure = 505.38;
exit_density = 2.4917e-3;
Mach_Soltn = 5.3344;
speed_of_sound = 479.73;
gamma = 1.37; %could be 1.37, 1.13, 1.4

throat_density = 0.7479;
throat_velocity = 1038.7;
throat_pressure = 5.8778e5;
exit_velocity = Mach_Soltn * speed_of_sound;

%% Boundary Layer Thickness

Reynould_Number = (exit_density * exit_velocity * nozzle_length) / gas_mix_viscosity;

boundary_thickness.legge = nozzle_length * 6.25 / (sqrt(Reynould_Number));
boundary_thickness.laminar = nozzle_length * 5 / (sqrt(Reynould_Number));
boundary_thickness.turbulent = nozzle_length * 0.37 / (nthroot(Reynould_Number,5));

Boundary_Thickness = boundary_thickness.(Boundary_Layer_Model);
%% limiting velocities, and angles

u_lim = sqrt((2*gamma)/(gamma-1) *R_gas_mix*t0);

mean_bound_layer_u_lim = 0.75 * u_lim;  %simon's paper quotes this value to be between 0.5 and 1

prandtl_mayer_limit = pi/2 *(sqrt((gamma+1)/(gamma-1))-1);
prandtl_mayer_exit = sqrt((gamma+1)/(gamma-1)) * atan(sqrt(((gamma-1)/(gamma+1))*((Mach_Soltn^2) -1))) - atan(sqrt((Mach_Soltn^2) -1));

phi_lim = rad2deg(prandtl_mayer_limit - prandtl_mayer_exit);

phi_boundary_edge = phi_lim * (1 - 2/pi* (2*Boundary_Thickness/ exit_radius)^((gamma-1)/(gamma+1)));

% Radial and angular distance (defining polar grid)
%distance = 0:0.001:50; % Radial distance
%theta_core = 0:0.1:theta_boundary_edge; % Angular range from 0 to the boundary edge

% Create the f_theta_core function (as a function of theta)
%f_theta_core_func = @(theta) (cos(pi * deg2rad(theta) / (2 * deg2rad(theta_lim)))) .^ (2 / (gamma - 1));

% Define the integrand function (compute f_theta_core inside the integrand)
%integrand = @(theta) f_theta_core_func(theta) .* sin(theta);

% Calculate the integral using MATLAB's integral function
%A_P = (throat_velocity / (2 * u_lim)) / integral(integrand, 0, deg2rad(theta_lim));

%% Constants A_P and c_rho

% Define the function f(phi)
f_core = @(phi) (cos(pi * phi / (2 * deg2rad(phi_lim)))).^(2 / (gamma - 1));

Ap_integral_core = integral(@(phi) f_core(phi) .* sin(phi), 0, deg2rad(phi_lim));

A_P_core = (throat_velocity) / (2 * u_lim * Ap_integral_core);


% c_rho_func = @(A_P_boundary) A_P_boundary * sqrt((gamma+1)/(gamma-1)) * ((2 * mean_bound_layer_u_lim) / u_lim) * (exit_radius / (2 * Boundary_Thickness))^((gamma-1)/(gamma+1));
% 
% % Define function f_boundary(phi)
% f_boundary = @(phi, A_P_boundary) (cos((pi * deg2rad(phi_boundary_edge)) / (2 * deg2rad(phi_lim))))^(2 / (gamma-1)) .* exp(-c_rho_func(A_P_boundary) * (phi - deg2rad(phi_boundary_edge)));
% 
% % Define function to solve for A_P_boundary
% Ap_to_solve = @(A_P_boundary) (throat_velocity) / (2 * u_lim * integral(@(phi) f_boundary(phi, A_P_boundary) .* sin(phi), 0, deg2rad(phi_lim))) - A_P_boundary;
% 
% % Solve for A_P_boundary using fzero
% A_P_boundary = fzero(Ap_to_solve, 2 ); %second value is an initial guess
% c_rho = A_P_boundary * sqrt((gamma+1)/(gamma-1)) * ((2* mean_bound_layer_u_lim)/(u_lim)) * (exit_radius/(2*Boundary_Thickness))^((gamma-1)/(gamma+1));


%% define flowfield grid


% Define Variables
r = logspace(-3, 2, 1000);  % Distance from 0.001 to 100
phi_core = linspace(0, deg2rad(phi_boundary_edge), 1000);  % Theta from 0 to theta_boundary_edge
phi_boundary = linspace(deg2rad(phi_boundary_edge), deg2rad(phi_lim), 1000);
% Create Meshgrid
[R1, Phi_core] = meshgrid(r, phi_core);
[R2, Phi_boundary] = meshgrid(r, phi_boundary);

%% core f(phi), density, and mach, and T, and u


% Compute f_theta_core
f_theta_core = (cos((pi * Phi_core) ./ (2 * deg2rad(phi_lim)))).^(2/(gamma-1));


% Compute Plume Density
plume_rho_core = throat_density * A_P_core .* (throat_radius ./ R1).^2 .* f_theta_core;

%Compute Plume Mach
plume_Mach_core = sqrt((2/(gamma-1)) .* ((rho0./plume_rho_core).^(gamma-1) -1));

plume_Temp_core = t0 .* (1 + ((gamma-1)/2).*plume_Mach_core.^2 ).^-1;

plume_u_core = sqrt((2*gamma*R_gas_mix*t0 .* (1 - (plume_rho_core./rho0).^(gamma-1)))./(gamma-1));


figure()
semilogx(r, plume_u_core(1,:))
xlabel("centreline distance [m]")
ylabel("velocity [m/s]")
yline(u_lim,'--')

figure()
semilogx(r, plume_Temp_core(1,:),r, plume_Mach_core(1,:))
xlabel("centreline distance [m]")
legend("temp", "Mach")
xlim([10e-2 10e0])

x = 1;
y = 801;
S = sqrt(gamma/2)*plume_Mach_core(x,y);
T_w = 300;
heat_flux = plume_rho_core(x,y) * R_gas_mix * plume_Temp_core(x,y) * sqrt((R_gas_mix*plume_Temp_core(x,y))/2*pi) * ((S^2 + (gamma/(gamma-1)) - ((gamma+1)/(2*(gamma-1)))*(T_w/plume_Temp_core(x,y))) * (exp(-(S*sin(y))^2) + sqrt(pi)*(S*sin(y))*(1 + erf((S*sin(y))))) - 0.5*exp(-(S*sin(y))^2)  );


%function handle versions of key parameters

func_rho_core = @(r, phi) throat_density * A_P_core * (throat_radius / (r+exit_radius))^2 * f_core(phi);
func_Mach_core = @(r, phi) sqrt((2 / (gamma - 1)) * ((rho0 / func_rho_core(r, phi))^(gamma - 1) - 1));
func_temp_core = @(r, phi) t0 * (1 + ((gamma - 1) / 2) * func_Mach_core(r, phi)^2 )^(-1);
func_u_core = @(r, phi) sqrt((2 * gamma * R_gas_mix * t0 * (1 - (func_rho_core(r, phi) / rho0)^(gamma - 1))) / (gamma - 1));



%% Define Boundary layer Effective Conditions, f(phi), density, Mach

c_po = (-log(exit_pressure/throat_pressure))/(deg2rad(phi_lim- phi_boundary_edge));


cu_to_solve = @(c_u) (exp(-c_u * (deg2rad(phi_lim - phi_boundary_edge))) - 1) / (-c_u) - 0.75;
c_u = -fzero(cu_to_solve, 0.1);
c_u = 0.811;
rho0_boundary = (rho0 * exp(-(c_po - 2*c_u).*( phi_boundary- deg2rad(phi_boundary_edge))));

% Define function for variable u_lim in the boundary layer
%u_lim_boundary = @(phi) u_lim .* exp(-c_u .* (phi - deg2rad(phi_boundary_edge)));

u_lim_boundary = u_lim .* exp(- c_u .*( phi_boundary -deg2rad(phi_boundary_edge)));

% Initialize arrays for A_P_boundary and c_rho
A_P_boundary = zeros(1, length(phi_boundary));
c_rho = zeros(1, length(phi_boundary));

% Loop through each boundary angle and solve for A_P_boundary and c_rho
for i = 1:length(phi_boundary)
    % Extract the current boundary angle and corresponding u_lim
    phi_current = phi_boundary(i);
    u_lim_current = u_lim_boundary(i); % Directly access precomputed values

    % Define c_rho function at this specific phi
    c_rho_func = @(A_P) A_P * sqrt((gamma+1)/(gamma-1)) .* ...
        ((2 * mean_bound_layer_u_lim) / u_lim) .* ...
        (exit_radius / (2 * Boundary_Thickness))^((gamma-1)/(gamma+1));

    % Define f_boundary at this phi
    f_boundary = @(phi, A_P) (cos((pi * phi) / ...
        (2 * deg2rad(phi_lim)))).^(2 / (gamma-1)); ... 
        %.* exp(-c_rho_func(A_P) .* (phi - deg2rad(phi_boundary_edge)));

    % Define function to solve for A_P_boundary at this phi
    Ap_to_solve = @(A_P) (throat_velocity) ./ ...
        (2 * u_lim_current * ...
        integral(@(phi) f_boundary(phi, A_P) .* sin(phi), 0, deg2rad(phi_lim))) ...
        - A_P;

    % Solve for A_P_boundary at this specific phi
    if i == 1
        A_P_initial_guess = A_P_core; % Ensure first value is A_P_core
    else
        A_P_initial_guess = A_P_boundary(i-1); % Use previous value as initial guess
    end

    A_P_boundary(i) = fzero(Ap_to_solve, A_P_initial_guess);

    % Compute c_rho for this specific phi
    c_rho(i) = c_rho_func(A_P_core);
end


f_theta_boundary = (cos((pi * deg2rad(phi_boundary_edge))/ (2*deg2rad(phi_lim)))^(2/(gamma-1))) .* exp(-c_rho .* (Phi_boundary-deg2rad(phi_boundary_edge)));
plume_rho_boundary = throat_density * A_P_boundary .* (throat_radius ./ R2).^2 .* f_theta_boundary;
plume_Mach_boundary = sqrt((2/(gamma-1)) .* ((rho0_boundary./plume_rho_boundary).^(gamma-1) -1));

func_f_theta_boundary = @(phi) (cos((pi * deg2rad(phi_boundary_edge))/ (2*deg2rad(phi_lim)))^(2/(gamma-1))) * exp(-c_rho(1) * (phi-deg2rad(phi_boundary_edge)));
func_rho0_boundary = @(phi) rho0* exp(-(c_po- 2*c_u)*(phi-deg2rad(phi_boundary_edge)));
func_temp0_boundary = @(phi) t0* exp(-2*c_u*(phi-deg2rad(phi_boundary_edge)));
func_ulim_boundary = @(phi) u_lim * exp(-c_u *(phi-deg2rad(phi_boundary_edge)));

func_rho_boundary = @(r, phi) throat_density * A_P_core * (throat_radius / (r+exit_radius))^2 * func_f_theta_boundary(phi);
func_Mach_boundary = @(r, phi) sqrt((2 / (gamma - 1)) * ((func_rho0_boundary(phi) / func_rho_boundary(r, phi))^(gamma - 1) - 1));
func_temp_boundary = @(r, phi) func_temp0_boundary(phi) * (1 + ((gamma - 1) / 2) * func_Mach_boundary(r, phi)^2 )^(-1);
func_u_boundary = @(r, phi) sqrt((2 * gamma * R_gas_mix * func_temp0_boundary(phi) * (1 - (func_rho_boundary(r, phi) / func_rho0_boundary(phi))^(gamma - 1))) / (gamma - 1));


%% Bird's Parameter

particle_mass = 4.65e-26;
particle_sigma = 2.3e-16;

lambda = particle_mass./ (sqrt(2) .* plume_rho_core .* (particle_sigma.*(plume_Temp_core./298).^0.75));

v_local = sqrt((8*1.38e-23.*plume_Temp_core)./(pi * particle_mass));
plume_Bird_core = (2.* plume_u_core)./( (v_local ./lambda).*R1);

figure()
loglog(r, plume_Bird_core(1,:),r, plume_Bird_core(500,:),r, plume_Bird_core(1000,:))
yline(0.02)
yline(2)
legend(num2str(rad2deg(phi_core(1))),num2str(rad2deg(phi_core(500))), num2str(rad2deg(phi_core(1000))))


%defining function versions of lambda, local v, and BIrd

func_lambda_core = @(r, phi) particle_mass ./ (sqrt(2) .* func_rho_core(r, phi) .* (particle_sigma .* (func_temp_core(r, phi) ./ 298).^0.75));

func_v_local_core = @(r,phi) sqrt((3 * 1.38e-23 .* func_temp_core(r, phi)) ./ ( particle_mass));

func_nu_core= @(r, phi) func_v_local_core(r,phi) / func_lambda_core(r,phi);

func_bird_core = @(r, phi) 2* func_u_core(r,phi) / (func_nu_core(r,phi) * r);

func_Bird_core = @(r, phi) sqrt(pi / (2 * R_gas_mix * func_temp_core(r, phi))) * (func_u_core(r, phi) * func_lambda_core(r,phi))/r;

func_lambda_boundary = @(r, phi) particle_mass ./ (sqrt(2) .* func_rho_boundary(r, phi) .* (particle_sigma .* (func_temp_boundary(r, phi) ./ 298).^0.75));

func_v_local_boundary = @(r,phi) sqrt((8 * 1.38e-23 .* func_temp_boundary(r, phi)) ./ (pi * particle_mass));

func_nu_boundary = @(r, phi) func_v_local_boundary(r,phi) / func_lambda_boundary(r,phi);

func_bird_boundary = @(r, phi) 2* func_u_boundary(r,phi) / (func_nu_boundary(r,phi) * r);

func_Bird_boundary = @(r, phi) sqrt(pi / (2 * R_gas_mix * func_temp_boundary(r, phi))) * (func_u_boundary(r, phi) * func_lambda_boundary(r,phi))/r;


%% Determining Mach, density, and Bird contours

% Convert Polar Coordinates to Cartesian for theta_boundary_edge
[X_boundary, Y_boundary] = pol2cart(ones(1,length(r))* deg2rad(phi_boundary_edge), r);

% Convert Polar Coordinates to Cartesian for theta_lim
[X_lim, Y_lim] = pol2cart(ones(1,length(r))*deg2rad(phi_lim), r);


% Compute density at different r for each of the given density values
densities = [1e-5, 1e-6, 1e-7, 1e-8];  % Given density values (log10 scale)
%Compute density at different r for each of the given density values
machs = [40,55, 65, 70,80 ];  % Given density values (log10 scale)

birds = [2];

% Compute corresponding r for each density
r_densities_core = zeros(length(densities), length(phi_core));

for i = 1:length(densities)
    for j = 1:length(phi_core)
        % Solve for r based on the density formula
        % Re-arranging: r = throat_radius * (throat_density * A_P * f_theta_core / density)^(1/2)
        r_densities_core(i, j) = throat_radius * (throat_density * A_P_core * f_theta_core(j) / densities(i))^(1/2);
    end
end

r_densities_boundary = zeros(length(densities), length(phi_core));

for i = 1:length(densities)
    for j = 1:length(phi_boundary)
        % Solve for r based on the density formula
        % Re-arranging: r = throat_radius * (throat_density * A_P * f_theta_core / density)^(1/2)
        r_densities_boundary(i, j) = throat_radius * (throat_density * A_P_core * f_theta_boundary(j) / densities(i))^(1/2);
    end
end


% Compute corresponding r for each density
r_machs_core = zeros(length(machs), length(phi_core));


for i = 1:length(machs)
    for j = 1:length(phi_core)
        % Solve for r based on the density formula
        % Re-arranging: r = throat_radius * (throat_density * A_P * f_theta_core / density)^(1/2)
        r_machs_core(i, j) = throat_radius * (throat_density * A_P_core * f_theta_core(j) / (rho0/ ((1+ ((gamma-1)/2)*machs(i)^2)^(1/(gamma-1)))))^(1/2);
    end
end

% Compute corresponding r for each density
r_machs_boundary = zeros(length(machs), length(phi_core));



for i = 1:length(machs)
    for j = 1:length(phi_boundary)
        % Solve for r based on the density formula
        % Re-arranging: r = throat_radius * (throat_density * A_P * f_theta_core / density)^(1/2)
        r_machs_boundary(i, j) = throat_radius * (throat_density * A_P_core * f_theta_boundary(j) / (rho0_boundary(j)/ ((1+ ((gamma-1)/2)*machs(i)^2)^(1/(gamma-1)))))^(1/2);
    end
end


r_birds_core = zeros(length(birds), length(phi_core));

for i = 1:length(birds)
    for j = 1:length(phi_core)
        % Set an initial guess for r
        r_guess = 1;  % This is just a starting guess for r
        
        % Solve for r iteratively
        tol = 1e-9;  % Set tolerance for convergence
        max_iter = 1000;  % Set maximum number of iterations
        iter = 0;
        
        while iter < max_iter
            % Compute u_core, lambda, and v_local at the current guess of r
            %u_core_val = func_u_core(r_guess, phi_core(j));  % Note the correct function for u_core
            %lambda_core_val = func_lambda_core(r_guess, phi_core(j));
            %v_local_core_val = func_v_local_core(r_guess, phi_core(j));
            
            % Compute Birds(i) based on the guess of r
            %birds_computed_core = (2 * u_core_val * lambda_core_val) / (r_guess * v_local_core_val);
            birds_core = func_Bird_core(r_guess,phi_core(j));
            % Log values for debugging
            %disp(['i = ', num2str(i), ', j = ', num2str(j), ', r_guess = ', num2str(r_guess)]);
            %disp(['u_core = ', num2str(u_core_val), ', lambda = ', num2str(lambda_val), ', v_local = ', num2str(v_local_val), ', birds_computed = ', num2str(birds_computed)]);
            
            % Compute the error between computed Birds(i) and the actual Birds(i)
            error = birds_core - birds(i);
            
            % If the error is within tolerance, break out of the loop
            if abs(error) < tol
                break;
            end
            
            % Update r_guess using a simple update rule (e.g., gradient descent)
            r_guess = r_guess - error * 0.1;  % Update rule (simple, could be more sophisticated)
            
            iter = iter + 1;
        end
        
        % Store the solved r value
        r_birds_core(i, j) = r_guess;
    end
end

r_birds_boundary = zeros(length(birds), length(phi_boundary));

for i = 1:length(birds)
    for j = 1:length(phi_boundary)
        % Set an initial guess for r
        r_guess = 1;  % This is just a starting guess for r
        
        % Solve for r iteratively
        tol = 1e-9;  % Set tolerance for convergence
        max_iter = 1000;  % Set maximum number of iterations
        iter = 0;
        
        while iter < max_iter
            % Compute u_core, lambda, and v_local at the current guess of r
            %u_boundary_val = func_u_boundary(r_guess, phi_boundary(j));  % Note the correct function for u_core
            %lambda_boundary_val = func_lambda_boundary(r_guess, phi_boundary(j));
            %v_local_boundary_val = func_v_local_boundary(r_guess, phi_boundary(j));
            
            % Compute Birds(i) based on the guess of r
            %birds_computed_boundary = (2 * u_boundary_val * lambda_boundary_val) / (r_guess * v_local_boundary_val);
            birds_boundary = func_Bird_boundary(r_guess,phi_boundary(j));
            

            % Log values for debugging
            %disp(['i = ', num2str(i), ', j = ', num2str(j), ', r_guess = ', num2str(r_guess)]);
            %disp(['u_core = ', num2str(u_core_val), ', lambda = ', num2str(lambda_val), ', v_local = ', num2str(v_local_val), ', birds_computed = ', num2str(birds_computed)]);
            
            % Compute the error between computed Birds(i) and the actual Birds(i)
            error = birds_boundary - birds(i);
            
            % If the error is within tolerance, break out of the loop
            if abs(error) < tol
                break;
            end
            
            % Update r_guess using a simple update rule (e.g., gradient descent)
            r_guess = r_guess - error * 0.1;  % Update rule (simple, could be more sophisticated)
            
            iter = iter + 1;
        end
        
        % Store the solved r value
        r_birds_boundary(i, j) = r_guess;
    end
end

%% PLUME CONTOUR GRAPH


% --- Plotting the two lines ---
figure();
hold on;

% Plot the line corresponding to theta_boundary_edge
plot(X_boundary, Y_boundary, 'k-.', 'LineWidth', 1);
text(8.2, 4, ['\phi_0 = ' num2str(phi_boundary_edge, '%.1f') '°'], ...
'Color', 'k', 'FontSize', 8, 'VerticalAlignment', 'bottom');

% Plot the line corresponding to theta_lim
plot(X_lim, Y_lim, 'k-.', 'LineWidth', 1);
text(3.2, 4, ['\phi_{lim} = ' num2str(phi_lim, '%.1f') '°'], ...
'Color', 'k', 'FontSize', 8, 'VerticalAlignment', 'bottom');



% Labels and Formatting
xlabel('x [m]');
ylabel('y [m]');
axis equal;  % Equal scaling for both axes
grid on;  % Turn on grid

xlim([0 12]);
ylim([0 4.5]);


colors = lines(length(densities)+ length(machs)); % Get a colormap with `length(densities)` colors
colororder(colors); % Set the same color order for the figure

% Plot lines of constant density for core
for i = 1:length(densities)
    [X1, Y1] = pol2cart(phi_core, r_densities_core(i, :));
    plot(X1, Y1, 'b');
    [X1_Text, Y1_Text] = adjustTextPosition(X1(1), Y1(1), X1, Y1, xlim, ylim);
    text(X1_Text, Y1_Text, sprintf('10^{%d}', floor(log10(densities(i)))), ...
    'Color', 'b', 'FontSize', 7, 'VerticalAlignment', 'bottom');
end

% Plot lines of constant density for boundary
for i = 1:length(densities)
    [X2, Y2] = pol2cart(phi_boundary, r_densities_boundary(i, :));
    plot(X2, Y2, 'b');
end

% Plot lines of constant mach for core
for i = 1:length(machs)
    [X3, Y3] = pol2cart(phi_core, r_machs_core(i, :));
    plot(X3, Y3, 'r');
    [X3_Text, Y3_Text, vert_alignment] = adjustTextPosition(X3(1), Y3(1), X3, Y3, xlim, ylim);
    text(X3_Text, Y3_Text, num2str(machs(i)), ...
    'Color', 'r', 'FontSize', 7, 'VerticalAlignment', vert_alignment);
end

% Plot lines of constant mach for boundary
for i = 1:length(machs)
    [X4, Y4] = pol2cart(phi_boundary, r_machs_boundary(i, :));
   plot(X4, Y4,'r');
end

for i = 1:length(birds)
    [X5, Y5] = pol2cart( phi_core, r_birds_core(i, :));
   plot(X5, Y5,'g');
   [X5_Text, Y5_Text, vert_alignment] = adjustTextPosition(X5(1), Y5(1), X5, Y5, xlim, ylim);
    text(X5_Text, Y5_Text, num2str(birds(i)), ...
    'Color', 'G', 'FontSize', 7, 'VerticalAlignment', vert_alignment);
end

for i = 1:length(birds)
    [X6, Y6] = pol2cart( phi_boundary, r_birds_boundary(i, :));
   plot(X6, Y6,'g');
end



%for legend entries
h_rho = plot(NaN,NaN, 'b', 'DisplayName', '\rho','LineWidth', 1.5);
h_M = plot(NaN,NaN, 'r', 'DisplayName', 'M','LineWidth', 1.5);
h_P = plot(NaN,NaN, 'g', 'DisplayName', 'P','LineWidth', 1.5);
legend([h_rho, h_M, h_P], {'\rho', 'M', 'P'}, 'Location', 'northwest', 'FontSize', 8)

hold off


exportgraphics(gcf, "Report/Figures/MainSimonsPlot.pdf", 'ContentType','vector')


% figure(2)
% loglog(r, plume_rho_core(1,:))
% grid on
% xlabel("Distance along centre steamline [m]")
% ylabel("Density [kg/m^3]")



%% Mean free path & Knudsen


% Assuming 'hatchfill' is properly installed and on the path
r_Knudsen = logspace(-2,3,1000);
debris_characteristic_length_1 = 2;
debris_characteristic_length_2 = 4;

N2_particle_mass = 4.65e-26;
N2_sigma = 1.6e-20;
NH3_particle_mass = 2.83e-26;
NH3_sigma = 2e-20;

% Compute N2 and NH3 lambda and Knudsen numbers for both characteristic lengths
N2_lambda_1 = N2_particle_mass./(sqrt(2) * N2_sigma .* plume_rho_core./100);
N2_lambda_2 = N2_particle_mass./(sqrt(2) * N2_sigma .* plume_rho_core./100);
N2_Knudsen_1 = N2_lambda_1 ./ debris_characteristic_length_1;
N2_Knudsen_2 = N2_lambda_2 ./ debris_characteristic_length_2;

NH3_lambda_1 = NH3_particle_mass./(sqrt(2) * NH3_sigma .* plume_rho_core./100);
NH3_lambda_2 = NH3_particle_mass./(sqrt(2) * NH3_sigma .* plume_rho_core./100);
NH3_Knudsen_1 = NH3_lambda_1 ./ debris_characteristic_length_1;
NH3_Knudsen_2 = NH3_lambda_2 ./ debris_characteristic_length_2;

% Plotting the Knudsen number vs distance
figure()
hold on  % Ensure all elements are plotted on the same figure

% Define x-axis range for the fill regions
x_fill = [min(r_Knudsen), max(r_Knudsen), max(r_Knudsen), min(r_Knudsen)];

% Green Region: Below y = 0.01
y_fill_green = [1e-8, 1e-8, 0.01, 0.01]; % Ensure it covers the bottom
fill(x_fill, y_fill_green, [0, 1, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Yellow Region: Between y = 0.01 and y = 10
y_fill_yellow = [0.01, 0.01, 10, 10];
fill(x_fill, y_fill_yellow, [1, 1, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Set log-log scaling explicitly
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize',10)

% N2 plot for debris_characteristic_length = 2 and 4
h1 = loglog(r_Knudsen, N2_Knudsen_1(1,:), 'Color', [0.7, 0, 0], 'LineWidth', 1.5);  % Dark red for L=2
h2 = loglog(r_Knudsen, N2_Knudsen_2(1,:), 'Color', [1, 0.5, 0.5], 'LineWidth', 1.5);  % Light red for L=4

% NH3 plot for debris_characteristic_length = 2 and 4
h3 = loglog(r_Knudsen, NH3_Knudsen_1(1,:), 'Color', [0, 0, 0.7], 'LineWidth', 1.5);  % Dark blue for L=2
h4 = loglog(r_Knudsen, NH3_Knudsen_2(1,:), 'Color', [0.5, 0.7, 1], 'LineWidth', 1.5);  % Light blue for L=4

% Create patches for the region to hatch
% N2 patch
X_N2 = [r_Knudsen, flip(r_Knudsen)];
Y_N2 = [N2_Knudsen_1(1,:), flip(N2_Knudsen_2(1,:))];
patch_N2 = patch(X_N2, Y_N2, 'r', 'FaceAlpha', 0, 'EdgeColor', 'none');

% NH3 patch
X_NH3 = [r_Knudsen, flip(r_Knudsen)];
Y_NH3 = [NH3_Knudsen_1(1,:), flip(NH3_Knudsen_2(1,:))];
patch_NH3 = patch(X_NH3, Y_NH3, 'b', 'FaceAlpha', 0, 'EdgeColor', 'none');

% Apply hatchfill on the created patches
hatch_N2 = hatchfill(patch_N2, 'cross', 45, 2, 'r');  % N2 shaded region
hatch_NH3 = hatchfill(patch_NH3, 'cross', 45, 2, 'b');  % NH3 shaded region

% Set hatch color explicitly (in case it isn't applied)
set(hatch_N2, 'Color', 'r');  % Set N2 hatch color to red
set(hatch_NH3, 'Color', 'b');  % Set NH3 hatch color to green

% Adding horizontal lines for transition onset and free molecular onset
yline(0.01, "Label", "Kn=0.01", 'LineWidth', 1, 'Color', 'k', 'LineStyle','--','FontName','Times New Roman')
yline(10, "Label", "Kn=10", 'LineWidth', 1, 'Color', 'k', 'LineStyle','--','FontName','Times New Roman')

% Labels and titles
xlabel("Distance r along plume centreline [m]")
ylabel("Knudsen Number [-]")

ylim([1e-3 1e2])
xlim([1e-2 1e2])

% Add a legend to the plot
legend([h1, h2, h3, h4],'N_2, L=2', 'N_2, L=4', 'NH_3, L=2', 'NH_3, L=4 ', 'Location', 'northwest')

hold off
exportgraphics(gcf, "Report/Figures/KnudsenRange.pdf", 'ContentType','vector')



%% Knudsen Number

debris_characteristic_length = 10;

N2_particle_mass = 5.32e-26;
N2_sigma = 2.83e-19;
N2_lambda = N2_particle_mass./ (sqrt(2) * N2_sigma .* plume_rho_core);
N2_Knudsen.Throat_R = N2_lambda/ throat_radius;
N2_Knudsen.Throat_D = N2_lambda/ (throat_radius*2);
N2_Knudsen.Exit_R = N2_lambda/exit_radius;
N2_Knudsen.Exit_D = N2_lambda/(exit_radius*2);
N2_Knudsen.Nozz_L = N2_lambda/nozzle_length;
N2_Knudsen.Debris = N2_lambda/debris_characteristic_length;


%% Rarefaction effects on plume_u and plume_T

%introduce uLim and T_f at P=2

%% Finalise Plume Struct
%%%ADD Knudsen number

plume.rho = @(r,phi) (phi <= phi_boundary_edge) .* func_rho_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_rho_boundary(r,phi);
plume.fPhi = @(phi) (phi <= phi_boundary_edge) .* f_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_f_theta_boundary(r,phi);
plume.u = @(r,phi) (phi <= phi_boundary_edge) .* func_u_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_u_boundary(r,phi);
plume.T = @(r,phi) (phi <= phi_boundary_edge) .* func_temp_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_temp_boundary(r,phi);
plume.M = @(r,phi) (phi <= phi_boundary_edge) .* func_Mach_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_Mach_boundary(r,phi);
plume.uLim = @(phi) (phi <= phi_boundary_edge) .* u_lim + (phi > phi_boundary_edge & phi <= phi_lim) .* func_ulim_boundary(phi);
plume.rho0 = @(phi) (phi <= phi_boundary_edge) .* rho0 + (phi > phi_boundary_edge & phi <= phi_lim) .* func_rho0_boundary(r,phi);
plume.T0 = @(phi) (phi <= phi_boundary_edge) .* t0 + (phi > phi_boundary_edge & phi <= phi_lim) .* func_temp0_boundary(r,phi);
plume.lambda = @(r,phi) (phi <= phi_boundary_edge) .* func_lambda_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_lambda_boundary(r,phi);
plume.vLocal = @(r,phi) (phi <= phi_boundary_edge) .* func_v_local_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_v_local_boundary(r,phi);
plume.nu = @(r,phi) (phi <= phi_boundary_edge) .* func_nu_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_nu_boundary(r,phi);
plume.Bird = @(r,phi) (phi <= phi_boundary_edge) .* func_Bird_core(r,phi) + (phi > phi_boundary_edge & phi <= phi_lim) .* func_Bird_boundary(r,phi);
plume.phi0 = phi_boundary_edge;
plume.phiLim = phi_lim;
plume.ReynouldsExit = Reynould_Number;
plume.BoundaryThicknessExit = Boundary_Thickness;
plume.uRarefied = @(r,phi) (plume.Bird(r,phi) <= 2) .* plume.u(r,phi) +  (plume.Bird(r,phi) > 2) .* plume.uLim(phi);
plume.TRarefied = @(r,phi) (plume.Bird(r,phi) <= 2) .* plume.T(r,phi) +  (plume.Bird(r,phi) > 2) .* plume.T(getRf(phi, plume),phi);
plume.SpeedRatio = @(r,phi) plume.uRarefied(r,phi) / sqrt(2* R_gas_mix* plume.TRarefied(r,phi));

%% Plotting Centreline T and u with rarefaction effect at P=2 (could also show angular decay of T and u at r=10)
uRarefied_Centreline = zeros(1, length(r));
for i=1:length(r)
    uRarefied_Centreline(i) = plume.uRarefied(r(i),0);

end
TRarefied_Centreline = zeros(1, length(r));
for i=1:length(r)
    TRarefied_Centreline(i) = plume.TRarefied(r(i),0);

end

figure()
% Normalize data
u_norm = uRarefied_Centreline(:) ./ 2381.6;
T_norm = TRarefied_Centreline(:) ./ 1.7673;

% Create stacked plot
semilogx(r, T_norm, 'r');
xlabel("Distance r along plume centreline [m]")
ylabel("$T / T_{f}$", 'Interpreter', 'latex')
grid on;
exportgraphics(gcf, "Report/Figures/NormalisedT.pdf", 'ContentType', 'vector')

figure()

% Create stacked plot
semilogx(r, u_norm, 'b');
xlabel("Distance r along plume centreline [m]")
ylabel("$u / u_{lim}$", 'Interpreter', 'latex')
grid on;
exportgraphics(gcf, "Report/Figures/NormalisedU.pdf", 'ContentType', 'vector')



%% Retrieve Surface Nodal Locations

debris_distance = 10;
Incidence = 20;
surface_width = 1;
surface_height = 2;
number_of_nodes_width = 3;
number_of_nodes_height = 5;
generate_Node_Plot = true;

nodes = SurfacePosition(debris_distance, Incidence, surface_width, surface_height, number_of_nodes_width, number_of_nodes_height, generate_Node_Plot);
%% Compute Thermal Network

N = length(nodes);
surf_temp = 300;
Initial_surf_temp = surf_temp * ones(N,1);
% Define material and environment properties
conductivity = 1.5;  % Thermal conductivity (W/m·K) (example: aluminum)
rho = 2700; % Density (kg/m³)
cp = 712.8; % Specific heat capacity (J/kg·K)
epsilon = 0.8; % Emissivity (assumed for the panel)
sigma = 5.67e-8; % Stefan-Boltzmann constant (W/m²K⁴)
dt = 1; % Time step (seconds)
T_f = 1.8;%was 26.63
Incidence = deg2rad(Incidence);
sigma_e = 1;
speed_ratio = u_lim./sqrt(2*R_gas_mix.*T_f);
q_plume = sigma_e*R_gas_mix .* T_f .* func_rho_core(debris_distance,0) .* sqrt(R_gas_mix .* T_f ./ (2 * pi)) ...
    .* ( (speed_ratio.^2 + (gamma / (gamma - 1))) - ((gamma + 1) / (2 * (gamma - 1))) .* (surf_temp ./ T_f) ) ...
    .* ( exp(-(speed_ratio .* sin(Incidence)).^2) + sqrt(pi) .* (speed_ratio .* sin(Incidence)) .* (1 + erf(speed_ratio .* sin(Incidence))) ) ...
    - (1/2) .* exp(-(speed_ratio .* sin(Incidence)).^2);

% Calculate node spacing
dx = surface_width / (number_of_nodes_width - 1);
dy = surface_height / (number_of_nodes_height - 1);
A = dx * dy; % Area of a node
surface_thickness = [0.00003, 0.001, 0.01, 0.02, 0.05];


%% q with theta
wall_temp_initial = [220,300, 380];
for j=1:length(wall_temp_initial)
    T_i = wall_temp_initial(j);
    for i=1:1000
        theta = pi*i/1000;
        q_plume_theta(i,j) = sigma_e*R_gas_mix .* T_f .* func_rho_core(debris_distance,0) .* sqrt(R_gas_mix .* T_f ./ (2 * pi)) ...
        .* ( (speed_ratio.^2 + (gamma / (gamma - 1))) - ((gamma + 1) / (2 * (gamma - 1))) .* (T_i ./ T_f) ) ...
        .* ( exp(-(speed_ratio .* sin(theta)).^2) + sqrt(pi) .* (speed_ratio .* sin(theta)) .* (1 + erf(speed_ratio .* sin(theta))) ) ...
        - (1/2) .* exp(-(speed_ratio .* sin(theta)).^2);
    end

end

AoA = linspace(pi/1000,pi,1000);

figure()
for j=1:length(wall_temp_initial)
    plot(rad2deg(AoA), q_plume_theta(:,j))
    hold on
end
grid on
xlabel('Angle of Attack \theta [\circ]')
% Modify ylabel to reflect that it's heat flux, q_dot
ylabel('Heat flux [W/m^2]')
legend('T_w = 220K', 'T_w = 300K','T_w = 380K')
% Set x-axis limits in degrees (0 to 180°)
xlim([0 180])
hold off
exportgraphics(gcf, "Report/Figures/HeatFluxWithTw.pdf", 'ContentType','vector')


%% 


% Define time iterations
time_steps = 1800; 


figure();
hold on;

% Initialize Plot Handles
h = gobjects(length(surface_thickness), 1); % Preallocate graphics objects
legend_labels = cell(length(surface_thickness), 1); % Store legend labels

xlabel('Time [s]');
ylabel('Peak Temperature [K]');
grid on;
xlim([1, time_steps]);
ylim([surf_temp, 315]); % Adjust as necessary




% Loop over surface thickness values
for k = 1:length(surface_thickness)
    V = dx * dy * surface_thickness(k);
    Ci = rho * cp * V; % Heat capacity per node
    
    % Initialize the plot for this surface thickness
    h(k) = plot(nan, nan, 'LineWidth', 2); % Start with empty plot
    legend_labels{k} = sprintf('%.2f mm', surface_thickness(k)*1000); % Label for legend
    
    avg_temp_history = zeros(1, time_steps); % Store average temperatures
    avg_temp_history_degC = zeros(1, time_steps); % Store average temperatures

    % Loop over time steps
    for j = 1:time_steps
        Q_conduction = zeros(N, 1);
        Q_radiation = zeros(N, 1);
        Q_input = q_plume * A .* ones(N, 1); % Plume heating applied to all nodes

        for i = 1:N
            m = nodes(i).m; % Width index
            n = nodes(i).n; % Height index
            T_i = Initial_surf_temp(i);
            q_plume = sigma_e*R_gas_mix .* T_f .* func_rho_core(debris_distance,0) .* sqrt(R_gas_mix .* T_f ./ (2 * pi)) ...
            .* ( (speed_ratio.^2 + (gamma / (gamma - 1))) - ((gamma + 1) / (2 * (gamma - 1))) .* (T_i ./ T_f) ) ...
            .* ( exp(-(speed_ratio .* sin(Incidence)).^2) + sqrt(pi) .* (speed_ratio .* sin(Incidence)) .* (1 + erf(speed_ratio .* sin(Incidence))) ) ...
            - (1/2) .* exp(-(speed_ratio .* sin(Incidence)).^2);
            disp(q_plume)
            Q_input(i) = q_plume * A;

            % Conductive heat exchange with neighbors
            if m > 1  % Left neighbor
                idx_left = find([nodes.m] == m-1 & [nodes.n] == n);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dx * (Initial_surf_temp(idx_left) - T_i);
            end
            if m < number_of_nodes_width  % Right neighbor
                idx_right = find([nodes.m] == m+1 & [nodes.n] == n);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dx * (Initial_surf_temp(idx_right) - T_i);
            end
            if n > 1  % Bottom neighbor
                idx_bottom = find([nodes.m] == m & [nodes.n] == n-1);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dy * (Initial_surf_temp(idx_bottom) - T_i);
            end
            if n < number_of_nodes_height  % Top neighbor
                idx_top = find([nodes.m] == m & [nodes.n] == n+1);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dy * (Initial_surf_temp(idx_top) - T_i);
            end

            % Radiative heat loss
            Q_radiation(i) = epsilon * sigma *2* A * (T_i^4 - 300^4);
        end

        % Compute temperature update
        dT = dt * (Q_conduction - Q_radiation + Q_input) / Ci;
        Initial_surf_temp = Initial_surf_temp + dT; % Update node temperatures

        % Store the average temperature
        avg_temp_history(j) = max(Initial_surf_temp);
        disp(avg_temp_history(j))

        % Update the plot with new data
        set(h(k), 'XData', 1:j, 'YData', (avg_temp_history(1:j)));

        pause(1e-5); % Pause for visualization
    end

    % Reset initial temperature after each surface thickness iteration
    Initial_surf_temp = surf_temp * ones(N, 1);
end

% Add the legend after the loop
legend(h, legend_labels);

exportgraphics(gcf, "Report/Figures/TempIncreaseWithThicknessPhi20.pdf", 'ContentType','vector')

%% Scenario
dt=0.1;
active_range = [20 160];


secsToStop = 1800 / ((active_range(2)-active_range(1))/360) ;
time_steps=round(secsToStop/dt,0);

%% 

SpinRate = 5;
SpinDecay = SpinRate/1800;
Incidence = pi/2;
SurfTemp = surf_temp * ones(N,1);

V = dx * dy * 0.001;
Ci = rho * cp * V; % Heat capacity per node


peak_temp_history = zeros(1, time_steps); % Store average temperatures
Incidence_wrapped = Incidence* ones(1, time_steps); % Store average temperatures
Rate = zeros(1, time_steps);


 % Loop over time steps
    for j = 1:time_steps
        nodes = SurfacePosition(debris_distance, rad2deg(Incidence), surface_width, surface_height, number_of_nodes_width, number_of_nodes_height, false);
        Q_conduction = zeros(N, 1);
        Q_radiation = zeros(N, 1);
        %Q_input = q_plume * A .* ones(N, 1); % Plume heating applied to all nodes
        %SpinRate = SpinRate - SpinDecay*dt;

        Incidence = Incidence + deg2rad(SpinRate*dt);
        Incidence_wrapped(j) = mod(Incidence, 2*pi);
        if Incidence_wrapped(j) >= deg2rad(active_range(1)) && Incidence_wrapped(j) <= deg2rad(active_range(2))
            SpinRate = SpinRate - SpinDecay*dt;
        else
            SpinRate = SpinRate;
        end
        Rate(j) = SpinRate;
        %disp(rad2deg(Incidence))
        fprintf('%d/%d\n', j, time_steps);

        for i = 1:N
            m = nodes(i).m; % Width index
            n = nodes(i).n; % Height index
            T_i = SurfTemp(i);
            q_plume = sigma_e*R_gas_mix .* T_f .* func_rho_core(debris_distance,0) .* sqrt(R_gas_mix .* T_f ./ (2 * pi)) ...
            .* ( (speed_ratio.^2 + (gamma / (gamma - 1))) - ((gamma + 1) / (2 * (gamma - 1))) .* (T_i ./ T_f) ) ...
            .* ( exp(-(speed_ratio .* abs(sin(Incidence))).^2) + sqrt(pi) .* (speed_ratio .* abs(sin(Incidence))) .* (1 + erf(speed_ratio .* abs(sin(Incidence)))) ) ...
            - (1/2) .* exp(-(speed_ratio .* abs(sin(Incidence))).^2);
            %disp(q_plume)
            if Incidence_wrapped(j) >= deg2rad(active_range(1)) && Incidence_wrapped(j) <= deg2rad(active_range(2))
                Q_input(i) = q_plume * A;
            else
                Q_input(i) = 0;
            end

            % Conductive heat exchange with neighbors
            if m > 1  % Left neighbor
                idx_left = find([nodes.m] == m-1 & [nodes.n] == n);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dx * (SurfTemp(idx_left) - T_i);
            end
            if m < number_of_nodes_width  % Right neighbor
                idx_right = find([nodes.m] == m+1 & [nodes.n] == n);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dx * (SurfTemp(idx_right) - T_i);
            end
            if n > 1  % Bottom neighbor
                idx_bottom = find([nodes.m] == m & [nodes.n] == n-1);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dy * (SurfTemp(idx_bottom) - T_i);
            end
            if n < number_of_nodes_height  % Top neighbor
                idx_top = find([nodes.m] == m & [nodes.n] == n+1);
                Q_conduction(i) = Q_conduction(i) + conductivity * A / dy * (SurfTemp(idx_top) - T_i);
            end

            % Radiative heat loss
            Q_radiation(i) = epsilon * sigma *2* A * (T_i^4 - 300^4);
        end

        % Compute temperature update
        dT = dt * (Q_conduction - Q_radiation + Q_input) / Ci;
        SurfTemp = SurfTemp + dT; % Update node temperatures

        % Store the average temperature
        peak_temp_history(j) = max(SurfTemp);
        %disp(peak_temp_history(j));

        % Update the plot with new data
        %set(h(k), 'XData', 1:j, 'YData', (avg_temp_history(1:j)));

    end

%% 

figTempScenerio = tiledlayout(1,1,'Padding','tight');
figTempScenerio.Units = 'inches';
figTempScenerio.OuterPosition = [0 0 5.5 2.5];
nexttile;
%splot = stackedplot(0+dt:0.1:time_steps/10, [peak_temp_history(:), rad2deg(Incidence_wrapped(:))]);
plot(0+dt:0.1:time_steps/10, peak_temp_history(:));
grid on
% Label x-axis
xlabel('Time [s]');
ylabel('Peak Temperature [K]')
xlim([0 time_steps*0.1])
% Set y-axis labels for each row
%splot.DisplayLabels = {'Peak Temperature [K]', 'Incidence Angle [\circ]'};

exportgraphics(figTempScenerio, "Report/Figures/TempWithScenario.pdf",'ContentType','vector')

figure()
plot(0+dt:1:time_steps, Rate)

%% Pressure and Shear Stress

sigma_n = 0.9;
sigma_t = 0.9;
surf_temp = 300;
plume.Pressure = @(r,phi,theta) 0.5* plume.rho(r,phi) * plume.uRarefied(r,phi)^2 * (1/(plume.SpeedRatio(r,phi))^2) * (   ( (2-sigma_n)/sqrt(pi) * plume.SpeedRatio(r,phi) * sin(theta) + sigma_n/2 * sqrt(surf_temp/plume.TRarefied(r,phi))*exp(-(plume.SpeedRatio(r,phi) * sin(theta))^2)) + (   (2-sigma_n)*((plume.SpeedRatio(r,phi) * sin(theta))^2 + 0.5) + sigma_n/2*sqrt(pi * surf_temp/plume.TRarefied(r,phi))* plume.SpeedRatio(r,phi) * sin(theta)) * (1 + erf(plume.SpeedRatio(r,phi) * sin(theta))));
plume.ShearStress = @(r,phi,theta) 0.5* plume.rho(r,phi) * plume.uRarefied(r,phi)^2 * sigma_t *cos(theta) * 1/sqrt(pi) * ( exp(-(plume.SpeedRatio(r,phi)*sin(theta))^2)  + sqrt(pi)*(plume.SpeedRatio(r,phi)*sin(theta))*(1 + erf(plume.SpeedRatio(r,phi) * sin(theta))));

AoA = linspace(pi/1000,pi,1000);
pressure=zeros(1,length(AoA));
shearstress=zeros(1,length(AoA));

for i= 1:length(AoA)
    pressure(i) = plume.Pressure(10,0,AoA(i));
    shearstress(i) = plume.ShearStress(10,0,AoA(i));
end

figure()
plot(rad2deg(AoA), pressure, rad2deg(AoA), shearstress, rad2deg(AoA), pressure+shearstress)
xlabel("Angle of Attack \theta [\circ]")
ylabel('Stress \sigma [Pa]')
legend("Pressure", "Shear Stress", "Pressure + Shear Stress")
grid on
exportgraphics(gcf, "Report/Figures/StressWithAngle.pdf",'ContentType','vector')
%% TEMP SECTION

centreline.T = zeros(1,length(r));
for i = 1:length(r)
    centreline.T(i) = plume.TRarefied(r(i),0);
end
centreline.u = zeros(1,length(r));
for i = 1:length(r)
    centreline.u(i) = plume.uRarefied(r(i),0);
end
figure();
splot = stackedplot(r, [ centreline.u(:),centreline.T(:)]);

% Label axes
xlabel('r along centreline [m]');

% Remove DisplayLabels to avoid horizontal text
splot.DisplayLabels = {'', ''};

% Add proper y-axis labels manually
ax = findobj(splot.NodeChildren, 'Type', 'Axes');
ylabel(ax(1), 'T [K]', 'Rotation', 90, 'VerticalAlignment', 'bottom','FontSize',12);
ylabel(ax(2), 'u [m/s]', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'FontSize',12);

% Format x-axis to show 10^-3, 10^-2, etc.
set(ax, 'XScale', 'log','FontSize', 12);
% Format x-axis to show 10^-3, 10^-2, etc.
set(ax(1), 'YScale', 'lin','FontSize', 12);
% Format x-axis to show 10^-3, 10^-2, etc.
set(ax(2), 'YScale', 'lin','FontSize', 12);

% Manually set the y-limits for each axis
ylim(ax(2), [2100, 2400]);  % For 'u [m/s]'
ylim(ax(1), [0, 300]);  % For 'T [K]'
xlim([0,100])

grid on

% Set tick label font size
set(ax(1).XAxis, 'FontSize', 20);
set(ax(1).YAxis, 'FontSize', 20);
set(ax(2).XAxis, 'FontSize', 20);
set(ax(2).YAxis, 'FontSize', 20);

exportgraphics(gcf, "Report/Figures/CentrelineTandU.pdf", 'ContentType','vector')

%% Thermal Stress
thermal_stress = zeros(1,length(GlassProps.Type));
max_temp_change = [9];
for i=1:length(max_temp_change)
    for j=1:length(GlassProps.Type)
        thermal_stress(i,j) = max_temp_change(i).*GlassProps.YoungsMod(j)*GlassProps.ThermExpan(j);
        fprintf('Thermal Stress: %.2f MPa\n', thermal_stress ./ 1e6);
    end
end
%% 


a=linspace(0.00001, 0.01,11000);
Glass_Fracture_Toughness = 0.73e6;
Fracture_Energy = GlassProps.FracTough.^2 ./ GlassProps.YoungsMod;
critical_fracture_stress = zeros(length(GlassProps.Type),length(a));
for i = 1:length(GlassProps.Type)
    critical_fracture_stress(i,:) = sqrt((Fracture_Energy(i)*GlassProps.YoungsMod(i))./(pi.*a));
end

figCritStress = tiledlayout(1,1,'Padding','tight');
figCritStress.Units = 'inches';
figCritStress.OuterPosition = [0 0 5.5 2.5];
nexttile;

place(1) = loglog(a.*1000,critical_fracture_stress(1,:)./10^6)
hold on
place(2) = loglog(a.*1000,critical_fracture_stress(2,:)./10^6)
place(3) = loglog(a.*1000,critical_fracture_stress(3,:)./10^6)
defaultColors = get(gca, 'ColorOrder');
yline(thermal_stress(1)/1e6, 'Color', defaultColors(1, :),'LineStyle', '--')
yline(thermal_stress(2)/1e6, 'Color', defaultColors(2, :),'LineStyle','--')
yline(thermal_stress(3)/1e6, 'Color', defaultColors(3, :),'LineStyle','--')
xlabel('Initial crack length a [mm]')
ylabel('Critical Fracture Stress [MPa]')
ylim([0.2 200])
%legend([yline.1, yline.2, yline.3, yline.4])
grid on
legend([place(1), place(2), place(3)], {GlassProps.Type{1},GlassProps.Type{2},GlassProps.Type{3} })

exportgraphics(figCritStress, "Report/Figures/CritFractStress.pdf", 'ContentType','vector')

%% functions

function [xAdj, yAdj, vert_alignment] = adjustTextPosition(x, y, X, Y, xLimits, yLimits)
    % Check if the text is outside x and y limits
    if x < xLimits(1) || x > xLimits(2) || y < yLimits(1) || y > yLimits(2)
        % Find the last valid point within the limits
        validIdx = (X >= xLimits(1) & X <= xLimits(2) & Y >= yLimits(1) & Y <= yLimits(2));
        
        if any(validIdx)
            xAdj = X(find(validIdx, 1, 'first')+50); % Last valid X point
            yAdj = Y(find(validIdx, 1, 'first')+50); % Last valid Y point
            vert_alignment = 'bottom';
        else
            xAdj = x; % Default to original if no valid points exist
            yAdj = y;
            vert_alignment = 'bottom';
        end
    else
        % If within limits, keep original position
        xAdj = x;
        yAdj = y;
        vert_alignment = 'bottom';
    end
end


%% 

function R_f = getRf(phi, plume)
    freeze_r = @(r) plume.Bird(r,phi) -2;

    r_guess = 5;

    R_f = fzero(freeze_r , r_guess);
end


% Clearing workspace, command window, and closing all figures
clear all;
close all;
clc;
tic;  % Starting the timer

% Constants
epsilon_null = 8.8541878128e-12; % Permittivity of free space

% Geometrical Properties
w = 6e-3; % width in meters
l = 14e-3; % length in meters
h = 1.6e-3; % height in meters
A = w * h; % cross-sectional area in square meters
g = 1e-3; % gap in meters

%Actuation Properties
voltage_AC = 180; % AC Voltage in V
voltage_DC = 360; % DC Voltage in V

% Material Properties
epsilon_m = 1; % Relative permittivity of the material
epsilon_d = 1.006;
m = 4.2e-6; % mass in kilograms
rho = 1.3; % density in kg/m^3
EI = 2.79e-8; % Bending stiffness
k = 3 * EI / l^3; % Spring constant
Q = 10; % Quality factor
zheta = 1 / (2 * Q); % Damping ratio

% Natural Frequency Computations
omega_non = 3.51602;
omega_coeff = sqrt(EI / (rho * A * l^4));
omega_r = omega_non * omega_coeff;

% Spatial Coordinates
x_vector = linspace(0, l, 1000);
x_vector_non = x_vector / l;

% Time Coordination
time = [0 1];
time_non = time / sqrt((rho * A * l^4) / EI);

% Mode shape function (normalized)
phi = @(x) cosh(sqrt(omega_non) * x) - cos(sqrt(omega_non) * x) ...
    - 0.7341 * (sinh(sqrt(omega_non) * x) - sin(sqrt(omega_non) * x));

phi_vector = phi(x_vector_non);
phi_max = phi(1);
phi_force = integral(phi, 0, 1);

% Initial Conditions
U_initial = @(x) 0;
U_dot_initial = @(x) 0;
U_initial_non = integral(@(x) U_initial(x) .* phi(x), 0, 1);
U_dot_initial_non = integral(@(x) U_dot_initial(x) .* phi(x), 0, 1);
IC = [U_initial_non; U_dot_initial_non];

% Simulation Setup for Response vs Frequency
jj = 1;
Amp_vector = zeros(length(700:1:1000), 2);

% Sweeping over a range of frequencies
for omega = 700:1:1000
    syms y(t)
    
    % Formulating the differential equation
    [V] = odeToVectorField(diff(y, 2) + 2 * zheta * omega_non * diff(y) ...
        + (omega_non^2) * y == ((l^3 / EI) * phi_force * (1 / (2 * pi)) ...
        * epsilon_null * (epsilon_d - epsilon_m) * h / ((g + y) * ((g + y) + h)) ...
        * (voltage_DC + voltage_AC * sin((omega / omega_coeff) * t))^2));
    
    M = matlabFunction(V, 'vars', {'t', 'Y'});
    opts = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-14);
    sol = ode45(M, time_non, IC, opts);
    
    % Amplitude computation
    Amp_vector(jj, :) = [0.5 * l * (max(sol.y(1, end/2:end)) ...
        - min(sol.y(1, end/2:end))) * phi_max, omega];
    jj = jj + 1;
    
    % Plotting for specific frequency close to resonance
    if omega == round(omega_r)
        plotResponseAtResonance(l, sol, phi_max, x_vector, phi_vector, rho, A, EI);
    end
end

% Plotting Frequency Response
figure;
plot(1 / (2 * pi) * Amp_vector(:, 2), Amp_vector(:, 1) * 1e6, 'color', 'r', 'LineWidth', 2.5);
ylabel('Amplitude (\mu m)');
xlabel('Frequency (Hz)');
saveas(gcf, 'Frq.png');

toc % Ending the timer

function plotResponseAtResonance(l, sol, phi_max, x_vector, phi_vector, rho, A, EI)
    % This function plots the Time response of the cantilever at a specific frequency 
    % close to resonance and creates a video of the displacement over time.

    time_r = (sqrt((rho * A * l^4) / EI)) * sol.x;
    amp_r = l * sol.y(1, :) * phi_max;

    figure;
    plot(time_r, amp_r * 1e6);
    ylabel('Displacement (\mu m)');
    xlabel('Time (sec)');

    % Video creation
    videoFileName = 'my_video.mp4';
    writerObj = VideoWriter(videoFileName, 'MPEG-4');
    writerObj.FrameRate = 30;
    open(writerObj);
    max_amp_time = 0;

    for ii = 1:10:2800
        
        plot(x_vector * 1e3, l * sol.y(1, ii) * phi_vector * 1e6, 'color', 'k', 'LineWidth', 8);
        hold on;
        title(['Time = ', num2str(round(time_r(ii) * 1e3)), ' ms']);

        if max_amp_time < max(l * sol.y(1, ii) * phi_vector * 1e6)
            max_amp_time = max(l * sol.y(1, ii) * phi_vector * 1e6);
        end
        
        plot([0, 20], [max_amp_time, max_amp_time], '--', 'Color', 'b');
        fill([0, 0, 25, 25, 1, 1], [0, -3, -3, -1.5, -1.5, 0], 'r');
        ylim([-2, 2]);
        ylabel('Y (\mu m)');
        xlim([0, 20]);
        xlabel('x (mm)');
        hold off;

        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end
    
    close(writerObj);
end

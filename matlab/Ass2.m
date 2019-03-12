%% Question 1: Environment/Acceleration Models

path = "C:\tudatBundle\tudatApplications\PropagationOptimisation\SimulationOutput\ShapeOptimization\";

benchmark = dlmread(strcat(path, 'stateHistory_0.dat'));
sph_harm = dlmread(strcat(path, 'stateHistory_interpolated_1.dat'));
nrlmsis = dlmread(strcat(path, 'stateHistory_interpolated_2.dat'));
moon = dlmread(strcat(path, 'stateHistory_interpolated_3.dat'));
sun = dlmread(strcat(path, 'stateHistory_interpolated_4.dat'));
rad_pres = dlmread(strcat(path, 'stateHistory_interpolated_5.dat'));
oblateness = dlmread(strcat(path, 'stateHistory_interpolated_6.dat'));
jupiter = dlmread(strcat(path, 'stateHistory_interpolated_7.dat'));


error_sph_harm = sph_harm - benchmark;
error_nrlmsis = nrlmsis - benchmark;
error_moon = moon - benchmark;
error_sun = sun - benchmark;
error_jupiter = jupiter - benchmark;
error_rp = rad_pres - benchmark;
error_oblateness = oblateness - benchmark;

error_sph_harm = sqrt(error_sph_harm(:, 2).^2 + error_sph_harm(:, 3).^2 + error_sph_harm(:, 4).^2);
error_nrlmsis = sqrt(error_nrlmsis(:, 2).^2 + error_nrlmsis(:, 3).^2 + error_nrlmsis(:, 4).^2);
error_moon = sqrt(error_moon(:, 2).^2 + error_moon(:, 3).^2 + error_moon(:, 4).^2);
error_sun = sqrt(error_sun(:, 2).^2 + error_sun(:, 3).^2 + error_sun(:, 4).^2);
error_jupiter = sqrt(error_jupiter(:, 2).^2 + error_jupiter(:, 3).^2 + error_jupiter(:, 4).^2);
error_rp = sqrt(error_rp(:, 2).^2 + error_rp(:, 3).^2 + error_rp(:, 4).^2);
error_oblateness = sqrt(error_oblateness(:, 2).^2 + error_oblateness(:, 3).^2 + error_oblateness(:, 4).^2);

time = benchmark(:, 1);
time = time - time(1);

semilogy(time(10:end-10), error_sph_harm(10:end-10), time(10:end-10), error_nrlmsis(10:end-10), ...
    time(10:end-10), error_moon(10:end-10), time(10:end-10), error_sun(10:end-10), time(10:end-10), ...
    error_jupiter(10:end-10), time(10:end-10), error_rp(10:end-10), ...
    time(10:end-10), error_oblateness(10:end-10), 'Linewidth', 1.5);
xlabel('Time [s]');
ylabel('Position error [m]');
grid on;
legend('Spherical harmonics up to 2/2', 'NRLMSISE00 atmospheric model', 'Incl. Lunar gravity', ...
    'Incl. Solar gravity', 'Incl. Jovian gravity', 'Solar radiation pressure', 'Oblate Earth');
title('Question 1');

%% Question 1

% High-fidelity orbit propagation

path = "C:\tudatBundle\tudatApplications\PropagationOptimisation\SimulationOutput\ShapeOptimization\";

state_0_0 = dlmread(strcat(path, 'stateHistory_0_0.dat'));
state_1_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_1_0.dat'));
state_2_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_2_0.dat'));
state_3_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_3_0.dat'));
state_4_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_4_0.dat'));
state_5_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_5_0.dat'));
state_6_0_interpolated = dlmread(strcat(path, 'stateHistory_interpolated_6_0.dat'));

depvar_0_0 = dlmread(strcat(path, 'dependentVariables_0_0.dat'));
depvar_1_0 = dlmread(strcat(path, 'dependentVariables_interpolated_1_0.dat'));
depvar_2_0 = dlmread(strcat(path, 'dependentVariables_interpolated_2_0.dat'));
depvar_3_0 = dlmread(strcat(path, 'dependentVariables_interpolated_3_0.dat'));
depvar_4_0 = dlmread(strcat(path, 'dependentVariables_interpolated_4_0.dat'));
depvar_5_0 = dlmread(strcat(path, 'dependentVariables_interpolated_5_0.dat'));
depvar_6_0 = dlmread(strcat(path, 'dependentVariables_interpolated_6_0.dat'));

error_1 = state_1_0_interpolated - state_0_0;
error_2 = state_2_0_interpolated - state_0_0;
error_3 = state_3_0_interpolated - state_0_0;
error_4 = state_4_0_interpolated - state_0_0;
error_5 = state_5_0_interpolated - state_0_0;
error_6 = state_6_0_interpolated - state_0_0;

alt_error_1 = abs(depvar_1_0(:,2) - depvar_0_0(:,2));
alt_error_2 = abs(depvar_2_0(:,2) - depvar_0_0(:,2));
alt_error_3 = abs(depvar_3_0(:,2) - depvar_0_0(:,2));
alt_error_4 = abs(depvar_4_0(:,2) - depvar_0_0(:,2));
alt_error_5 = abs(depvar_5_0(:,2) - depvar_0_0(:,2));
alt_error_6 = abs(depvar_6_0(:,2) - depvar_0_0(:,2));

error_1_pos = sqrt(error_1(:, 2).^2 + error_1(:, 3).^2 + error_1(:, 4).^2);
error_2_pos = sqrt(error_2(:, 2).^2 + error_2(:, 3).^2 + error_2(:, 4).^2);
error_3_pos = sqrt(error_3(:, 2).^2 + error_3(:, 3).^2 + error_3(:, 4).^2);
error_4_pos = sqrt(error_4(:, 2).^2 + error_4(:, 3).^2 + error_4(:, 4).^2);
error_5_pos = sqrt(error_5(:, 2).^2 + error_5(:, 3).^2 + error_5(:, 4).^2);
error_6_pos = sqrt(error_6(:, 2).^2 + error_6(:, 3).^2 + error_6(:, 4).^2);

time = state_0_0(:, 1);
time = time - time(1);

% figure;
plot(time(10:end-10), depvar_0_0(10:end-10, 2), time(10:end-10), depvar_1_0(10:end-10, 2), time(10:end-10), depvar_2_0(10:end-10, 2),...
    time(10:end-10), depvar_3_0(10:end-10, 2), time(10:end-10), depvar_5_0(10:end-10, 2), time(10:end-10), depvar_5_0(10:end-10, 2));
grid;
xlabel('Time [s]');
ylabel('Altitude [m]');
legend('RK78', 'RK45', 'RK56', 'RK78', 'DP87', 'ABM', 'BS');

figure;
%subplot(1, 2, 1);
semilogy(time(10:end-10), error_1_pos(10:end-10), time(10:end-10), error_2_pos(10:end-10), time(10:end-10),...
    error_3_pos(10:end-10), time(10:end-10), error_4_pos(10:end-10), time(10:end-10), error_5_pos(10:end-10),...
    time(10:end-10), error_6_pos(10:end-10), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Position error [m]');
grid on;
legend('RK4(5)', 'RK5(6)', 'RK7(8)', 'DP87', 'ABM', 'BS');
title('Position integration error for different integrators');

figure;
%subplot(1, 2, 2);
semilogy(time(10:end-10), alt_error_1(10:end-10), time(10:end-10), alt_error_2(10:end-10),  time(10:end-10), alt_error_3(10:end-10),...
     time(10:end-10), alt_error_4(10:end-10), time(10:end-10), alt_error_5(10:end-10),...
     time(10:end-10), alt_error_6(10:end-10), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Altitude error [m]');
grid on;
legend('RK4(5)', 'RK5(6)', 'RK7(8)', 'DP87', 'ABM', 'BS');
title('Altitude error');



%% Question 2

state_0_0 = dlmread(strcat(path, 'stateHistory_0_0.dat'));
state_7_0 = dlmread(strcat(path, 'stateHistory_interpolated_7_0.dat'));
state_7_1 = dlmread(strcat(path, 'stateHistory_interpolated_7_1.dat'));
state_7_2 = dlmread(strcat(path, 'stateHistory_interpolated_7_2.dat'));
state_7_3 = dlmread(strcat(path, 'stateHistory_interpolated_7_3.dat'));
state_7_4 = dlmread(strcat(path, 'stateHistory_interpolated_7_4.dat'));
state_7_5 = dlmread(strcat(path, 'stateHistory_interpolated_7_5.dat'));
state_7_6 = dlmread(strcat(path, 'stateHistory_interpolated_7_6.dat'));

time = state_0_0(:, 1);
time = time - time(1);
% 
% figure;
% plot3(state_0_0(:, 2), state_0_0(:, 3), state_0_0(:, 4));
% xlabel('x [m]');
% ylabel('y [m]');
% xlabel('z [m]');
% grid;

error_0 = state_7_0 - state_0_0;
error_1 = state_7_1 - state_0_0;
error_2 = state_7_2 - state_0_0;
error_3 = state_7_3 - state_0_0;
error_4 = state_7_4 - state_0_0;
error_5 = state_7_5 - state_0_0;
error_6 = state_7_6 - state_0_0;

error_0_pos = sqrt(error_0(:, 2).^2 + error_0(:, 3).^2 + error_0(:, 4).^2);
error_1_pos = sqrt(error_1(:, 2).^2 + error_1(:, 3).^2 + error_1(:, 4).^2);
error_2_pos = sqrt(error_2(:, 2).^2 + error_2(:, 3).^2 + error_2(:, 4).^2);
error_3_pos = sqrt(error_3(:, 2).^2 + error_3(:, 3).^2 + error_3(:, 4).^2);
error_4_pos = sqrt(error_4(:, 2).^2 + error_4(:, 3).^2 + error_4(:, 4).^2);
error_5_pos = sqrt(error_5(:, 2).^2 + error_5(:, 3).^2 + error_5(:, 4).^2);
error_6_pos = sqrt(error_6(:, 2).^2 + error_6(:, 3).^2 + error_6(:, 4).^2);

figure;
semilogy(time(10:end-10), error_0_pos(10:end-10), time(10:end-10), error_1_pos(10:end-10), time(10:end-10),...
    error_2_pos(10:end-10), '-.', time(10:end-10), error_3_pos(10:end-10), '-.', time(10:end-10), error_4_pos(10:end-10), '--',...
    time(10:end-10), error_5_pos(10:end-10), '--', time(10:end-10), error_6_pos(10:end-10), '--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Position error [m]');
grid on;
legend('Cowell', 'Encke', 'Kepler', 'MEE', 'USM4', 'USM MRP', 'USM EM');
title('Position integration error for different propagators');

%% Question 3

benchmark = dlmread(strcat(path, 'stateHistory_0_0.dat'));
time = benchmark(:, 1);
time = time - time(1);

% Cowell
cow45 = dlmread(strcat(path, 'stateHistory_interpolated_8_0.dat'));
cow56 = dlmread(strcat(path, 'stateHistory_interpolated_9_0.dat'));
cow78 = dlmread(strcat(path, 'stateHistory_interpolated_10_0.dat'));

error_cow45 = cow45 - benchmark;
error_cow56 = cow56 - benchmark;
error_cow78 = cow78 - benchmark;

cow45_pos = sqrt(error_cow45(:, 2).^2 + error_cow45(:, 3).^2 + error_cow45(:, 4).^2);
cow56_pos = sqrt(error_cow56(:, 2).^2 + error_cow56(:, 3).^2 + error_cow56(:, 4).^2);
cow78_pos = sqrt(error_cow78(:, 2).^2 + error_cow78(:, 3).^2 + error_cow78(:, 4).^2);

% Encke
encke45 = dlmread(strcat(path, 'stateHistory_interpolated_8_1.dat'));
encke56 = dlmread(strcat(path, 'stateHistory_interpolated_9_1.dat'));
encke78 = dlmread(strcat(path, 'stateHistory_interpolated_10_1.dat'));

error_encke45 = encke45 - benchmark;
error_encke56 = encke56 - benchmark;
error_encke78 = encke78 - benchmark;

encke45_pos = sqrt(error_encke45(:, 2).^2 + error_encke45(:, 3).^2 + error_encke45(:, 4).^2);
encke56_pos = sqrt(error_encke56(:, 2).^2 + error_encke56(:, 3).^2 + error_encke56(:, 4).^2);
encke78_pos = sqrt(error_encke78(:, 2).^2 + error_encke78(:, 3).^2 + error_encke78(:, 4).^2);

% Plot data
figure;
semilogy(time(10:end-10), cow45_pos(10:end-10), time(10:end-10), cow56_pos(10:end-10), time(10:end-10),...
    cow78_pos(10:end-10), time(10:end-10), encke45_pos(10:end-10), '-.', time(10:end-10), encke56_pos(10:end-10), '-.',...
    time(10:end-10), encke78_pos(10:end-10), '-.', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Position error [m]');
grid on;
legend('Cowell RK45', 'Cowell RK56', 'Cowell RK78', 'Encke RK45', 'Encke RK56', 'Encke RK78');
title('Position integration error Cowell/Encke with tolerance of 10^{-5}, min. 0.1 s, max 0.1 s');
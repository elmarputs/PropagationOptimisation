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

%% Question 3: Analysis of uncertainty in initial state

path = "C:\tudatBundle\tudatApplications\PropagationOptimisation\SimulationOutput\ShapeOptimization\";

reference    = dlmread(strcat(path, 'MC_0_dep.dat'));
uncert_alt   = dlmread(strcat(path, 'MC_1_dep.dat'));
uncert_vel   = dlmread(strcat(path, 'MC_2_dep.dat'));
uncert_gamma = dlmread(strcat(path, 'MC_3_dep.dat'));

latlon_alt_deg = rad2deg(uncert_alt(:,7:8));
latlon_vel_deg = rad2deg(uncert_vel(:,7:8));
figure;
histogram2(latlon_alt_deg(:,1), latlon_alt_deg(:,2), 12);

figure;
hold on;
scatter(latlon_alt_deg(:,2), latlon_alt_deg(:,1), 50, '.');
scatter(rad2deg(reference(8)), rad2deg(reference(7)), 50, 'r.');
legend('Monte Carlo results', 'Nominal run');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
grid on;

figure;
hold on;
scatter(latlon_vel_deg(:,2), latlon_vel_deg(:,1), 50, '.');
scatter(rad2deg(reference(8)), rad2deg(reference(7)), 50, 'r.');
legend('Monte Carlo results', 'Nominal run');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
grid on;

h = heatscatter(latlon_alt_deg(:,2), latlon_alt_deg(:,1), '', 'heatmap.eps', '20');

covar_lat_lon = cov(latlon_alt_deg(:,1), latlon_alt_deg(:,2));
ellipse_axes = sqrt(eig(covar_lat_lon));
sigma_lat = sqrt(covar_lat_lon(1,1));
sigma_lon = sqrt(covar_lat_lon(2,2));

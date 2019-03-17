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
set(gca, 'FontSize', 14);

%% Question 3: Analysis of uncertainty in initial state

path = "C:\tudatBundle\tudatApplications\PropagationOptimisation\SimulationOutput\ShapeOptimization\";
r_e  = 6378+25; %km

reference    = dlmread(strcat(path, 'MC_0_dep.dat'));
uncert_alt   = dlmread(strcat(path, 'MC_1_dep.dat'));
uncert_vel   = dlmread(strcat(path, 'MC_2_dep.dat'));
uncert_gamma = dlmread(strcat(path, 'MC_3_dep.dat'));
uncert_lat = dlmread(strcat(path, 'MC_4_dep.dat'));
uncert_lon = dlmread(strcat(path, 'MC_5_dep.dat'));
uncert_hdg = dlmread(strcat(path, 'MC_6_dep.dat'));
uncert_mass = dlmread(strcat(path, 'MC_7_dep.dat'));
uncert_flat = dlmread(strcat(path, 'MC_8_dep.dat'));

latlon_alt   = uncert_alt(:,7:8);
latlon_vel   = uncert_vel(:,7:8);
latlon_gamma = uncert_gamma(:,7:8);
latlon_lat = uncert_lat(:,7:8);
latlon_lon = uncert_lon(:,7:8);
latlon_hdg = uncert_hdg(:,7:8);
latlon_mass = uncert_mass(:,7:8);
latlon_flat = uncert_flat(:,7:8);

r_alt = vecnorm(latlon_alt, 2, 2).*r_e;
r_vel = vecnorm(latlon_vel, 2, 2).*r_e;
r_gamma = vecnorm(latlon_gamma, 2, 2).*r_e;
r_lat = vecnorm(latlon_lat, 2, 2).*r_e;
r_lon = vecnorm(latlon_lon, 2, 2).*r_e;
r_hdg = vecnorm(latlon_hdg, 2, 2).*r_e;
r_mass = vecnorm(latlon_mass, 2, 2).*r_e;
r_flat = vecnorm(latlon_flat, 2, 2).*r_e;

figure;

subplot(2, 3, 1);
histogram(r_alt, 15);
title('Altitude');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

subplot(2, 3, 2);
histogram(r_vel, 15);
title('Velocity');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

subplot(2, 3, 3);
histogram(r_gamma, 15);
title('Flight path angle');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

subplot(2, 3, 4);
histogram(r_lat, 15);
title('Latitude');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

subplot(2, 3, 5);
histogram(r_lon, 15);
title('Longitude');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

subplot(2, 3, 6);
histogram(r_hdg, 15);
title('Heading');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid

figure;
subplot(1,2,1);
histogram(r_mass, 25);
title('Vehicle density');
xlabel('Distance from (? = 0, ? = 0) [km]')
ylabel('Frequency [#/bin]')
set(gca, 'FontSize', 14);
grid;

subplot(1,2,2);
histogram(r_flat, 25);
title('Earth ellipsoid flatness');
xlabel('Distance from (? = 0, ? = 0) [km]');
ylabel('Frequency [#/bin]');
set(gca, 'FontSize', 14);
grid;

%%% 3D Histogram

figure;
subplot(1, 2, 1);
histogram2(rad2deg(latlon_alt(:,2)), rad2deg(latlon_alt(:,1)), 15);
xlabel('Longitude [\circ]')
ylabel('Latitude [\circ]')
grid on;
set(gca, 'FontSize', 14);

subplot(1, 2, 2);
histogram2(rad2deg(latlon_alt(:,2)), rad2deg(latlon_alt(:,1)), 15);
xlabel('Longitude [\circ]')
ylabel('Latitude [\circ]')
grid on;
set(gca, 'FontSize', 14);


%%% Scatter plot

figure;
hold on;
scatter(rad2deg(latlon_alt(:,2)), rad2deg(latlon_alt(:,1)), 50, '.');
scatter(rad2deg(reference(8)), rad2deg(reference(7)), 80, 'r.');
legend('Monte Carlo results', 'Nominal run');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
grid on;
set(gca, 'FontSize', 14);
% 
% figure;
% hold on;
% scatter(latlon_vel_deg(:,2), latlon_vel_deg(:,1), 50, '.');
% scatter(rad2deg(reference(8)), rad2deg(reference(7)), 50, 'r.');
% legend('Monte Carlo results', 'Nominal run');
% xlabel('Longitude [deg]');
% ylabel('Latitude [deg]');
% grid on;
% 
% h = heatscatter(latlon_alt_deg(:,2), latlon_alt_deg(:,1), '', 'heatmap.eps', '20');
% 
% covar_lat_lon = cov(latlon_alt_deg(:,1), latlon_alt_deg(:,2));
% ellipse_axes = sqrt(eig(covar_lat_lon));
% sigma_lat = sqrt(covar_lat_lon(1,1));
% sigma_lon = sqrt(covar_lat_lon(2,2));

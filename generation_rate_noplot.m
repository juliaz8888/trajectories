clear;
clc;

%% Load data
load('CdSe_trajs.mat');

%% Set up grid 
max_length = 150; % maximum length of grid/grid boundary, in um, make sure this includes all the data points 
N_grid = 51; % number of grid elements in each dimension 

% Set up grids 
x = linspace(0, max_length, N_grid);
y = linspace(0, max_length, N_grid);
z = linspace(0, max_length, N_grid);

% Energy lost in each grid 
E_lost = zeros(N_grid, N_grid, N_grid);

%% Assign energy lost to each grid element 
% Combine trajs

for i = 1:100
     try        
        traj_all = vertcat(traj_all, eval(['traj' num2str(i)]));
        trajs{i} = eval(['traj' num2str(i)]);
     catch
     end
end

size = length(traj_all);

for i = size:-1:1
    if abs(traj_all(i, 1)) > 10
        traj_all(i,:) = [];
    end
end

for itr = 1:length(trajs) % change to however many generation rates, but this will plot a bunch of figures
    traj = cell2mat(trajs(itr));
    % Process traj data so that every dimension starts from zero 
    % Can offset this later 
    % Originally in cm, convert this to um 
    x_traj = traj(:,1) * 1e4;
    y_traj = traj(:,2) * 1e4;
    z_traj = traj(:,3) * 1e4;

    x_traj = x_traj - min(traj(:,1)) * 1e4; % change for each trajectory
    y_traj = y_traj - min(traj(:,2)) * 1e4;
    z_traj = z_traj - min(traj(:,3)) * 1e4;
    
    % I want to offset the x, y, z coordinates to put them at the center of the
    % simulation region. For the simulation region, I will go all dimensions
    % from 0-150 um. 
    x_traj = x_traj + 30;
    y_traj = y_traj + 5;
    z_traj = z_traj + 55;

    % Energy of electron at each location, in eV 
    E_left = traj(:,4);

    % The energy of electron at the previous location, starting at the first
    % value of E_left 
    E_prev = E_left(1);
    
    % Iterate over the electrons and assign the delta_energy to each grid 
    for i = 1:length(x_traj)
        [x_val, x_ind] = min(abs(x - x_traj(i)));
        [y_val, y_ind] = min(abs(y - y_traj(i)));
        [z_val, z_ind] = min(abs(z - z_traj(i)));
        % Energy always decreases over rows 
        delta_E = E_prev - E_left(i);
        E_lost(x_ind, y_ind, z_ind) = E_lost(x_ind, y_ind, z_ind) + delta_E; 
        E_prev = E_left(i);
    end
% end

%% For plot only
x_plot = zeros(N_grid * N_grid * N_grid, 1);
y_plot = zeros(N_grid * N_grid * N_grid, 1);
z_plot = zeros(N_grid * N_grid * N_grid, 1);
E_plot = zeros(N_grid * N_grid * N_grid, 1);
for i = 1:N_grid
    for j = 1:N_grid
        for k = 1:N_grid 
            ind = (i - 1) * N_grid * N_grid + ...
                (j - 1) * N_grid + k;
            x_plot(ind) = x(i);
            y_plot(ind) = y(j);
            z_plot(ind) = z(k);
            E_plot(ind) = E_lost(i, j, k);
        end
    end
end
%% For plotting, show the energy depostion profile 
% eps = 1e-10;
% scale = 4000;
%  E_plot = (E_plot/ max(E_plot) + eps) * scale;
% figure();
% scatter3(x_plot, y_plot, z_plot, E_plot, '.');
% xlabel('x (\mum)');
% ylabel('y (\mum)');
% zlabel('z (\mum)');
% title('(a). CdSe');
% axis([20 160 0 160 45 160]);
% set(gca, 'FontSize',22);
%% Convert energy lost to generation rate 
% Generation rate has a unit of m^-3*s^-1
% For CdSe, the electron-hole pair production energy is 5.22 eV
eh_gen = 5.22; 
% The volume for a grid, convert to m^3 
grid_volume = (x(2) - x(1)) * 1e-6 * ...
    (y(2) - y(1)) * 1e-6 * ...
    (z(2) - z(1)) * 1e-6; 
% The generation time for the electrons, say 1 ps, convert to s 
gen_time = 1e-12; 
G = E_lost / eh_gen / grid_volume / gen_time; 

% x, y, z have units of m
x = x * 1e-6;
y = y * 1e-6;
z = z * 1e-6;

% Finally need to save (x, y, z, G) in a .mat file 
% save('all_charge_track.mat', 'x', 'y', 'z', 'G')
save(['generation_rate' num2str(itr) '.mat'], 'x', 'y', 'z', 'G')
end











clear; clc; close all;

%%% TODO: find a meaningful case study
P = [ 20 5 4 20 4 1
    2 3 1 10 3 3
    10 4 5 3 10 9
    2 5 4 2 4 10]; % Processing time job j on machine m (JxM matrix)
G_init0 = [1 2 3 5 0
      1 2 3 4 1
      2 4 4 0 0
      2 1 5 4 0
      1 2 3 6 0
      1 4 5 2 5]; % graph path (AxN matrix, N = max(length(Jobs))
G_init = G_init0;
G_j0 = [1
      1
      2
      3
      4
      4]; % alternatives related to the jobs (Ax1 vector)
G_j = G_j0;
BigOmega = 1:5:11;
Release_planned = [0 0 5 7]';
S0 = Release_planned;
max_delay = 3; % max advance/delay on a job w.r.t. planned release time
horizon = 2; % prediction horizon for MPC-scheduling
Release_real = Release_planned + randi([-max_delay max_delay], [length(Release_planned) 1]);
%%% Run the offline and dynamic scheduling for comparison
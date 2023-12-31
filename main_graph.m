%%% MPC-BASED SCHEDULING
clear; clc; close all;

% Parameters of the problem
%%% Case study
P = [9 5 7 10 4 12
     4 7 3 7 1 10
     5 7 6 3 10 1
     4 3 10 6 4 5
     2 4 7 3 5 2
     1 6 5 3 6 8];

G_init0 = [1 2 3 4 5
           1 3 5 6 0
           1 2 3 4 6
           1 2 4 5 0
           1 2 5 6 0
           1 2 3 5 0
           1 3 4 6 0
           1 2 6 5 3
           1 4 6 5 3
           1 2 4 5 6
           1 3 4 5 0
           1 2 3 4 6
           1 3 4 6 0];
G_j0 = [1 1 2 2 3 3 3 4 5 5 6 6 6]';
Release_planned = [0 2 4 7 10 12]';
% Set planned release time and real release time
max_delay = 3; % max advance/delay on a job w.r.t. planned release time
horizon = 2; % prediction horizon for MPC-scheduling
Release_real = Release_planned + randi([-max_delay max_delay], [length(Release_planned) 1]);
Release_real = [0 1 2 4 7 14]';
Release_real(Release_real < 0) = 0; % Set release time to 0 as minimum release time
Release_real = sort(Release_real); % Avoid possible job swap


%%% TODO: Sort data by release time considering possible job swap
% R=[];
% for i=1:length(G_j0)
%     R=[R; Release_real(G_j0(i))];
% end
% sorted = table(G_init0,G_j0,R);
% sorted = sortrows(sorted,3);
% G_init0 = table2array(sorted(:,1));
% G_j0 = table2array(sorted(:,2));
% Release_real = sort(Release_real);
% P = P(unique(G_j0, 'stable'),:);
%%%

% Pre processing of data
M0 = max(max(G_init0));
[G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init0, P);
J = length(unique(G_j0)); %jobs
M = max(max(G)); %machines
A = size(G_j0,1);%alternatives
D = compute_D_from_graph(G_init0,G_j0); % disjunctive connections (2 constraints for each connection)


%% SOLVE PROBLEM
sol_noNoise = [];
events = unique(Release_real); % Find the events (i.e. release of products)
% Loop through each "event" (i.e. arrival of a new job)
paths_atm = zeros(length(events),J);
for t=1:length(events)
    % Consider the jobs currently present in the shopfloor
    idx_job_in_shop = find(Release_real <= events(t)); % Index of jobs in the shop
    jobs_in_shop = Release_real(idx_job_in_shop)'; % Jobs in the shop
    idx_job_predicted = find(Release_planned >= events(t) & ...
            Release_planned <= events(t) + horizon); % Index of job predicted to be in the shop
    job_diff = setdiff(idx_job_predicted,idx_job_in_shop);
    job_predicted = Release_planned(job_diff)'; % Do not consider twice the jobs already in the shop
    S0 = [jobs_in_shop job_predicted];
    % Update open-shop graph structures to consider only jobs in the shop
    % or in the prediction horizon
    G_j = G_j0(G_j0<=length(S0)); 
    G_init = G_init0(G_j0 <= length(S0),:);

    %BigOmega =1:3:10; % Test with different values of noises
    tStart = tic;
    BigOmega = 4;
    % Bouncing algorithm --> Find the best trade-off solution
    for i=1:length(BigOmega)
        u=1;
        vec_gamma = [];
        sol_noNoise = Graph_minimization(G_init,G_j,P, S0, sol_noNoise,M0, events(t));
       solMax(u) = Graph_maximization(G_init,G_j,P,sol_noNoise,BigOmega(i), S0);
        for j=1:20
            u=u+1;
            new_omega = zeros(size(P));
            new_omega(1:size(solMax(u-1).omega,1), ...
                1:size(solMax(u-1).omega,2)) = solMax(u-1).omega; % Fill the new omega with solMax.omega, and the remaining are zeros
           solMin(j) = Graph_minimization(G_init,G_j,P+P.*new_omega, S0, [], M0, events(t));
             solaux = Graph_minimization(G_init,G_j,P+P.*new_omega, S0, [], M0, events(t));
            if(j > 1 & ismember(int8(solaux.gamma)',vec_gamma, 'rows') )
                % if the solution found is already in the pool of
                % solutions, exit the loop after finding the best one in
                % the current pool
                [minVal, minIdx ] = min(vertcat(solMin.C));
                solOpt(i)=solMin(minIdx);
                % Update the disturbances structure to match with the
                % subpart of the scheduling considered (jobs currently in
                % the shop or in the prediction horizon) and the original
                % structure of the shop (all jobs considered)
                new_disturbances = zeros(size(P));
                new_disturbances(1:size(solMax(minIdx).omega,1), ...
                    1:size(solMax(minIdx).omega,2)) = solMax(minIdx).omega;
                disturbances(i,:,:) = new_disturbances;
                break
            else
                solMin(j) = solaux;
                vec_gamma = [vec_gamma; int8(solaux.gamma)'];
            end
             if(j>1 & sum(int8(solMin(j).gamma) == int8(sol_noNoise.gamma)) == length(solMin(j).gamma))
                % If the solution is equivalent to the solution without
                % disturbances, exit the loop (this is the best solution)
                 solOpt(i)=solMin(j);
                new_disturbances = zeros(size(P));
                new_disturbances(1:size(solMax(minIdx).omega,1), ...
                    1:size(solMax(minIdx).omega,2)) = solMax(u-1).omega;
                disturbances(i,:,:) = new_disturbances;
                break
            end
            if (j>=10)
                % After a predefined number of iterations, exit the loop
                [minVal, minIdx ] = min(vertcat(solMin.C));
                solOpt(i)=solMin(minIdx);
                disturbances(i,:,:) = solMax(u-1).omega;
                flag =1;
                break
            end
            solMax(u) = Graph_maximization(G_init,G_j,P,solMin(j),BigOmega(i), S0);
        end
        clear solMax
        clear solMin
    end
    paths_atm(1:size(sol_noNoise.gamma),t) = sol_noNoise.gamma;
end
tEnd = toc(tStart);
    %% Robust analysis
    for i=1:length(BigOmega)
        for j=1:1000
        [robustCompletion(i,j), ~] = graph_minimization_robust(G_init,G_j,P,S0, solOpt(i).gamma, solOpt(i).delta, BigOmega(i));
        end
    end
  %% Plot robust analysis
boxplot(robustCompletion')
hold on
for i=1:length(BigOmega)
    plot(i,solOpt(i).C, "b*")
    %plot(i,min(robustCompletion(i,:)), "ko")
end
title("Robust analysis")
ylabel("Completion time")
xlabel("\Omega")
xticklabels({'1', '4','7','10'})
%Plot test
%% Get noise-free solution
figure()
graph_Gantt(sol_noNoise, G_init, G_j, P, sol_noNoise.gamma, M0, "Noise-free solution");
xlabel("Time")
% Get best solution between noise-free and robust
i = randi(length(solOpt)); % Print random noise solution
i=4;
figure()
graph_Gantt(solOpt(i), G_init, G_j, P, solOpt(i).gamma, M0, "Best solution (trade-off optimal/robust), \Omega= " + BigOmega(i));
xlabel("Time")
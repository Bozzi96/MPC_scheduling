clear; clc; %close all;

% Parameters of the problem
%%% EXAMPLE 1: Manufacturing
G_init0 = [1 2 2 5 1
     1 2 3 5 1
     1 4 2 5 1
     1 4 3 5 1
     1 2 3 5 1
     1 2 4 5 1
     1 4 3 5 1
     1 4 4 5 1
     1 2 2 5 1
     1 2 4 5 1
     1 3 2 5 1
     1 3 4 5 1
    ];
G_j0 = [1
    1
    1
    1
    2
    2
    2
    2
    3
    3
    3
    3
    ];
P = [ 10 20 30 20 10
    10 20 30 20 10
    10 20 30 20 10
    ];

Release_planned = [0 0 30]';

%%% EXAMPLE 2: Quick
P = [ 20 5 4 20 4 1
    2 3 1 10 3 3
    10 4 5 3 10 9
    2 5 4 2 4 10]; % Processing time job j on machine m (JxM matrix)
G_init0 = [1 2 3 5 0
      1 4 1 3 1
      2 4 4 0 0
      2 1 5 4 0
      1 2 3 6 0
      1 4 5 2 5]; % graph path (AxN matrix, N = max(length(Jobs))
G_j0 = [1
      1
      2
      3
      4
      4]; % alternatives related to the jobs (Ax1 vector)
Release_planned = [0 0 5 7]';
max_delay = 3;
horizon = 1;
Release_real = Release_planned + randi([-max_delay max_delay], [length(Release_planned) 1]);
Release_real(Release_real < 0) = 0;
    % Pre processing of data
    M0 = max(max(G_init0));
    [G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init0, P);
    J = length(unique(G_j0)); %jobs
    M = max(max(G)); %machines
    A = size(G_j0,1);%alternatives
    D = compute_D_from_graph(G_init0,G_j0); % disjunctive connections (2 constraints per each connection)
    
%% SOLVE PROBLEM
sol_noNoise = [];
for t=1:length(unique(Release_real))
    idx = unique(Release_real); % Find the events (i.e. release of products)
    % Consider the jobs currently present in the shopfloor
    S0 = Release_planned(Release_planned <= idx(t) + horizon);
    G_j = G_j0(G_j0<=length(S0)); 
    G_init = G_init0(G_j0 <= length(S0),:);

    %graph_plots(solMin, G_init, G_j, P, solMin.gamma);
    %BigOmega =1:2:11;
    BigOmega = 5;
    for i=1:length(BigOmega)
        u=1;
        vec_gamma = [];
        sol_noNoise = Graph_minimization(G_init,G_j,P, S0, sol_noNoise,M0);
%        solMax(u) = Graph_maximization(G_init,G_j,P,sol_noNoise,BigOmega(i), S0);
%         for j=1:20
%             u=u+1;
%            solMin(j) = Graph_minimization(G_init,G_j,P+P.*solMax(u-1).omega, S0);
%              solaux = Graph_minimization(G_init,G_j,P+P.*solMax(u-1).omega, S0);
%             if(j > 1 & ismember(int8(solaux.gamma)',vec_gamma, 'rows') )
%                 % if the solution found is already in the pool of
%                 % solutions, exit the loop after finding the best one in
%                 % the current pool
%                 [minVal, minIdx ] = min(vertcat(solMin.C));
%                 solOpt(i)=solMin(minIdx);
%                 disturbances(i,:,:) = solMax(minIdx).omega;
%                 break
%             else
%                 solMin(j) = solaux;
%                 vec_gamma = [vec_gamma; int8(solaux.gamma)'];
%             end
%              if(j>1 & sum(int8(solMin(j).gamma) == int8(sol_noNoise.gamma)) == length(solMin(j).gamma))
%                 % If the solution is equivalent to the solution without
%                 % disturbances, exit the loop (this is the best solution)
%                  solOpt(i)=solMin(j);
%                 disturbances(i,:,:) = solMax(u-1).omega;
%                 break
%             end
%             if (j>=20)
%                 % After a predefined number of iterations, exit the loop
%                 [minVal, minIdx ] = min(vertcat(solMin.C));
%                 solOpt(i)=solMin(minIdx);
%                 disturbances(i,:,:) = solMax(u-1).omega;
%                 break
%             end
%             solMax(u) = Graph_maximization(G_init,G_j,P,solMin(j),BigOmega(i), S0);
%         end
        clear solMax
        clear solMin
        %graph_plots(solMax(i), G_init, G_j, P, solMin.gamma);
    end
end
    %% Robust analysis
%     for i=1:length(BigOmega)
%         for j=1:1000
%         [robustCompletion(i,j), ~] = graph_minimization_robust(G_init,G_j,P,S0, solOpt(i).gamma, solOpt(i).delta, BigOmega(i));
%         end
%     end
%     boxplot(robustCompletion')

clear; clc; close all;

% Parameters of the problem
G_init = [1 2 2 5 1
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
G_j = [1
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

S0 = [0 0 0]';

% P = [ 2 5 4 2 4 1
%     2 3 1 1 3 3
%     1 4 5 3 1 9
%     2 5 4 2 4 1]; % Processing time job j on machine m (JxM matrix)
% G_init = [1 2 3 5 0
%       1 4 1 3 1
%       2 4 4 0 0
%       2 1 5 4 0
%       1 2 3 6 0
%       1 4 5 2 5]; % graph path (AxN matrix, N = max(length(Jobs))
% G_j = [1
%       1
%       2
%       3
%       4
%       4]; % alternatives related to the jobs (Ax1 vector)
% S0 = [0 0 0 0]';
    % Pre processing of data
    [G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init, P);
    J = length(unique(G_j)); %jobs
    M = max(max(G)); %machines
    A = size(G_j,1);%alternatives
    D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints per each connection)
    

%% SOLVE PROBLEM

    %graph_plots(solMin, G_init, G_j, P, solMin.gamma);
    BigOmega =1:2:11;
    for i=1:length(BigOmega)
        u=1;
        vec_gamma = [];
        sol_noNoise = Graph_minimization(G_init,G_j,P, S0);
        solMax(u) = Graph_maximization(G_init,G_j,P,sol_noNoise,BigOmega(i), S0);
        for j=1:20
            u=u+1;
           solMin(j) = Graph_minimization(G_init,G_j,P+P.*solMax(u-1).omega, S0);
             solaux = Graph_minimization(G_init,G_j,P+P.*solMax(u-1).omega, S0);
            if(j > 1 & ismember(int8(solaux.gamma)',vec_gamma, 'rows') )
                % if the solution found is already in the pool of
                % solutions, exit the loop after finding the best one in
                % the current pool
                [minVal, minIdx ] = min(vertcat(solMin.C));
                solOpt(i)=solMin(minIdx);
                disturbances(i,:,:) = solMax(minIdx).omega;
                break
            else
                solMin(j) = solaux;
                vec_gamma = [vec_gamma; int8(solaux.gamma)'];
            end
             if(j>1 & sum(int8(solMin(j).gamma) == int8(sol_noNoise.gamma)) == length(solMin(j).gamma))
                % If the solution is equivalent to the solution without
                % disturbances, exit the loop (this is the best solution)
                 solOpt(i)=solMin(j);
                disturbances(i,:,:) = solMax(u-1).omega;
                break
            end
            if (j>=20)
                % After a predefined number of iterations, exit the loop
                [minVal, minIdx ] = min(vertcat(solMin.C));
                solOpt(i)=solMin(minIdx);
                disturbances(i,:,:) = solMax(u-1).omega;
                break
            end
            solMax(u) = Graph_maximization(G_init,G_j,P,solMin(j),BigOmega(i), S0);
        end
        clear solMax
        clear solMin
        %graph_plots(solMax(i), G_init, G_j, P, solMin.gamma);
    end

    %% Robust analysis
    for i=1:length(BigOmega)
        for j=1:1000
        [robustCompletion(i,j), ~] = graph_minimization_robust(G_init,G_j,P,S0, solOpt(i).gamma, solOpt(i).delta, BigOmega(i));
        end
    end
    boxplot(robustCompletion')

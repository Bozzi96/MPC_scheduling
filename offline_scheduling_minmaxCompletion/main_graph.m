%%% OFFLINE SCHEDULING
clear; clc; close all;
%%% Case study
P = [9 5 7 10 4 12
     4 7 3 7 1 10
     5 7 6 3 10 1
     4 3 10 6 4 5
     2 4 7 3 5 2
     1 6 5 3 6 8];

G_init = [1 2 3 4 6
           1 3 5 6 0
           1 2 3 4 5
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
G_j = [1 1 2 2 3 3 3 4 5 5 6 6 6]';
S0= [0 2 4 7 10 12]'; % planned
%S0 = [= [0 1 2 4 7 14]' ; % real
% Set planned release time and real release time
    % Pre processing of data
    [G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init, P);
    J = length(unique(G_j)); %jobs
    M = max(max(G)); %machines
    A = size(G_j,1);%alternatives
    D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints for each connection)
    

%% SOLVE PROBLEM
%    BigOmega =1:3:10;
BigOmega = 5;
    tStart = tic;
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
            if (j>=10)
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
    end
tEnd = toc(tStart);
    %% Robust analysis
%     for i=1:length(BigOmega)
%         for j=1:1000
%         [robustCompletion(i,j), ~] = graph_minimization_robust(G_init,G_j,P,S0, solOpt(i).gamma, solOpt(i).delta, BigOmega(i));
%         end
%     end
%%  Plot robust analysis 
boxplot(robustCompletion')
hold on
for i=1:length(BigOmega)
    plot(i,solOpt(i).C, "b*")
    %plot(i,min(robustCompletion(i,:)), "ko")
end
ylabel("Completion time")
xlabel("\Omega")
xticklabels({'1', '4','7','10'})
%Plot test
%% Get noise-free solution
figure()
graph_Gantt(sol_noNoise, G_init, G_j, P, sol_noNoise.gamma, M_init, "Noise-free solution");
xlabel("Time")
% Get best solution between noise-free and robust
i = randi(length(solOpt)); % Print random noise solution
i=4;
figure()
graph_Gantt(solOpt(i), G_init, G_j, P, solOpt(i).gamma, M_init, "Best solution (trade-off optimal/robust), \Omega= " + BigOmega(i));
xlabel("Time")
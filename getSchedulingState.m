function [startDates, endDates, path] =  getSchedulingState(sol, G_init, G_j, P, gamma, J_current)
[G,~, M_init, aux, aux_alt] = pre_processing_graph(G_init, P);
    %J_current = length(unique(G_j)); %jobs
    M = max(max(G)); %machines
    A = size(G_j,1);%alternatives
    D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints per each connection)
    
    mySol = G(gamma  > 0.1,:);
    mySol_init = mySol;
    for i=1:size(mySol,1)
        for j=1:size(mySol,2)
            if(mySol(i,j) > M_init)
                [mySol_init(i,j), ~] = find(mySol(i,j) == aux);
            end
        end
    end
    startDates = [];
    endDates = [];
    for j=1:J_current-1
        %TODO: filter to take into account only jobs that have already 
        % arrived in the shopfloor. J_current-1 is wrong!!!
         startDates(j,:)=sol.s(j,mySol(j,:));
         endDates(j,:)=sol.c(j,mySol(j,:));
    end
    path = G(gamma>0.1,:);
end

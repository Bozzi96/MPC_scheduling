function [startDates, endDates, path] =  getSchedulingState(sol, G_init, G_j, P, gamma, M_init)
[G,~, ~, aux, aux_alt] = pre_processing_graph(G_init, P, M_init);
    %J_current = length(unique(G_j)); %jobs
    %M = max(max(G)); %machines
    
    mySol = G(gamma  > 0.1,:);
    mySol_init = mySol;
    for i=1:size(mySol,1)
        for j=1:size(mySol,2)
            if(mySol(i,j) > M_init)
                [mySol_init(i,j), ~] = find(mySol(i,j) == aux);
            end
        end
    end
    startDates = {};
    endDates = {};
    for j=1:size(mySol,1)
        current_row = mySol(j,mySol(j,:)~=0); % Remove zero elements
        startDates{j} =sol.s(j,current_row);
        endDates{j}=sol.c(j,current_row);
    end
    path = G(gamma>0.1,:);
end

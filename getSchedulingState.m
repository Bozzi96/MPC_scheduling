function [startTime, completionTime, path] =  getSchedulingState(sol, G_init, G_j, P, gamma, M_init)
%%% Needed to accomplish dynamic scheduling --
%%% Get the state of the system in order to save the operations already
%%% performed and re-schedule jobs considering their past operations
[G,~, ~, aux, aux_alt] = pre_processing_graph(G_init, P, M_init);
    
    mySol = G(gamma  > 0.1,:);
    mySol_init = mySol;
    for i=1:size(mySol,1)
        for j=1:size(mySol,2)
            if(mySol(i,j) > M_init)
                [mySol_init(i,j), ~] = find(mySol(i,j) == aux);
            end
        end
    end
    startTime = {};
    completionTime = {};
    for j=1:size(mySol,1)
        current_row = mySol(j,mySol(j,:)~=0); % Remove zero elements
        startTime{j} =sol.s(j,current_row); % Get starting time
        completionTime{j}=sol.c(j,current_row); % Get completion time
    end
    path = G(gamma>0.1,:); % Save the path
end

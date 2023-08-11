function D = compute_D_from_graph(G,G_j)
% Find disjunctive connections from graph G    
D = 0;
for i=1:size(G,1)
    for j=1:size(G,2)
        if(G(i,j) == 0)
            continue;
        end
        m = G(i,j); %find machine m
        j_m = G_j(i); % job of machine m
        alt_other_j = find(G_j~=G_j(i)); % index of other jobs
        shared_machines = 0;
            shared_machines = sum(sum(ismember(G(alt_other_j,:),m)));
                    if(~isempty(shared_machines))
                        D = D+ shared_machines;
                    end
    end
end
end


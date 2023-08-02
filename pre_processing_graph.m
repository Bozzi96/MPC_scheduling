function [G_final, P_final, M, aux, aux_alt] = pre_processing_graph(G, P, M0)
    if nargin == 3
        M = M0;
    else
        M = max(max(G)); % find max number of machines
    end
    P_final = P; %initialize final processing matrix as the initial one
    delta = 1;
    aux = zeros(M,1);
    aux_alt = zeros(M,1);
    for i=1:size(G,1)
        for j=1:size(G,2)
            for k=j+1:size(G,2)
                if(G(i,j) == G(i,k) && G(i,j)~=0)
                    G(i,k) = M + delta;
                    idx = find(aux(G(i,j)) == 0, 1 ); % min(find(..))
                    idx_alt = find(aux_alt(G(i,j)) == 0, 1 ); % min(find(..))
                    if(~isempty(idx))
                        aux(G(i,j),idx) = M + delta;
                    else
                        aux(G(i,j),end+1) = M + delta;
                    end
                    if(~isempty(idx_alt))
                        aux_alt(G(i,j),idx_alt) = i; 
                    else
                        aux_alt(G(i,j),end+1) = i;
                    end
                    P(:,end+1) = P(:,G(i,j));
                    delta = delta + 1;
                end
            end
        end
        
    end
    G_final = G;
    P_final = P;
end
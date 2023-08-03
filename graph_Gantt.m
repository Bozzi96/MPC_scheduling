function graph_Gantt(sol, G_init, G_j, P, gamma, M0, graph_title)
    [G,~, M_init, aux, aux_alt] = pre_processing_graph(G_init, P, M0);
    J = length(unique(G_j)); %jobs
%     M = max(max(G)); %machines
%     A = size(G_j,1);%alternatives
%     D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints per each connection)
    
    mySol = G(gamma  > 0.1,:);
    mySol_init = mySol;
    for i=1:size(mySol,1)
        for j=1:size(mySol,2)
            if(mySol(i,j) > M_init)
                [mySol_init(i,j), ~] = find(mySol(i,j) == aux);
            end
        end
    end
    
%     t= 1:sol.C+1;
%     start_time = sol.s;
%     occupancy = zeros(M, length(t),J);
%     utilz = zeros (M_init,1);
    
%     Number_of_tasks=M_init;
    col = maxdistcolor(J, @sRGB_to_OKLab); % From matlab file exchange
    for j=1:J
        current_row = mySol(j,mySol(j,:)~=0); % Remove zero elements
        startDates{j} =sol.s(j,current_row);
        endDates{j}=sol.c(j,current_row);
%          startDates=sol.s(j,mySol(j,:));
%          endDates=sol.c(j,mySol(j,:));
         %width=sol.c(j,mySol(j,:))-sol.s(j,mySol(j,:));
         varname(j) = "job " + num2str(j); 
         i=1;
         for i=1:length(startDates{1,j})
             plot([startDates{1,j}(i),endDates{1,j}(i)],[mySol_init(j,i),mySol_init(j,i)],'b','Linewidth',10,'Color',col(j,:), DisplayName=varname(j))
            hold on
        end
    end
    title(graph_title)
   legendUnq();  % From matlab file exchange
   legend('-DynamicLegend')
    ylabel('Machine');
    ylim([0 M_init+1]);



%     for j=1:size(mySol,1)
%         Occupancy = [];
%         for m=1:size(mySol,2)
%             idx = mySol(j,m);
%             if(idx > 0 &&  ...
%                     start_time(j,idx)+P(j,idx) <= sol.C+1)
%                 
%                  processing = start_time(j,idx)+1:start_time(j,idx)+P(j,idx)+1;
%                 if(idx > M_init) 
%                     [row_m, ~] = find(aux == idx);
%                     occupancy(row_m, start_time(j,idx)+1:start_time(j,idx)+P(j,idx)+1,j) = 1;
%                     Occupancy = [Occupancy [processing;
%                                  row_m*ones(1,length(processing))]];
%                     utilz(row_m) = utilz(row_m) + P(j,row_m);
%                 else
%                     occupancy(idx, start_time(j,idx)+1:start_time(j,idx)+P(j,idx)+1,j) = 1;
%                    occupancy(idx, start_time(j,idx)+1:start_time(j,idx)+P(j,idx)+1,j) = 1;
%                                    Occupancy = [Occupancy [processing;
%                                  idx*ones(1,length(processing))]];
%                    utilz(idx) = utilz(idx) + P(j,idx);
%                 end
%             end
%         end
%         occ_final{j} = Occupancy;
%     end
%     figure();
%     xlabel('Time')
%     ylabel('Machine')
%     hold on
%     grid on
%     axis([0 sol.C+3 0 M_init+1])
%     for j=1:J
%         varname = "job " + num2str(j); 
%         plot(occ_final{1,j}(1,:),occ_final{1,j}(2,:), LineWidth=2, LineStyle="--", DisplayName=varname)
%        % ganttChart(occ_final{1,j}(2,:),duration([0 sol.C+3 0]))
%         legend('-DynamicLegend')
%     end
    
    %set(gca, 'ColorOrder', colormap(gray(J+1))) % black and white plot
    
    %% Plot the utilization rate of each machine
%     utilz = utilz/sol.C*100;
%     figure()
%     xlabel('Machine')
%      ylabel('Utilization rate (%)')
%     X = categorical({'M1','M2','M3','M4', 'M5'});
%     X = reordercats(X,{'M1','M2','M3','M4', 'M5'});
%     bar(X,utilz)
end
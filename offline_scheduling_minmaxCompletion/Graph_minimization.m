function sol = Graph_minimization(G,G_j,P, S0)% Parameters: 
    % G = graph 
    % G_j = number of alternatives (rows in the flow-shop graph)
    % P = matrix with processing time of job j on machine m (jobs x machines)
    % S0 = arrival time of jobs in the shop
    % SETS: (sets are re-computed from G and G_j automatically!)
    % J = jobs; M = machines; A = alternatives; D = disjunctive connections
    
    % Params
    BigM = 10000; % Big-M
    % Set computation
    G_init = G ;
    % Pre processing dei dati
    [G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init, P);
    J = length(unique(G_j)); %jobs
    M = max(max(G)); %machines
    A = size(G_j,1);%alternatives
    D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints per each connection)    
    
    % Optimization problem
    prob = optimproblem('ObjectiveSense','min');
    
    % Decision variables
    % s [j,m] = Start time of job j on machine m
    % c [j,m] = Completion time of job j on machine m
    % C = last completion time 
    % delta [D,1] = Disjunctive variables
    % gamma [A,1] = Choice variables
    s = optimvar('s', J, M, 'LowerBound', 0);
    c = optimvar('c', J, M, 'LowerBound', 0);
    C = optimvar('C', 1, 'LowerBound', 0);
    gamma = optimvar('gamma', A, 1, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
    delta = optimvar('delta', D, 1, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
    
    %%% Constraints %%%
    
    % Start time > S0
    cons_startTime = optimconstr(J, M);
    for j=1:J
        cons_startTime(j,:) = s(j,:) >= S0(j)*ones(1,M);
    end
    prob.Constraints.cons_startTime = cons_startTime;
    
    % Start time > Completion time previous machine conditioned to the choice
    % of that alternative in the graph
    cons_alternatives = optimconstr(sum(sum(G(:,2:end)~=0)),1);
    i = 1;
    for g1=1:size(G,1)
        for g2=2:size(G,2)
            if(G(g1,g2) ~=0)
            cons_alternatives(i) = s(G_j(g1),G(g1,g2)) >= c(G_j(g1),G(g1,g2-1)) - (1 - gamma(g1))*BigM;
            i = i+1;
            end
        end
    end
    prob.Constraints.cons_alternatives = cons_alternatives;
    
    % Completion time = start time + processing time on the same machine *
    % summation on alternatives in which job j passes through that machine
    cons_processingTime = optimconstr(J*M,1);
    i = 1;
    gamma_aux = 0;
    for j=1:J
        idx = find(G_j==j); % find alternatives of job j in G_j
        for m=1:M
            gamma_aux = 0;
            [idx_m,~] = find(G(:,:) == m); % find the row in which machine m is present among the rows of job j
            if(~isempty(idx_m))
                shared_idx = intersect(idx, idx_m); % find intersection between wors of job j and rowa in which machine m is present
                if(~isempty(shared_idx))
                    for k=1:size(shared_idx)
                        gamma_aux = gamma_aux+ gamma(shared_idx(k)); % add the alternatives
                    end
                    cons_processingTime(i) = c(j,m) == s(j,m) + P(j,m)*gamma_aux;
                else
                    cons_processingTime(i) = c(j,m) == s(j,m) + P(j,m);
                end
            else 
                cons_processingTime(i) = c(j,m) == s(j,m) + P(j,m);
            end
            i = i+1;
        end
    end
    prob.Constraints.cons_processingTime = cons_processingTime;
    
    % Disjunctive constraints 
    cons_disjunctive = optimconstr(D,1);
    idx_constraint = 1;
    idx_delta=1;
    aux_disj = zeros(D,6);
    for j=1:A
        other_jobs = find(G_j~=G_j(j)); % index of other jobs
        for i=1:length(other_jobs)
            idx_m_other_jobs = G(other_jobs(i),:); % index of machines of other jobs
                    shared_machines = intersect(idx_m_other_jobs, G(j,:));
                    shared_machines = shared_machines(shared_machines~=0); %remove zero from shared machines
                    for kkk=1:length(shared_machines)
                        if(sum(ismember(aux_disj,[G_j(j) G_j(other_jobs(i)) shared_machines(kkk) shared_machines(kkk) j other_jobs(i)], 'rows'))==0 && ...
                                sum(ismember(aux_disj,[G_j(other_jobs(i)) G_j(j) shared_machines(kkk) shared_machines(kkk) other_jobs(i) j], 'rows'))==0) % do not add two times the same disjunctive constraint
                                % If here, the double disjunctive constraint has not been added yet --> add it
                                cons_disjunctive(idx_constraint) = s(G_j(j),shared_machines(kkk)) >= (c(G_j(other_jobs(i)),shared_machines(kkk)) - (delta(idx_delta)*BigM ));
                                % aux_disj has the following structure: [job1 job2 machine1 alternative1 alternative2]
                                % it is needed to keep track of the constraints already inserted
                                aux_disj(idx_constraint,:) = [G_j(j) G_j(other_jobs(i)) shared_machines(kkk) shared_machines(kkk) j other_jobs(i) ];
                                idx_constraint = idx_constraint + 1;
                                cons_disjunctive(idx_constraint) = s(G_j(other_jobs(i)),shared_machines(kkk)) >= (c(G_j(j),shared_machines(kkk)) - ((1-delta(idx_delta))*BigM ));
                                idx_constraint = idx_constraint + 1;
                                idx_delta = idx_delta+1;
                        end
                    end
        end
    end
    
    cons_disjunctiveOnDuplicate = optimconstr(D,1);
    % Disjunctive constraints due to machine duplication
    for i=1:size(G,1)
        for j=1:size(G,2) 
            if(G(i,j)>M_init ) % if G(i,j) is a duplicated machine
                % Find the original machine
                [m_orig, col_orig] = find(G(i,j) == aux);
                % Find the alternative related to the duplicated machine and its job
                alt_m_orig = aux_alt(m_orig,col_orig);
                job_m_orig = G_j(alt_m_orig);
                % Find alternative in which m_orig is present
                [rows_alt, ~] = find(G_init == m_orig);
                rows_alt = unique(rows_alt);
                % If the alternative is different from the one of the
                % duplicated machine, and if the job is different --> add constraints
                rows_alt(G_j(rows_alt) == job_m_orig) = [];
                for a=1:length(rows_alt)
                    for aj=1:size(G,2)
                        if(G_init(i,j)==G_init(rows_alt(a),aj))
                            if(sum(ismember(aux_disj,[G_j(rows_alt(a)) job_m_orig G(i,j) G(rows_alt(a),aj) rows_alt(a) alt_m_orig], 'rows'))==0 && ...
                                sum(ismember(aux_disj,[job_m_orig G_j(rows_alt(a)) G(rows_alt(a),aj) G(i,j) alt_m_orig rows_alt(a)], 'rows'))==0) % non inserire due volte lo stesso vincolo disgiuntivo
                                cons_disjunctiveOnDuplicate(idx_constraint) = s(G_j(rows_alt(a)),G(rows_alt(a),aj)) >= (c(job_m_orig,G(i,j)) - (delta(idx_delta)*BigM ));%* (1-gamma(j))*BigM;
                                aux_disj(idx_constraint,:) = [G_j(rows_alt(a)) job_m_orig G(i,j) G(rows_alt(a),aj) rows_alt(a) alt_m_orig ]; % matrice per tenere traccia dei vincoli già aggiunti
                                idx_constraint = idx_constraint + 1;
                                cons_disjunctiveOnDuplicate(idx_constraint) = s(job_m_orig,G(i,j)) >= (c(G_j(rows_alt(a)),G(rows_alt(a),aj)) - ((1-delta(idx_delta))*BigM ));%* (1-gamma(j))*BigM;
                                idx_constraint = idx_constraint + 1;
                                idx_delta = idx_delta+1; 
                            end
                        end
                    end
                end
            end
        end
    end
    prob.Constraints.cons_disjunctive = cons_disjunctive;
    prob.Constraints.cons_disjunctiveOnDuplicate = cons_disjunctiveOnDuplicate;
    % Final completion time constraints
    cons_completionTime = optimconstr(A,1);
    for j=1:A
        cons_completionTime(j) = C >= c(G_j(j),G(j,find(G(j,:)~=0, 1, 'last' ))); % find(..'last') = max(find(..))
    end
    
    prob.Constraints.cons_completionTime = cons_completionTime;
    
    % Decision variables constraints (gamma)
    cons_gamma = optimconstr(J,1);
    for j=1:J
        idx = find(G_j == j); % for each possible alternative on job j
        cons_gamma(j) = sum(gamma(idx)) == 1; % choose only one alternative
    end
    
    prob.Constraints.cons_gamma = cons_gamma;
    
    % Cost function
    prob.Objective = C;%+sum(sum(s))+sum(sum(c));
    
    % Initial conditions
    x0.gamma = zeros(A,1);
    x0.delta = zeros(D,1);
    x0.C = 0;
    x0.c = zeros(J,M);
    x0.s = zeros(J,M);
    %% Solve problem
    %show(prob)
    tic
    [sol,val] = solve(prob,x0);
    toc
end
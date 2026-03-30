function T_opt = solveOptimization_CEI(T0, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Dataset_massflowrate, Dataset_cost, params)
% SOLVEOPTIMIZATION_CEI Global search using GA for CEI formulation.

    % Define the objective function handle specific to CEI
    objFun = @(T) objectiveFun_CEI(T, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Dataset_massflowrate, Dataset_cost, params);

    % Torque bounds
    T_nom = params.T_nom(:).';
    lb = T_nom + params.T_lower_offset(:).';
    ub = T_nom + params.T_upper_offset(:).';
    dim = numel(lb);

    % EA hyperparameters
    npop        = params.eaPopSize;         
    ngen        = params.eaNumGenerations;  
    cr          = params.eaCrossoverRate;   
    mr          = params.eaMutationRate;    
    sigma_mut   = params.eaMutationStd;     
    tournamentK = params.eaTournamentSize;  

    % Initialize population
    pop = zeros(npop, dim);
    pop(1,:) = min(max(T0(:).', lb), ub); 
    for i = 2:npop
        pop(i,:) = lb + rand(1,dim).*(ub - lb);
    end

    fitness = zeros(npop,1);
    for i = 1:npop, fitness(i) = objFun(pop(i,:)); end

    % Evolution loop
    for g = 1:ngen
        newpop = pop;
        for i = 1:npop
            idx1 = randi(npop, [tournamentK, 1]); [~, b1] = min(fitness(idx1)); p1 = pop(idx1(b1),:);
            idx2 = randi(npop, [tournamentK, 1]); [~, b2] = min(fitness(idx2)); p2 = pop(idx2(b2),:);

            child = p1;
            if rand < cr, child = rand*p1 + (1-rand)*p2; end
            if rand < mr, child = child + sigma_mut*randn(1,dim); end

            newpop(i,:) = min(max(child, lb), ub);
        end

        [~, bestIdxOld] = min(fitness);
        newpop(1,:) = pop(bestIdxOld,:); % Elitism

        pop = newpop;
        for i = 1:npop, fitness(i) = objFun(pop(i,:)); end
    end

    [~, bestIdx] = min(fitness);
    T_opt = pop(bestIdx,:);
end
function T_opt = solveOptimization_UCB(T0, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, params)
% SOLVEOPTIMIZATION Global search using a custom Evolutionary Algorithm (EA).
%
% Evolves a population of torque vectors to minimize the composite 
% acquisition function over the GP surrogate models.

    objFun = @(T) objectiveFun_UCB(T, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, params);

    T_nom = params.T_nom(:).';
    loff  = params.T_lower_offset(:).';
    uoff  = params.T_upper_offset(:).';
    lb = T_nom + loff;
    ub = T_nom + uoff;
    dim = numel(lb);

    npop        = params.eaPopSize;         
    ngen        = params.eaNumGenerations;  
    cr          = params.eaCrossoverRate;   
    mr          = params.eaMutationRate;    
    sigma_mut   = params.eaMutationStd;     
    tournamentK = params.eaTournamentSize;  

    pop = zeros(npop, dim);
    pop(1,:) = min(max(T0(:).', lb), ub); % First individual = current torque

    for i = 2:npop
        pop(i,:) = lb + rand(1,dim).*(ub - lb);
    end

    fitness = zeros(npop,1);
    for i = 1:npop
        fitness(i) = objFun(pop(i,:));
    end

    % --- Evolution loop ---
    for g = 1:ngen
        newpop = pop;
        for i = 1:npop

            idx1 = randi(npop, [tournamentK, 1]);
            [~, best1] = min(fitness(idx1));
            parent1 = pop(idx1(best1),:);
            idx2 = randi(npop, [tournamentK, 1]);
            [~, best2] = min(fitness(idx2));
            parent2 = pop(idx2(best2),:);

            child = parent1;
            if rand < cr
                alpha = rand; 
                child = alpha*parent1 + (1-alpha)*parent2;
            end

            if rand < mr
                child = child + sigma_mut*randn(1,dim);
            end

            newpop(i,:) = min(max(child, lb), ub);
        end

        [~, bestIdxOld] = min(fitness);
        newpop(1,:) = pop(bestIdxOld,:);

        pop = newpop;
        for i = 1:npop
            fitness(i) = objFun(pop(i,:));
        end
    end

    [~, bestIdx] = min(fitness);
    T_opt = pop(bestIdx,:);
end
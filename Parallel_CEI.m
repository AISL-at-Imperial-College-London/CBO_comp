function [sys,x0,str,ts] = Parallel_CEI(t,x,u,flag,params)
% parallel_cei_2
%
% MATLAB Level-1 S-function for real-time load-sharing optimization (RTO)
% of a 3-compressor parallel network using Constrained Expected Improvement (CEI).

persistent Dataset_T Dataset_Y ...
           GP_Power GP_Pdis GP_ds1 GP_ds2 GP_ds3...
           Current_setpoint DOE_points ParamsInit

if nargin < 5, params = struct(); end

switch flag
    case 0  % Initialization
        ParamsInit = initParams(params);
        [sys,x0,str,ts] = mdlInitializeSizes(ParamsInit);
        
        Dataset_T = []; Dataset_Y = [];
        GP_Power = []; GP_Pdis = []; GP_ds1 = []; GP_ds2 = []; GP_ds3 = [];
        
        Current_setpoint = ParamsInit.T_nom(:).';
        DOE_points = buildDOEpoints(ParamsInit);
        
    case 3  % Outputs
        [T_next, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Current_setpoint] = ...
            mdlOutputsCore(u, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Current_setpoint, DOE_points, ParamsInit);
        sys = T_next(:); 
        
    case {1, 2, 4, 9}
        sys = []; x0 = []; str = []; ts = [];
    otherwise
        sys = []; x0 = []; str = []; ts = [];
end
end 

%% ------------------------------------------------------------------------
function [sys,x0,str,ts] = mdlInitializeSizes(params)
    sizes = simsizes;
    sizes.NumContStates = 0; sizes.NumDiscStates = 0;
    sizes.NumOutputs = 3; sizes.NumInputs = 10;   
    sizes.DirFeedthrough = 1; sizes.NumSampleTimes = 1;
    sys = simsizes(sizes); x0 = []; str = [];
    ts  = [params.sampleTime 0];
end

%% ------------------------------------------------------------------------
function DOE_points = buildDOEpoints(params)
    T0 = params.T_nom(:).';   
    DOE_points = [T0; T0 + [0 0.25 0.25]; T0 + [0.25 0.25 0]];
    for i = 1:size(DOE_points,1)
        DOE_points(i,:) = enforceTorqueBounds(DOE_points(i,:), params);
    end
end

%% ------------------------------------------------------------------------
function [T_next, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Current_setpoint] = ...
    mdlOutputsCore(u, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Current_setpoint, DOE_points, params)
    
    u = u(:);
    P_tot = u(1); massflow_meas = u(2);
    ds1_meas = u(3); ds2_meas = u(4); ds3_meas = u(5);
    T_cur = [u(6), u(7), u(8)];
    massflow_ref = u(9); SS_flag = u(10);  
    
    T_next = Current_setpoint;
    dataAdded = false;
    nPoints = size(Dataset_T, 1);

    if SS_flag >= 0.5
        if nPoints == 0
            Dataset_T = T_cur;
            Dataset_Y = [P_tot, massflow_meas, ds1_meas, ds2_meas, ds3_meas];
            dataAdded = true;
        else
            if min(sqrt(sum((Dataset_T - T_cur).^2, 2))) > params.NoveltyRadius
                Dataset_T = [Dataset_T; T_cur];
                Dataset_Y = [Dataset_Y; P_tot, massflow_meas, ds1_meas, ds2_meas, ds3_meas];
                dataAdded = true;
            end
        end
        nPoints = size(Dataset_T, 1);

        if dataAdded && nPoints >= params.MinPointsForGP
            D = struct('U', Dataset_T, 'cost', Dataset_Y(:,1), 'constr', Dataset_Y(:,2:5));
            [GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3] = fit_all_gps(D, params.T_nom + params.T_lower_offset, params.T_nom + params.T_upper_offset);
        end

        if nPoints < params.MinPointsForGP
            T_next = chooseNextDOEPoint(DOE_points, Dataset_T, params);
        elseif ~isempty(GP_Power)
            % Pass Dataset constraint (massflow) and cost (power) to CEI solver
            T_next = solveOptimization_CEI(Current_setpoint, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Dataset_Y(:,2), Dataset_Y(:,1), params);
        end
    end
    
    T_next = enforceTorqueBounds(T_next, params);
    Current_setpoint = T_next;
end

%% ------------------------------------------------------------------------
function T_next = chooseNextDOEPoint(DOE_points, Dataset_T, params)
    if isempty(Dataset_T), T_next = DOE_points(1,:); return; end
    nDOE = size(DOE_points,1); used = false(nDOE,1);
    for i = 1:nDOE
        if min(sqrt(sum((Dataset_T - DOE_points(i,:)).^2, 2))) < params.NoveltyTol, used(i) = true; end
    end
    idx = find(~used, 1, 'first');
    if ~isempty(idx)
        T_next = DOE_points(idx,:);
    else
        jitter = params.dT_init * (rand(1, numel(params.T_nom)) - 0.5); 
        T_next = params.T_nom(:).' + jitter;
    end
end

%% ------------------------------------------------------------------------
function T_sat = enforceTorqueBounds(T, params)
    T = T(:).';  
    T_lower = params.T_nom(:).' + params.T_lower_offset(:).';
    T_upper = params.T_nom(:).' + params.T_upper_offset(:).';
    T_sat = min(max(T, T_lower), T_upper);
end

%% ------------------------------------------------------------------------
function p = initParams(p)
    if ~isfield(p, 'sampleTime'), p.sampleTime = 50; end
    if ~isfield(p, 'T_nom'), p.T_nom = [15; 16; 18]; end
    if ~isfield(p, 'T_lower_offset'), p.T_lower_offset = [-6; -7; -9]; end
    if ~isfield(p, 'T_upper_offset'), p.T_upper_offset = [ 6;  5;  3]; end
    if ~isfield(p, 'dT_init'), p.dT_init = 3.0; end
    if ~isfield(p, 'MinPointsForGP'), p.MinPointsForGP = 3; end
    if ~isfield(p, 'NoveltyRadius'), p.NoveltyRadius = 0.25; end
    if ~isfield(p, 'NoveltyTol'), p.NoveltyTol = 1e-3; end
    
    % CEI Safety & Tracking
    if ~isfield(p, 'surge_safe_margin'), p.surge_safe_margin = 0.01; end 
    if ~isfield(p, 'betaSafe'), p.betaSafe = 3; end 
    if ~isfield(p, 'cei_track_tol'), p.cei_track_tol = 0.001; end 

    % EA Config
    if ~isfield(p, 'eaPopSize'), p.eaPopSize = 25; end
    if ~isfield(p, 'eaNumGenerations'), p.eaNumGenerations = 20; end
    if ~isfield(p, 'eaCrossoverRate'), p.eaCrossoverRate = 0.6; end
    if ~isfield(p, 'eaMutationRate'), p.eaMutationRate = 0.4; end
    if ~isfield(p, 'eaMutationStd'), p.eaMutationStd = 0.5; end
    if ~isfield(p, 'eaTournamentSize'), p.eaTournamentSize = 2; end
end
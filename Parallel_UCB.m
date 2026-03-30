function [sys,x0,str,ts] = Parallel_UCB(t,x,u,flag,params)
% Parallel_ucb
%
% MATLAB Level-1 S-function for real-time load-sharing optimization (RTO)
% of a 3-compressor parallel network using Gaussian Process surrogates 
% and an evolutionary algorithm.
%
% INPUTS (u):
%   u(1) : P_tot         (total power consumption)
%   u(2) : massflow_meas (measured massflowrate)
%   u(3) : ds1_meas      (distance to surge compressor 1)
%   u(4) : ds2_meas      (distance to surge compressor 2)
%   u(5) : ds3_meas      (distance to surge compressor 3)
%   u(6) : T1_meas       (measured torque 1)
%   u(7) : T2_meas       (measured torque 2)
%   u(8) : T3_meas       (measured torque 3)
%   u(9) : massflow_ref  (reference massflowrate)
%   u(10): SS_flag       (steady-state indicator: 0 = no SS, 1 = SS)
%
% OUTPUTS:
%   sys(1) : T1_cmd (commanded torque 1)
%   sys(2) : T2_cmd (commanded torque 2)
%   sys(3) : T3_cmd (commanded torque 3)

persistent Dataset_T Dataset_Y ...
           GP_Power GP_Pdis GP_ds1 GP_ds2 GP_ds3...
           Current_setpoint DOE_points ParamsInit

if nargin < 5
    params = struct();  
end

switch flag
    case 0  % Initialization
        ParamsInit = initParams(params);
        [sys,x0,str,ts] = mdlInitializeSizes(ParamsInit);
        
        Dataset_T = [];
        Dataset_Y = [];
        GP_Power  = []; GP_Pdis = []; GP_ds1 = []; GP_ds2 = []; GP_ds3 = [];
        
        Current_setpoint = ParamsInit.T_nom(:).';
        DOE_points = buildDOEpoints(ParamsInit);
        
    case 3  % Outputs
        [T_next, Dataset_T, Dataset_Y, ...
            GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Current_setpoint] = ...
            mdlOutputsCore(u, Dataset_T, Dataset_Y, ...
                           GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3,...
                           Current_setpoint, DOE_points, ParamsInit);
        sys = T_next(:); 
        
    case {1, 2, 4, 9} % Unused flags
        sys = []; x0 = []; str = []; ts = [];
    otherwise
        sys = []; x0 = []; str = []; ts = [];
end
end 

%% ------------------------------------------------------------------------
function [sys,x0,str,ts] = mdlInitializeSizes(params)
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 3;   
    sizes.NumInputs      = 10;   
    sizes.DirFeedthrough = 1;   
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = []; str = [];
    ts  = [params.sampleTime 0];
end

%% ------------------------------------------------------------------------
function DOE_points = buildDOEpoints(params)
    T0 = params.T_nom(:).';   
    DOE_points = [
        T0 + [0 0 0];
        T0 + [0 0.25 0.25];
        T0 + [0.25 0.25 0]];
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
        % Novelty Detection
        if nPoints == 0
            Dataset_T = T_cur;
            Dataset_Y = [P_tot, massflow_meas, ds1_meas, ds2_meas, ds3_meas];
            dataAdded = true;
        else
            dists = sqrt(sum((Dataset_T - T_cur).^2, 2));
            if min(dists) > params.NoveltyRadius
                Dataset_T = [Dataset_T; T_cur];
                Dataset_Y = [Dataset_Y; P_tot, massflow_meas, ds1_meas, ds2_meas, ds3_meas];
                dataAdded = true;
            end
        end
        nPoints = size(Dataset_T, 1);

        % Train GPs
        if dataAdded && nPoints >= params.MinPointsForGP
            D = struct('U', Dataset_T, 'cost', Dataset_Y(:,1), 'constr', Dataset_Y(:,2:5));
            [GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3] = fit_all_gps(D, params.T_nom + params.T_lower_offset, params.T_nom + params.T_upper_offset);
            
            if isfield(params, 'debugMode') && params.debugMode
                assignin('base','RTO_Dataset', D);
            end
        end

        % Decide Next Action
        if nPoints < params.MinPointsForGP
            T_next = chooseNextDOEPoint(DOE_points, Dataset_T, params);
        elseif ~isempty(GP_Power)
            T_next = solveOptimization_UCB(Current_setpoint, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, params);
        end
    end
    
    T_next = enforceTorqueBounds(T_next, params);
    Current_setpoint = T_next;
end

%% ------------------------------------------------------------------------
function T_next = chooseNextDOEPoint(DOE_points, Dataset_T, params)
    if isempty(Dataset_T), T_next = DOE_points(1,:); return; end
    nDOE = size(DOE_points,1);
    used = false(nDOE,1);
    
    for i = 1:nDOE
        dists = sqrt(sum((Dataset_T - DOE_points(i,:)).^2, 2));
        if min(dists) < params.NoveltyTol, used(i) = true; end
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
    if ~isfield(p, 'T_nom'), p.T_nom = [16; 16; 16]; end
    if ~isfield(p, 'T_lower_offset'), p.T_lower_offset = [-5; -5; -5]; end
    if ~isfield(p, 'T_upper_offset'), p.T_upper_offset = [ 5;  5;  5]; end
    if ~isfield(p, 'dT_init'), p.dT_init = 3.0; end
    if ~isfield(p, 'MinPointsForGP'), p.MinPointsForGP = 3; end
    if ~isfield(p, 'NoveltyRadius'), p.NoveltyRadius = 0.25; end
    if ~isfield(p, 'NoveltyTol'), p.NoveltyTol = 1e-3; end
    
    % Core objective weights injected from paper logic
    if ~isfield(p, 'w_power'), p.w_power = 0.01; end
    if ~isfield(p, 'w_pdis'), p.w_pdis = 2000000; end
    if ~isfield(p, 'w_surge'), p.w_surge = 0.01; end
    if ~isfield(p, 'KappaUCB'), p.KappaUCB = 1.0; end % Note: Usually > 0 for exploration
    if ~isfield(p, 'surge_safe_margin'), p.surge_safe_margin = 0.01; end
    if ~isfield(p, 'betaSafe'), p.betaSafe = 3; end
    if ~isfield(p, 'debugMode'), p.debugMode = false; end

    % EA Config
    if ~isfield(p, 'eaPopSize'), p.eaPopSize = 25; end
    if ~isfield(p, 'eaNumGenerations'), p.eaNumGenerations = 20; end
    if ~isfield(p, 'eaCrossoverRate'), p.eaCrossoverRate = 0.6; end
    if ~isfield(p, 'eaMutationRate'), p.eaMutationRate = 0.4; end
    if ~isfield(p, 'eaMutationStd'), p.eaMutationStd = 0.5; end
    if ~isfield(p, 'eaTournamentSize'), p.eaTournamentSize = 2; end
end
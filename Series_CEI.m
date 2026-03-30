function [sys,x0,str,ts] = Series_CEI(t,x,u,flag,params)
% series_cei
%
% MATLAB Level-1 S-function for real-time load-sharing optimization (RTO)
% of a 2-compressor series network using Constrained Expected Improvement (CEI).

persistent Dataset_T Dataset_Y ...
           GP_Power GP_Pdis GP_ds1 GP_ds2 ...
           Current_setpoint DOE_points ParamsInit

if nargin < 5, params = struct(); end

switch flag
    case 0  
        ParamsInit = initParams(params);
        [sys,x0,str,ts] = mdlInitializeSizes(ParamsInit);
        Dataset_T = []; Dataset_Y = [];
        GP_Power = []; GP_Pdis = []; GP_ds1 = []; GP_ds2 = []; 
        Current_setpoint = ParamsInit.T_nom(:).';
        DOE_points = buildDOEpoints(ParamsInit);
        
    case 3  
        [T_next, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Current_setpoint] = ...
            mdlOutputsCore(u, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Current_setpoint, DOE_points, ParamsInit);
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
    sizes.NumOutputs = 2; sizes.NumInputs = 8;   
    sizes.DirFeedthrough = 1; sizes.NumSampleTimes = 1;
    sys = simsizes(sizes); x0 = []; str = [];
    ts  = [params.sampleTime 0];
end

%% ------------------------------------------------------------------------
function DOE_points = buildDOEpoints(params)
    T0 = params.T_nom(:).';   
    DOE_points = [T0; T0 + [0 0.3]; T0 + [0.3 0]];
    for i = 1:size(DOE_points,1)
        DOE_points(i,:) = enforceTorqueBounds(DOE_points(i,:), params);
    end
end

%% ------------------------------------------------------------------------
function [T_next, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Current_setpoint] = ...
    mdlOutputsCore(u, Dataset_T, Dataset_Y, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Current_setpoint, DOE_points, params)
    
    u = u(:);
    P_tot = u(1); P_dis_meas = u(2);
    ds1_meas = u(3); ds2_meas = u(4); 
    T_cur = [u(5), u(6)];
    P_dis_ref = u(7); SS_flag = u(8);  
    
    T_next = Current_setpoint;
    dataAdded = false;
    nPoints = size(Dataset_T, 1);

    if SS_flag >= 0.5
        if nPoints == 0
            Dataset_T = T_cur;
            Dataset_Y = [P_tot, P_dis_meas, ds1_meas, ds2_meas];
            dataAdded = true;
        else
            if min(sqrt(sum((Dataset_T - T_cur).^2, 2))) > params.NoveltyRadius
                Dataset_T = [Dataset_T; T_cur];
                Dataset_Y = [Dataset_Y; P_tot, P_dis_meas, ds1_meas, ds2_meas];
                dataAdded = true;
            end
        end
        nPoints = size(Dataset_T, 1);

        if dataAdded && nPoints >= params.MinPointsForGP
            D = struct('U', Dataset_T, 'cost', Dataset_Y(:,1), 'constr', Dataset_Y(:,2:4));
            [GP_Power, GP_Pdis, GP_ds1, GP_ds2] = fit_all_gps_series(D, params.T_nom + params.T_lower_offset, params.T_nom + params.T_upper_offset);
        end

        if nPoints < params.MinPointsForGP
            T_next = chooseNextDOEPoint(DOE_points, Dataset_T, params);
        elseif ~isempty(GP_Power)
            T_next = solveOptimization_CEI_series(Current_setpoint, P_dis_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Dataset_Y(:,2), Dataset_Y(:,1), params);
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
        T_next = params.T_nom(:).' + params.dT_init * (rand(1, 2) - 0.5); 
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
    if ~isfield(p, 'T_nom'), p.T_nom = [16.5; 16.5]; end
    if ~isfield(p, 'T_lower_offset'), p.T_lower_offset = [-5.5; -5.5]; end
    if ~isfield(p, 'T_upper_offset'), p.T_upper_offset = [ 4.5;  4.5]; end
    if ~isfield(p, 'dT_init'), p.dT_init = 3.0; end
    if ~isfield(p, 'MinPointsForGP'), p.MinPointsForGP = 3; end
    if ~isfield(p, 'NoveltyRadius'), p.NoveltyRadius = 0.25; end
    if ~isfield(p, 'NoveltyTol'), p.NoveltyTol = 1e-3; end
    
    if ~isfield(p, 'surge_safe_margin'), p.surge_safe_margin = 0.01; end % 0.01 as per series CEI
    if ~isfield(p, 'betaSafe'), p.betaSafe = 3; end 
    if ~isfield(p, 'cei_track_tol'), p.cei_track_tol = 1; end 

    if ~isfield(p, 'eaPopSize'), p.eaPopSize = 25; end
    if ~isfield(p, 'eaNumGenerations'), p.eaNumGenerations = 20; end
    if ~isfield(p, 'eaCrossoverRate'), p.eaCrossoverRate = 0.6; end
    if ~isfield(p, 'eaMutationRate'), p.eaMutationRate = 0.4; end
    if ~isfield(p, 'eaMutationStd'), p.eaMutationStd = 0.5; end
    if ~isfield(p, 'eaTournamentSize'), p.eaTournamentSize = 2; end
end
function [sys,x0,str,tss]=comp_plant(t,x,u,flag,Param,X_ss)

switch flag,

case 0,	% Initialize the states and sizes
    [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
    % ****************
  	% Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states

   sys = mdlDerivatives(t,x,u,Param);

otherwise,
   sys = [];
end

% ******************************************
% Sub-routines or Functions
% ******************************************

% ******************************************
% Initialization
% ******************************************

function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);

% This handles initialization of the function.
% Call simsize of a sizes structure.

sizes = simsizes;
sizes.NumContStates  = 6;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 11;     % outputs of model (add eta, T_out, POWER)
sizes.NumInputs      = 6;     % inputs of model (add T)
sizes.DirFeedthrough = 1;     % System is causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the states 


str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************
% Outputs
% ******************************************

function [sys] = mdlOutputs(t,x,u,Param);
global Out_pres m_out pr2

Valve_in_gain = 1/(sqrt(1.05e5-1e5));   % opening gain
Valve_rec_gain = 1/(sqrt(1.2891e5-1e5)); % opening gain
Valve_out_gain = 1/(sqrt(1.2891e5-1.2e5)); % opening gain

% Inputs

torque_drive = u(1); % 
Inflow_opening = u(2); % 
Outflow_opening = u(3); % 
Recycle_opening = u(4); % 

P_int    = u(5);  % input pressure

T_in = u(6); % question:do i need to put T_in as a step function as the sixth input as the paraellel one?
% States

p3 = x(1);%
p4 = x(2);%
m_comp = x(3);%
pr2 = x(4);%
omega_comp = x(5);%
m_rec = x(6);%

% Parameters

SpeedSound = 340;     % speed of sound
% In_pres    = 1.05e5;  % input pressure  
% Out_pres   = u(5); % output pressure


Out_pres   = 2.35e5; % output pressure


Valve_in_gain = 1/(sqrt(1.05e5-1e5));   % opening gain a
Valve_rec_gain = 1/(sqrt(1.2891e5-1e5)); % opening gain
Valve_out_gain = 1/(sqrt(1.2891e5-1.2e5)); % opening gain
 %In_pres
m_in = Valve_in_gain*Inflow_opening*sqrt(abs(P_int - p3))*(P_int - p3)/abs(P_int - p3);% Inflow valve
m_out = Valve_out_gain*Outflow_opening*sqrt(abs(p4 - Out_pres))*(p4 - Out_pres)/abs(p4 - Out_pres);% Outflow valve
m_out = max(m_out, 0);
m_rec_ss = Valve_rec_gain*Recycle_opening*sqrt(abs(p4 - p3))*(p4 - p3)/abs(p4 - p3);% Recycle valve

% newly add about the efficiency line (start)

% Step 1: 计算 omega_rpm 和百分比
omega_rpm = omega_comp / (2*pi/60);
omega_percent = omega_rpm / 8370 * 100;  % 注意8370为额定转速，如有不同需替换

% Step 2: 用多项式拟合压力比
%pr_coeff = [2.690940401290303  -0.013878128060951  -0.040925719808930   0.000986961896765  -0.000418575028867   0.000024527875520];
%pr2 = pr_coeff * [1; m_comp; omega_percent; m_comp*omega_percent; m_comp*m_comp; omega_percent*omega_percent];
pr2=x(4);


% Step 3: 用多项式拟合效率
eta = 0.90 - (5.5*(m_comp - 0.52).^2 - 1.8*(m_comp - 0.52).*(pr2 - 1.7) + 1.5*(pr2 - 1.7).^2);


%if ~isreal(eta)
    % clear eta
    % eta = 0.00001
% end

% Step 4: 计算 T_out 和 Power（假设 u(5) 是 T_in）
kappa = 1.27;
T_out = u(6) * (pr2)^((kappa - 1)/(kappa * eta));

if ~isreal(T_out) || isnan(T_out) || isinf(T_out)
    T_out = 300;  % 设置一个安全默认值
end

%if ~isreal(T_out)
 %    clear T_out
  %   T_out = 0;
%end

% 单位质量所需功
yp = ((0.9 * 8314 * (20+273))/16.04)*(kappa / (kappa - 1))*(pr2^((kappa - 1)/kappa) - 1);

% 总功率
POWER = yp / eta * m_comp;

if ~isreal(POWER) || isnan(POWER) || isinf(POWER)
    POWER = 100;  % 或者其他默认值
end

% newly add about the efficiency line (end)
surge_d = 100;

sys(1) = x(1); % p3
sys(2) = x(2); % p4
sys(3) = x(3); % m_comp 
sys(4) = x(4); % pr2
sys(5) = x(5); % omega_comp
sys(6) = m_in; % 
sys(7) = m_out; % 
sys(8) = m_rec; % 

sys(9)  = eta;    % Step 5: 将这三个新变量加入输出
sys(10) = T_out;
sys(11) = POWER;
% ******************************************
% Derivatives
% ******************************************

function sys = mdlDerivatives(t,x,u,Param)

% Inputs

torque_drive = u(1); % 
Inflow_opening = u(2); % 
Outflow_opening = u(3); % 
Recycle_opening = u(4); % 
P_int    = u(5);  % input pressure
T_in = u(6); % question:do i need to put T_in as a step function as the sixth input as the paraellel one?
% States

p3 = x(1);%
p4 = x(2);%
m_comp = x(3);%
pr2 = x(4);%
omega_comp = x(5);%
m_rec = x(6);%

% Parameters

SpeedSound = 340;     % speed of sound
P_int    = u(5);  % input pressure
Out_pres   = 2.35e5; % output pressure


Valve_in_gain = 1/(sqrt(1.05e5-1e5));   % opening gain
Valve_rec_gain = 1/(sqrt(1.2891e5-1e5)); % opening gain
Valve_out_gain = 1/(sqrt(1.2891e5-1.2e5)); % opening gain

VolumeT1 = 1e-2; % volume
VolumeT2 = 1e-3; % volume
AdivL = 0.2e-3;  % Duct cross section divided by length

%%% Compressor
Dt         = 0.18;  % impeller diameter 0.7230 m
sigma_slip = 0.9;   % compressor slip 0.9
J          = 0.022;  % shaft inertia ^just dynamic^
tauComp    = 8;     % time constant of the pressure ratio ^just dynamic^
tauRecycle = 1;     % time constant of the recycle valve ^just dynamic^

%%% coefficients of the compressor map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C0 =[0.000000000061147  -0.000000603072143   0.002205681837355   -1.461782281888261];
C1 =[-0.000000000882663   0.000008394635394  -0.026129863772614    25.840599079312330];
C2 =1.4*[0.000000001205096  -0.000008916645219   0.018455555634821   -1.277068107070351 ];
C3 =[0.0000000150689  -0.0001472249440   0.4791436514340  -529.0322734770704 ];



% Algebraic equations

% valves
m_in = Valve_in_gain*Inflow_opening*sqrt(abs(P_int - p3))*(P_int - p3)/abs(P_int - p3); % Inflow valve
m_out = Valve_out_gain*Outflow_opening*sqrt(abs(p4 - Out_pres))*(p4 - Out_pres)/abs(p4 - Out_pres); % Outflow valve
m_out = max(m_out, 0);
m_rec_ss = Valve_rec_gain*Recycle_opening*sqrt(abs(p4 - p3))*(p4 - p3)/abs(p4 - p3); % Recycle valve

% compressor pressure ratio
omega_compp = omega_comp*0.72;
A0=C0(4)+C0(3)*omega_compp+C0(2)*omega_compp^2+C0(1)*omega_compp^3;
A1=C1(4)+C1(3)*omega_compp+C1(2)*omega_compp^2+C1(1)*omega_compp^3;
A2=C2(4)+C2(3)*omega_compp+C2(2)*omega_compp^2+C2(1)*omega_compp^3;
A3=C3(4)+C3(3)*omega_compp+C3(2)*omega_compp^2+C3(1)*omega_compp^3;

p_ratio_ss=A0+A1*m_comp+A2*m_comp^2+A3*m_comp^3;

% compressor pressure
p_comp = pr2*p3;

% compressor torque
torque_comp = sigma_slip * (Dt/2)^2 * omega_comp * abs(m_comp); % Torque_comp (N.m) omega_comp (rpm) m_comp (kg/s)

% sys = zeros(6,1);

sys(1) = SpeedSound/VolumeT1*(m_in + m_rec - m_comp) ; % 
sys(2) = SpeedSound/VolumeT2*(m_comp - m_rec - m_out) ; % 
sys(3) = AdivL*(p_comp-p4); %
sys(4) = tauComp*(p_ratio_ss-pr2); % 
sys(5) = 1/J*(torque_drive-torque_comp); % look for compressor inertia in similar-size compressor)
sys(6) = tauRecycle*(m_rec_ss - m_rec);%



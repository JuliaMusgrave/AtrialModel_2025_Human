% This function runs the typical WL simulation given appropriate input
% parameters
% Inputs:
%       - model: function handle to the muscle model (should output k3 and
%           k-3). CURRENTLY NOT USED
%       - params: structure representing the input parameters for the
%           muscle model (should have M_frac element)
%       - c_input: 2xn array representing the time and calcium data to 
%           input into the model
% Outputs:
%       - WL: output structure which contains 'Energy', 'Work', 'shortening' and 
%           'velocity' per afterload (vectors); 'after' - the afterloads (vector); and 
%           'length', 'stress' and 'time' per afterload (cell arrays)

% Following parameters are defined in the code:
% starting SL (2.2), Temp (37), viscosity and mass for isotonic,
% params.mode (both sarcomere and loop modes used), ADP (0.036), XB
% density (0.25 mol/m^3), muscle model function (@Mmodel_2025_Human)

% Author: Julia Musgrave
% Date: June 2024

function [WL]=WL_function(model,params,c_input)

global isotonic;



% set up parameters

dias=min(c_input(2,:));
t=c_input(1,:);
L=2.2;

G_ATP0 = -30; % kJ/mol
R = 8.314/1000; % kJ/(mol.K)
T=37+273; % K (experiments performed at 37 degC)


% get ICs and run to diastolic SS
params.mode='sarcomere';
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);
y_dias(9)=0;


% start with isosarcometric Ca to get isometric peak

[t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_dias,options);
[~,F_twitch]=model(t,y,L,c_input,params);


afterloads=linspace((min(F_twitch)*1.05),max(F_twitch),10);
%_______________________________________________________
% run basic loop
params.mode='loop';
tspan=linspace(0,1,1000);

for loop=1:10
    a=afterloads(loop);
    if loop==10 % boost the afterload 1% past isometric just to make sure we don't get small shortening from numerical issues
        a=a*1.01;
    end
    params.loop=[a 0.00001 0.00005]; %afterload, mass, visc
    y_last(7)=1000;
    y_curr=y_dias;
    %run the sim until steady state on N
    while abs(y_last(7)-y_curr(7))>2e-5
        y_last=y_curr;
        % need a trigger for isotonic
        isotonic = false;
        [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),tspan,y_last,options);
        y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
    end

    isotonic = false;
    % in workloop mode, we have to go one point at a time to get force
    F=zeros(length(t),1);
    Fp=F; dL=F; k3=F; k_3=F;
    for i=1:length(t)    
        [dy, F(i), Fp(i),k3(i),k_3(i)] = model(t(i),y(i,:),L,c_input,params);
        dL(i) = dy(8);
    end
    
    % energy calculations
    ADP=0.036;
    G_ATP=G_ATP0+R*T*log((ADP/1000*params.met(2)/1000)/(params.met(1)/1000)); % convert mets to mol/L
    A=1 + (y(:,8)/2.3-1)*params.xb(8)-y(:,2)-y(:,1)-y(:,7);
    ATPase_rate=(k3.*y(:,2)-k_3.*A)*0.25*params.M_frac;
    WL.Energy(loop)=abs(trapz(t,ATPase_rate))*-G_ATP; % unitless?
    
    % other workloop analysis
    WL.Work(loop)=abs(trapz(y(:,8)/L,F)); % kPa.um?
    WL.shortening(loop)=100-min(y(:,8))/L*100; % percent of Lo
    %WL.velocity(loop)=-min(movmean(dL,100))/L; % muscle length /s 
    % finding max slop of the SL is much more reliable re noise
    WL.velocity(loop)=-min(gradient(y(:,8),0.001))/L; % muscle length /s 

    
    % save time cours data
    WL.length{loop}=y(:,8);
    WL.stress{loop}=F;
    WL.time{loop}=t;

end
WL.after=afterloads;



end

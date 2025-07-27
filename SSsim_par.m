%

% Runs a XB ODE function until force reaches a steady state. Called as part
% of many of my simulation protocols
%
% This version is currently configured to exit after 100 s of sim time
%
% Inputs:   
%       - fun: function handle
%       - tspan: tspan input required for ode15s
%       - y0: initial conditions
%       - s: vector of length perturbation information
%       - Ca: calcium concentration
%       - params: parameters needed to run the model
% Outputs:  
%       - t and y: vectors of the final iteration
%
% Author: Julia Musgrave
% Date last updated: Dec 2020 

function [t,y] = SSsim_par(fun, tspan, y0,s,Ca,params)
prev=0;
curr=100;
r=0;
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);
while abs((prev-curr)/curr)>1e-5 && r*tspan(end)<100
[t,y]=ode15s(@(t,y)fun(t,y,s,Ca,params),tspan,y0,options);
y0=y(end,:);
prev=curr;
[~,curr]=fun(t(end),y0,s,Ca,params);
r=r+1;
end

end
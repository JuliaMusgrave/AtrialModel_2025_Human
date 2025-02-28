%

% Function describing the final ODE cross-bridge model parameterised human
% atrial data
%
% Model is fully activated so has no Ca2+ dependence, also no passive force
% Should be run with x_p parameters.
%
% This is a different permutation of the final model presented in Musgrave
% et al. 2024.  XB strain equations based on Razumova et al. (1999) 3 state stiffness 
% distortion model. Sarcomere length dependence on available cross-bridge 
% proportion. Strain dependence on cross-bridge detachment (k3 rate) and
% attachment (k1 rate). Rapid equilibrium-based dependence on ATP only. 
% 
% Set up to be solved numerically under a range of sarcomere length inputs,
% including perturbation protocols.
%
% Inputs:
%       - t: time (seconds)
%       - y: states in the model 
%       - s: sarcomere length (um) input. Either a constant or 2 x n matrix 
%           with matching time and length points (t in 1st row and L in 2nd)
%       - params: variable parameters fed into the model
%           params = [k1 k-1' k2' k-2' k3' phi_x phi_v phi_l K Ks phi_s-2 phi_s3 K_ATP K_Pi];
%           no_params = 14
%       - met: 1x2 array containing ATP and Pi concentrations (mM) -
%       optional 5 mM ATP and 1 mM Pi used as baseline
% Outputs:
%       - dydt: system of odes to solve 
%       - F: active stress solution of the model (kPa)
%
% Author: Julia Musgrave
% Date : October 2024



function [dydt, F] = XBmodel_2025_Human(t,y,s,params,met)

if nargin==0 
   dydt=[0.1 0.1 0 0.01];
  return
end
% --------------------------
% getting sarcomere length and velocity
[L, ~, dL] = get_sim_inputs(t,s,33);

%----------------------------
% assigning the parameters 
k1=params(1); % s^-1
k_1=params(2); % s^-1
k2=params(3); % s^-1
k_2=params(4); % s^-1
k3=params(5); % s^-1
%phi_x=params(6); % unitless (scaling factor removed for human model)
phi_v=params(7);  % unitless 
phi_l=params(8); % unitless 
K=params(9); % GPa/m (kPa/um)
%Ks=params(10); % zero in human model
phi_s1=params(11); % um^-1
phi_s3=params(12); % um^-1
kd_ATP=params(13); % mM
%kd_Pi=params(14); % not present in human model

if nargin>5
ATP=met(1); % (mM) 
Pi=met(2); % (mM) 
else
    ATP=5;
    Pi=1;
end

% model constants
xC0=0.01;  % powerstroke size (um) 
Lmax=2.3;  % (um)

if size(y,2)==1; y=y'; end

%----------------------------------
% state variables
B=y(:,1);
C=y(:,2);
xB=y(:,3);
xC=y(:,4);

% algebraic equations (depending on state variables)
k1=k1*exp(-phi_s1*xB);
k3=k3*exp(phi_s3*(xC-xC0));
Z=1 + (L<Lmax).*(L/Lmax-1)*phi_l;
A=Z-B-C;

% incorporating metabolite dependence into rates
k_1=k_1*Pi;

k3=k3*ATP/(ATP+kd_ATP);
k_2=k_2*kd_ATP/(ATP+kd_ATP);

% thermodynamic constraint for k-3 (as in Tran et al., (2010))
G_ATP0 = -30; % kJ/mol
R = 8.314/1000; % kJ/(mol.K)
T=37+273; % K (experiments performed at 37 degC)
ADP=36e-6; % M
ATP=ATP/1000;
Pi=Pi/1000; %(M for thermo constraint)

G_ATP=G_ATP0+R*T*log((ADP*Pi)/ATP);
k_3 = k1 .* k2 .* k3 ./ (k_1 * k_2 * exp(-G_ATP/(R*T)));

% ODEs
dydt=zeros(size(y));
dydt(:,1) = k1.*A + k_2.*C -(k_1+k2).*B; %dB/dt
dydt(:,2) = k2*B -(k_2+k3).*C + k_3.*A; %dC/dt
dydt(:,3) = -1./B.*xB.*(A.*k1 + C.*k_2) + dL*phi_v; %dxB/dt
dydt(:,4) = -1./C.*(xC-xC0).*(B*k2+A.*k_3) + dL*phi_v; %dxC/dt

dydt=dydt(:);

% Active stress (force) equation 
F=K*(B.*xB+C.*xC); % (kPa)

end

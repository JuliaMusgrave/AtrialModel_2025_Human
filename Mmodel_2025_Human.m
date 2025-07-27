%

% Function describing an ODE muscle model for human atria 
% published by Musgrave et al in 2025.
%
% Ca2+ dependencies are based closely on those of Land et al. (2017) and
%
% XB strain equations based on Razumova et al. (1999) 3 state stiffeness 
% distortion model. Sarcomere length dependence on available cross-bridge 
% proportion. Strain dependence on cross-bridge detachment (k3 rate) and
% attachment (k1 rate).
% 
% Set up to be solved numerically under a range of sarcomere length inputs,
% including perturbation protocols.
%
% Inputs:
%       - t: time (seconds)
%       - y: states in the model 
%       - s: sarcomere length (um) input. Either a constant or 2 x n matrix 
%           with matching time and length points (t in 1st row and L in 2nd)
%       - c: intracellular [Ca2+] (uM). Either a constant or 2 x n matrix 
%           with matching time and Ca points (t in 1st row and Ca in 2nd)
%       - params: variable parameters fed into the model. Follows a
%       structure including: 
%       .xb : cross-bridge params (vary for diabetic)
%               = [k1 k-1 k2 k-2 k3 phi_x phi_v phi_l K [] phi_s1 phi_s3 kd_ATP []]
%       .pass : passive model params (vary for diabetic)
%               = [kP1 kP2 eta phi_e]
%       .ca : Ca-related parameters (based on Land model, except Ca50)
%               = [Ca50 n_Tm ktrpn n_trpn k0 k_0]
%       .met : concentration of ATP and Pi (in mM)
%               =[ATP Pi]
%       .mode : mode for the simulation
%               = either 'sarcomere' - for normal simulations with 
%               sarcomere control or 'loop'- for work loops
%       .loop : parameters for the work-loop simulation
%               = [Fafter mass visc]

% Outputs:
%       - dydt: system of odes to solve 
%       - F: total (active + passive) stress solution of the model (kPa)
%       - Fp: passive stress of the muscle (kPa)
%       - k3: ATPase rate during simulation (s^-1)
%       - k_3: reverse ATPase rate during simulation (s^-1)
%
% Author: Julia Musgrave
% Date : Sep 2024



function [dydt, F, Fp,k3,k_3] = Mmodel_2025_Human(t,y,s,c,params)

% ------------------ initialising state variables -----------------------
if nargin==0 
   dydt=[1e-4 1e-4 0 0.01 0 0.001 0.79 2.2 0];
  return
end

%---------------------------- parameters --------------------------------
try 
    mode= params.mode; %'sarcomere'; 'loop' 
catch
    mode = 'sarcomere'; % SL is controlled by default
end

xb_params=params.xb;
p_params=params.passive;
ca_params=params.ca;

% cross-bridge parameters
k1=xb_params(1); % s^-1
k_1=xb_params(2); % s^-1
k2=xb_params(3); % s^-1
k_2=xb_params(4); % s^-1
k3=xb_params(5); % s^-1
phi_x=xb_params(6); % unitless (= 1 for human model)
phi_v=xb_params(7);  % unitless 
phi_l=xb_params(8); % unitless 
K=xb_params(9); % GPa/m (kPa/um)
%Ks=xb_params(10); % zero in human model
phi_s1=xb_params(11); % um^-1
phi_s3=xb_params(12); % um^-1
kd_ATP=xb_params(13); % mM
%kd_Pi=xb_params(14); % not present in human model

% passive force parameters
kP1=p_params(1);
kP2=p_params(2);
eta=p_params(3);
phi_e=p_params(4);

% Ca regulation parameters 
Ca50=ca_params(1);
n_Tm=ca_params(2);
ktrpn=ca_params(3);
n_trpn=ca_params(4);
k0=ca_params(5);
k_0=ca_params(6);

% model constants
xC0=0.01;  % powerstroke size (um) 
Lmax=2.3;  % (um)
Lr=2.2*0.85;
ADP=0.036; % mM (based on Tran et al. 2010)

% input metabolite concentrations
ATP=params.met(1);
Pi=params.met(2);



% ------------ getting sarcomere length, velocity and Ca -----------------
if size(y,2)==1; y=y'; end
[L, L0, dL, Ca] = get_sim_inputs(t,s,c,y(:,8));
if mode(1)=='s' && length(s)==1  % isosarcometric contraction
    L=L0;
end

%------------------ state variables -------------------------------
% original 5
B=y(:,1);
C=y(:,2);
xB=y(:,3);
xC=y(:,4);
F1=y(:,5);

CaT=y(:,6);
N=y(:,7);

intF=y(:,9);

%XB_frac=(B+C)/(0.0137+0.0325);

% ------------------ cross bridge model ----------------------------------
% algebraic equations (depending on state variables)
k1=k1*exp(-phi_s1*xB);
k3=k3*exp(phi_s3*(xC-xC0));
Z=1 + (L<Lmax).*(L/Lmax-1)*phi_l;
A=max(0,Z-B-C-N);

% incorporating metabolite dependence into rates
k_1=k_1*Pi;
k3=k3*ATP/(ATP+kd_ATP);
k_2=k_2*kd_ATP/(ATP+kd_ATP);

% thermodynamic constraint for k-3 (as in Tran et al. (2010))
G_ATP0 = -30; % kJ/mol
R = 8.314/1000; % kJ/(mol.K)
T=37+273; % K (experiments performed at 37 degC)
G_ATP=G_ATP0+R*T*log((ADP/1000*Pi/1000)/(ATP/1000)); % convert mets to mol/L
k_3 = k1 .* k2 .* k3 ./ (k_1 * k_2 * exp(-G_ATP/(R*T)));

% ------------------ force calculations ----------------------------------

% Active stress (force) equation 
Fa=K*(B.*xB+C.*xC); % (kPa)


% Passive stress
F2= kP2*(exp(phi_e*(L-Lr))-1);
Fp=F1+F2;

% Total stress (force)
F=Fa+Fp;

Fafter=0;

% finding dL when length perturbation not imposed
if dL(1) == 'i' % 'iso'
        
    dL=0; % for isosarcometric case (including start of workloop)

    if mode(1)=='l'

        % funky stuff for loops

        Fafter=params.loop(1); % kPa
        mass=params.loop(2); % kPa.s^2/(um) = GPa.s^2/m
        visc=params.loop(3); % kPa.s/(um) = GPa.s/m
        global isotonic ;

        % Check isotonic switch
        if isotonic 
            dL = (intF + visc*(Lr - L))/mass;
            if dL > 0 % need to switch back to isometric
                dL=0;
            end
        end
        
        if t>0.8 && L<L0% restretch
            % start with a gradual acceleration
            dL = sin(pi*(t-0.85)/.1)+1;
            if t>0.85 % pick a constant velocity to reach starting SL by end of loop
                dL = (L0 - L)/(0.95-t);
                if t>0.95
                    dL=0;
                end
            end
        end
    
        % Prevent over extension in restretch
        if t>0.8 && isotonic && (L > (L0*1.0001))
            dL = 0;
        end

        if ~isotonic && F > Fafter
            isotonic = true;
        end

    end
    

end
% ----------------------- ODEs -----------------------------------------
dydt=zeros(size(y));

%troponin
dydt(:,6) = ktrpn*((Ca./Ca50).^n_trpn.*(1-CaT)-CaT); % dCaT/dt

% XB states
dydt(:,1) = k1.*A + k_2.*C -(k_1+k2).*B; % dB/dt
dydt(:,2) = k2*B -(k_2+k3).*C + k_3.*A; % dC/dt
dydt(:,7) = k_0*CaT.^(-n_Tm/2).*A-k0*CaT.^(n_Tm/2).*N; %dN/dt

% XB strain
dydt(:,3) = -phi_x./B.*xB.*(A.*k1 + C.*k_2) + dL*phi_v; % dxB/dt
dydt(:,4) = -phi_x./C.*(xC-xC0).*(B.*k2+A.*k_3) + dL*phi_v; % dxC/dt

% passive stress
dydt(:,5) = kP1*(dL-F1/eta);   % dF1/dt

% whole muscle mechanics
dydt(:,8) = dL; % dL/dt
dydt(:,9) = -Fa-Fp+Fafter;
if mode(1)=='l' && ~isotonic
    dydt(:,9) = 0;
end

dydt=dydt(:);




end


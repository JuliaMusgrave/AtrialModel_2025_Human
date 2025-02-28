% Function to simplify the length and or calcium inputs to an ode muscle
% model. Gets the right type of variable depending on the simulation
%
% Inputs:
%       - t: the t value from the ODE (either a scalar or vector)
%       - length_info: the variable that I typically call s in my ode
%       functions (sarcomere length input. Either a constant or 2 x n matrix 
%       with matching time and length points (t in 1st row and L in 2nd)
%       - ca_info: a similar variable but instead representing the
%       intracellular calcium concentration (either constant or the array
%       form)
%       - y8(optional) : the value of L found from the dL ode in models
%       with internal shortening
% Outputs:
%       - L: sarcomere length value or array (depending on form of t and length_info)
%       - L0: the initial sarcomere (which is only required if we are doing
%       internal shortening)
%       - dL: sarcomere velocity found numerically from length input, array
%       or scalar. Takes the string value 'iso' if it needs to
%       be computed when solving the model
%       - Ca: calcium concentration, array or scalar
% Author: Julia Musgrave
% Date last updated: 

function [L, L0, dL, Ca] = get_sim_inputs(t,length_info,ca_info,y8)

if length(length_info)>1 % length perturbation

    if length(t)==1 % in solver mode (step by step)
        % The ODE solver won't take the same steps as the length input so 
        % we need to interpolate to find the estimated length
        L=interp1(length_info(1,:),length_info(2,:),t);

        % Use finite difference formulae to get velocity
        dL_full=gradient(length_info(2,:),length_info(1,2)-length_info(1,1));
        dL=interp1(length_info(1,:),dL_full,t);

    else % finding stress for the whole simulation (tvals should match L vals)
        L=length_info(2,:)';
        if length(L)~=length(t) % weird fitting situation where t is cropped
            L=L(1:length(t));
        end
        dL=gradient(L,t(2)-t(1));
    end
    L0=length_info(2,1);

else % isometric/isosarcometric simulation (single SL value)
    L0=length_info;
    
    if nargin==3    % for models without internal shortening
        L=length_info;
        dL=0;
    else
        L=y8; % L is determined by the ODE
        dL='iso';
    end
end
if length(ca_info)==1
    Ca=ca_info;
elseif length(t)==1 % in solver mode (step by step)
    % The ODE solver won't take the same steps as the Ca input so 
    % we need to find closest index and the corresponding Ca
    [~,idx]=min(abs(ca_info(1,:)-t)); 
    Ca=ca_info(2,idx);

else % finding stress for the whole simulation (tvals should match C vals)
    Ca=ca_info(2,:)';
end
end 

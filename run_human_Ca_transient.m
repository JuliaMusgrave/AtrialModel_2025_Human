% Script to run Musgrave et al. 2025 human atrial muscle model with 
% non-diabetic or diabetic parameters.
% Demonstrates isometric twitches and work-loops. 

% Author: Julia Musgrave
% Date: September 2024

clear
% set diabetic or non-diabetic model
diabetic=false;

model=@Mmodel_2025_Human;

% ca handling parameters are all same for diabetic vs ND, except Ca50
params.ca=load('thin_fil_ps.mat','ca').ca;

if diabetic
    params.xb=load('D_xb_fit','x_i').x_i;
    % note that adjusted x takes into account the proportion of XBs still
    % 'non-permissible' at maximal activation
    params.passive=load('D_pass_fit','xPFL').xPFL;
    params.ca(1)=0.408 % diabetic Ca50 (uM)
    ca_T=load('Ca_transients_paper.mat','d_Ca').d_Ca
    params.M_frac=0.292;
else
    params.xb=load('ND_xb_fit','x_i').x_i;
    params.passive=load('ND_pass_fit','xPFL').xPFL;
    params.ca(1)=0.33; % non-diabetic Ca50 (uM)
    ca_T=load('Ca_transients_paper.mat','nd_Ca').nd_Ca;
    params.M_frac=0.361;
end

params.met=[5 1]; % baseline of 5 mM ATP and 1 mM Pi
params.mode='sarcomere';

L=2.2; % (initial) sarcomere length

% ca transient input set up 
t=linspace(0,1,1001);
c_input = [t;ca_T];

%% run basic isometric twitch

% Get initial conditions
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

% Run the twitch to steady state
y_last(7)=1000;
y_curr=y0;
%run the sim until steady state on N (proportion non-permissive XBs)
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end

% calculating force trace from solution
[~,F_twitch,Fp]=model(t,y,L,c_input,params);


% simple plot
figure
plot(t,F_twitch,"LineWidth",2)
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')

%% run and plot basic work loops

% run WL simulation function
outputs=WL_function(model,params,c_input);

% plotting the full workloops
figure('Units','normalized','Position',[0.3 0.4 0.25 0.39]);
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

for i=1:10

nexttile(1)
hold on
plot(outputs.time{i},outputs.stress{i},'k')
if i==10
    plot(outputs.time{i},outputs.stress{i},'k','LineWidth',2)
    xlabel('Time (s)')
    ylabel('Stress (kPa)')
end

nexttile(2)
hold on
plot(outputs.time{i},outputs.length{i}/L,'k')
if i==10
    plot(outputs.time{i},outputs.length{i}/L,'k','LineWidth',2)
    xlabel('Time (s)')
    ylabel('L/L_o')
end

nexttile(3)
hold on
plot(outputs.length{i}/L,outputs.stress{i},'k')
if i==10
    plot(outputs.length{i}/2.2,outputs.stress{i},'k','LineWidth',2)
    xlabel('L/L_o')
    ylabel('Stress (kPa)')
end

end

% plotting the WL output metrics
figure('Position',[100 100 758 430]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact')

afterloads=outputs.after;
norm_after=outputs.after/outputs.after(end); % relative afterload normalised to isometric stress

nexttile(1)
plot(norm_after,outputs.Work,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('Work (kJ\cdotm^-^3)')
xlim([0 1])
box off


nexttile(2)
plot(norm_after,outputs.shortening,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('Shortening Extent (%)')
xlim([0 1])
box off

nexttile(4)
plot(norm_after,outputs.velocity,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('Shortening Velocity (s^-^1)')
xlabel ('Relative Afterload')
xlim([0 1])
box off

nexttile(5)
plot(norm_after,outputs.velocity.*afterloads,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('Short. Power (kW\cdotm^-^3)')
xlabel ('Relative Afterload')
xlim([0 1])
box off

nexttile(3)
plot(norm_after,outputs.Energy,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('XB Energy (kJ\cdotm^-^3)')
box off
xlim([0 1])

nexttile(6)
plot(norm_after,outputs.Work./outputs.Energy*100,":.",'MarkerSize',18,'LineWidth',1.5)
ylabel('XB Efficiency (%)')
xlabel ('Relative Afterload')
box off
xlim([0 1])

%% Example of running the XB only model - sinusoidal length perturbation
% Note use of 'x_p' rather than x_i parameters

model=@XBmodel_2025_Human;

if diabetic
    x=load('D_xb_fit','x_p').x_p;
else
    x=load('ND_xb_fit','x_p').x_p;
end

mets=[5 1]; % 5 mM ATP and 1 mM Pi

L0=2.2; % sarcomere length

% length change input set up 
t=linspace(0,1,1000);
L=L0+0.1*sin(t*4*pi); % 2 Hz sinusoidal around L0
l_input = [t;L];

% Get initial conditions
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

% Run the model
prev=100;
curr=0;
%run the sim until steady state force is reached
while abs((prev-curr)/curr)>1e-5
[t,y]=ode15s(@(t,y)model(t,y,l_input,x,mets),t,y0,options);
y0=y(end,:);
prev=curr;
[~,curr]=model(t(end),y0,l_input,x,mets);

end

% calculating force output from solution
[~,F_sine]=model(t,y,l_input,x);


% simple plot
figure
plot(t,F_sine,"LineWidth",2)
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')

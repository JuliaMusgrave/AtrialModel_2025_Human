% This script runs the muscle model through a sensitivity analysis of all the
% inputs vs all the key outputs

% Author: Julia Musgrave
% Date: May 2024

function[]=full_sensitivity_for_paper(diabetic)
model=@Mmodel_2024_Human;


global isotonic;
L=2.2; 
factor=1.1;

if diabetic
    xb_ps=load('D_xb_fit','x_i').x_i;
    M_frac=0.292;
    pass_ps=load('D_pass_fit','xPFL').xPFL;
    Ca=load('Ca_transients_paper.mat','d_Ca').d_Ca;
    Ca50=0.408;
else
    xb_ps=load('ND_xb_fit','x_i').x_i;
    M_frac=0.361;
    pass_ps=load('ND_pass_fit','xPFL').xPFL;
    Ca=load('Ca_transients_paper.mat','nd_Ca').nd_Ca;
    Ca50=0.33;
end

starting_params=[xb_ps([1:5 7:9 11:13]), pass_ps, Ca50 5 100 2 1000 326, 1,1,1];

row_names={'$k_1''$','$k_{-1}''$','$k_2$','$k_{-2}''$','$k_3''$','$\phi_\mathrm{v}$',...
    '$\phi_\mathrm{l}$','$K$','$\phi_\mathrm{s1}$','$\phi_\mathrm{s3}$','$k_\mathrm{dATP}$',... xb
    '$K_\mathrm{P1}$','$K_\mathrm{P2}$','$\eta$','$\phi_\mathrm{p}$',... passive
    'Ca$_{50}$','$n_\mathrm{Tm}$','$k_\mathrm{trpn}$','$n_\mathrm{trpn}$','$k_0$','$k_{-0}$',... ca
    'Ca$_\mathrm{d}$','Ca$_\mathrm{amp}$','$t_{95}$'}; %ca transient
column_names={'$F$ max','$EC_{50}$',... SS
    'diastolic stress','systolic stress','twitch amplitude','twitch duration',... twitch
    'max work','max shortening','max velocity','max power','max energy','max efficiency'}; % WL

sensitivity=zeros(length(starting_params),length(column_names));

%constant sim params
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);
y0=model();


for p=0:length(starting_params)

ps=starting_params;
if p>0
    ps(p)=ps(p)*factor;
end

params.xb=[ps(1:5) 1 ps(6:8) 0 ps(9:11) 0];
params.passive=ps(12:15);
params.ca=ps(16:21); 
Ca_d_scale=ps(22);
Ca_amp_scale=ps(23);
Ca_dur_scale=ps(24);
params.M_frac=M_frac;


params.met=[5 1];

params.mode='sarcomere';

% ca transient (baseline)
Ca_dias=Ca(1);
Ca_act=Ca-Ca_dias;
Ca_amp=max(Ca_act);
% C_smooth=interp1(t,mean_Ca,linspace(0,1,Fs));
% [~,idx]=findpeaks(-abs(F_smooth-(F_dias+0.05*F_amp)),Fs,'SortStr','descend');
% t95=abs((idx(2)-idx(1))*1000);

% adjustments
new_amp=(Ca_amp-Ca_dias*(Ca_d_scale-1));
Ca=Ca_dias*Ca_d_scale+Ca_act*(new_amp/Ca_amp)*Ca_amp_scale;

t=linspace(0,1*Ca_dur_scale,1001);
c_input=[t; Ca];

%_______________________________________________________
% set up/pre sims

% Adjust # of XBs based on true frac available at pCa 4.5
[~,~,Fa_max] = XBmodel_2024_linear_perms(params.xb,[1 0 0 0 1],3,params.met);


% get pCa 4.5
[~,y]=SSsim_par(model,[0 1],y0,L,32,params);
[~,F,Fp]=model(0,y(end,:),L,32,params);
Fa_ca=F-Fp;

params.xb(9)=params.xb(9)/Fa_ca*Fa_max;

% get diastolic SS
y0(8)=L;
[~,y]=SSsim_par(model,[0 0.1],y0,L,Ca_dias,params);
y_dias=y(end,:);
y_dias(9)=0; %IntF

%_______________________________________________________
% static F-Ca
n=50;
Cas=logspace(-1,1,n);
Fas=zeros(1,n);
for i=1:50
    c=Cas(i);
    [~,y]=SSsim_par(model,[0 1],y0,L,c,params);
    ySS=y(end,:);
    [~,F,Fp]=model(0,ySS,L,c,params);
    Fas(i)=F-Fp;
end

% find EC50
F50=max(Fas)/2;
Ca_interp=logspace(-1, 1, 5000);
Fas=interp1(Cas,Fas,Ca_interp);
[~,i50]=min(abs(Fas-F50));
EC50=Ca_interp(i50);

outputs(1:2)=[Fa_max,EC50];

% figure(2)
% semilogx(Cas,Fas)
% 
% hold on

%_______________________________________________________
% isosarcometric twitch
y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end

%calculating force trace
[~,F_twitch,~]=model(t,y,L,c_input,params);

%finding output parameters
F_dias=min(F_twitch);
F_sys=max(F_twitch);
F_amp=F_sys-F_dias;
Fs=2000;
F_smooth=interp1(t,F_twitch,linspace(0,1,Fs));
[~,idx]=findpeaks(-abs(F_smooth-(F_dias+0.05*F_amp)),Fs,'SortStr','descend');
t95=abs((idx(2)-idx(1))*1000);

outputs(3:6)=[F_dias,F_sys,F_amp, t95];

afterloads=linspace(F_dias*1.05,F_sys,10);
%_______________________________________________________
% run basic loop
params.mode='loop';
tspan=[0 1];

for loop=1:10
    a=afterloads(loop);
    params.loop=[a 0.00001 0.00005]; %afterload, mass, visc
    y_last(5)=1000;
    y_curr=y_dias;
    %run the sim until steady state on F1
    while abs(y_last(5)-y_curr(5))>1e-5
        y_last=y_curr;
        % need a trigger for isotonic
        isotonic = false;
        [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),tspan,y_last,options);
        y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
    end
    isotonic = false;
    % in workloop mode, we have to go one point at a time to get force
    F=zeros(length(t),1);
    Fp=F; dL=F; k3=F;
    for i=1:length(t)
    
        [dy, F(i), Fp(i),k3(i)] = model(t(i),y(i,:),L,c_input,params);
        dL(i) = dy(8);
    end
    energy(loop)=abs(trapz(t,k3.*y(:,2)*0.25*M_frac))*30; % unitless?
    work(loop)=abs(trapz(y(:,8)/L,F)); % kPa.um?
    shortening(loop)=100-min(y(:,8))/L*100; % percent of Lo
    velocity(loop)=-min(movmean(dL,100))/L; % muscle length /s 
end

% finding work loop outputs
power=velocity.*afterloads;
efficiency=work./energy*100;

max_work=max(interp1(afterloads,work,linspace(F_dias*1.05,F_sys,1000)));
max_short=max(shortening);
max_vel=max(velocity);
max_power=max(interp1(afterloads,power,linspace(F_dias*1.05,F_sys,1000)));
max_energy=max(energy);
max_eff=max(interp1(afterloads,efficiency,linspace(F_dias*1.05,F_sys,1000)));

outputs(7:12)=[max_work,max_short,max_vel,max_power,max_energy,max_eff];



if p>0
    sensitivity(p,:)=(outputs-baseline)./baseline*100;
else
    baseline=outputs;
end

end
%% Beautiful plot!
%removing SL from the plotting data to make others more visible
sens_WL=sensitivity;
map_res=100;

figure('Position',[208,216,591,633])
subplot('Position',[0.1 0.1 0.7 0.7])
image(sens_WL*map_res/20+map_res+1)
[rows, cols]=size(sens_WL);
red=[166, 33, 33]/255;
blue=[58, 97, 201]/255;
map=[[linspace(blue(1),1,map_res+1),linspace(1-(1-red(1))/map_res,red(1),map_res)]',...
    [linspace(blue(2),1,map_res+1),linspace(1-(1-red(2))/map_res,red(2),map_res)]',...
    [linspace(blue(3),1,map_res+1),linspace(1-(1-red(3))/map_res,red(3),map_res)]'];
c=colorbar('Limits',[0 map_res*2],'Position',[0.855 0.825 0.023 0.15]);
c.Ticks=[0 map_res/2 map_res 3*map_res/2 map_res*2];
c.TickLabels={'-20','-10','0','10','20'};
c.Label.String = 'Change in output (%)';
colormap(map);

set(gca,'TickLabelInterpreter','latex')
yticks(1:rows);
yticklabels(row_names)
xticks(1:cols)
xticklabels(column_names)

hold on
plot([0 13],[11.5 11.5],'k') % XB demarcation
plot([0 13],[15.5 15.5],'k') % passive
plot([0 13],[21.5 21.5],'k') %Ca
plot([0 13],[24.5 24.5],'k') %Cat

plot([2.5 2.5],[0 28],'k') % SS force stuff
plot([6.5 6.5],[0 28],'k') % twitch stuff

for i=1:length(column_names)
    [~,row]=max(abs(sens_WL(:,i)));
    plot([i-0.5 i+0.5],[row-0.45 row-0.48],'k','LineWidth',2)
    plot([i-0.5 i+0.5],[row+0.45 row+0.48],'k','LineWidth',2)
    plot([i-0.48 i-0.48],[row-0.5 row+0.5],'k','LineWidth',2)
    plot([i+0.48 i+0.48],[row-0.5 row+0.5],'k','LineWidth',2)
end
hold off

annotation('textarrow',[0.018 0.058],[0.7 0.7],'String','\bf{XB parameters}','LineStyle','none','HeadStyle','none','TextRotation',90)
annotation('textarrow',[0.018 0.058],[0.45 0.45],'String','\bf{Passive}','LineStyle','none','HeadStyle','none','TextRotation',90)
annotation('textarrow',[0.018 0.058],[0.33 0.33],'String','\bf{Thin filament}','LineStyle','none','HeadStyle','none','TextRotation',90)
annotation('textarrow',[0.018 0.058],[0.17 0.17],'String','\bf{[Ca^{2+}]}','LineStyle','none','HeadStyle','none','TextRotation',90)

xlim([2.5 12.5])

subplot('Position', [0.1 0.81 0.7 0.18])
mean_col=mean(abs(sensitivity),1);
bar(1:cols,mean_col,'k')
xlim([0.5 cols+0.5])
box off
xticks({})
hold on
plot([2.5 2.5],[0 max(ylim)],'k') % SS force stuff
plot([6.5 6.5],[0 max(ylim)],'k') % twitch stuff
hold off
ylabel('Mean abs change (%)')
xlim([2.5 12.5])

subplot('Position', [0.81 0.1 0.18 0.7])
mean_row=mean(abs(sensitivity),2);
barh(1:rows,mean_row,'k')
set(gca,'YDir','reverse')
ylim([0.5 rows+0.5])
%xlim([0 30])   % needed if length is there
yticks({})
box off
hold on
plot([0 max(xlim)],[11.5 11.5],'k') % XB demarcation
plot([0 max(xlim)],[15.5 15.5],'k') % passive
plot([0 max(xlim)],[21.5 21.5],'k') %Ca
plot([0 max(xlim)],[24.5 24.5],'k') %Cat

hold off
xlabel('Mean abs change (%)')

end
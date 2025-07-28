% This script plots all of the figures for the paper: "Metabolite-sensitive 
% cross-bridge models of human atria reveal the impact of diabetes on 
% muscle mechanoenergetics" - Musgrave et al. (2025). 

% Requires the following data files:
% - ave_human_fitting_data.mat - to plot the model fits to data
% - D_xb_fit.mat  \ cross-bridge model parameters for each group
% - ND_xb_fit.mat /
% - D_pass_fit.mat  \ passive model parameters for each group
% - ND_pass_fit.mat /
% - Ca_transients_paper.mat  - Ca transients for each group
% - thin_fil_ps.mat  - thin filament parameters

% and calls following functions (directly or indirectly):
% - XBmodel_2024_linear_perms
% - Mmodel_2025_Human
% - SSsim_par
% - get_sim_inputs
% - WL_function
% - full_sensitivity_for_paper


% Author: Julia Musgrave
% Date: February 2025

clear

%% Figure2 - Passive model fit
red=[0.85,0.325,0.098];

figure('Position',[149,110,792,205])
tiledlayout(1,3,'TileSpacing','compact','Padding','tight')

%CM first
load('ave_human_fitting_data.mat','D_pass_data','ND_pass_data','freqs')
fs=logspace(-1,2,100);
Yn=passive_model_linear(load('ND_pass_fit','xPFL').xPFL);
Yd=passive_model_linear(load('D_pass_fit','xPFL').xPFL);

nexttile(1)
errorbar(freqs,real(ND_pass_data{1,2}),ND_pass_data{2,2},'.','LineStyle','none','MarkerSize',12,'Color',[0 0 0])
hold on
ax=gca;
ax.XScale='log';
errorbar(freqs,real(D_pass_data{1,2}),D_pass_data{2,2},'.','LineStyle','none','MarkerSize',12,'Color',red)
semilogx(fs,real(Yn),'k-','LineWidth',1)
semilogx(fs,real(Yd),'LineWidth',1,'Color',red)
xlim([0.1 99.9])
xticks([0.1,1,10,99.9])
xticklabels({'0.1','1','10','100'})
xlabel('Frequency (Hz)')
ylabel('Elastic Modulus (MPa)      ')
box off
text(0.02,max(ylim),'A','FontSize',16,'FontWeight','bold')

nexttile(2)
errorbar(freqs,imag(ND_pass_data{1,2}),ND_pass_data{3,2},'.','LineStyle','none','MarkerSize',12,'Color',[0 0 0])
hold on
ax=gca;
ax.XScale='log';
errorbar(freqs,imag(D_pass_data{1,2}),D_pass_data{3,2},'.','LineStyle','none','MarkerSize',12,'Color',red)
semilogx(fs,imag(Yn),'k-','LineWidth',1)
semilogx(fs,imag(Yd),'LineWidth',1,'Color',red)
xlim([0.1 99.9])
xticks([0.1,1,10,99.9])
xticklabels({'0.1','1','10','100'})
xlabel('Frequency (Hz)')
ylabel('Viscous Modulus (MPa)      ')
box off
text(0.018,max(ylim),'B','FontSize',16,'FontWeight','bold')

% FL fit
SL=[1 0.95 0.9 0.85];
SLm=linspace(min(SL),max(SL),100);
for i=1:100
    Fpn(i)=passive_model_SS(SLm(i)*2.2,load('ND_pass_fit','xPFL').xPFL);
    Fpd(i)=passive_model_SS(SLm(i)*2.2,load('D_pass_fit','xPFL').xPFL);
end
nexttile(3)
errorbar(SL,ND_pass_data{4,2},ND_pass_data{5,2},'.','LineStyle','none','MarkerSize',14,'Color',[0 0 0])
hold on
errorbar(SL,D_pass_data{4,2},D_pass_data{5,2},'.','LineStyle','none','MarkerSize',14,'Color',red)

plot(SLm,Fpn,'Color',[0 0 0],'LineWidth',2)
plot(SLm,Fpd,'Color',red,'LineWidth',2)
ylabel('Passive stress (kPa)    ')
xlabel('L/L_o')
legend('Non-diabetic','Diabetic','','','Location','northwest')
box off
text(0.827,max(ylim),'C','FontSize',16,'FontWeight','bold')

%% Figure 3 - plot Ca inputs

figure('Position',[100, 100,411,285])

% get Ca data
ca_nd=load('Ca_transients_paper.mat','nd_Ca').nd_Ca;
ca_d=load('Ca_transients_paper.mat','d_Ca').d_Ca;
t=linspace(0,1,1001);

plot(t,ca_nd,'k',"LineWidth",2)
hold on
plot(t,ca_d,"LineWidth",2)
box off
ylabel('[Ca^2^+] (\muM)')
xlabel('Time (s)')
legend('Non-diabetic','Diabetic')
ylim([0 0.6])

%% Figure 4: plotting fits for both XB models
load('ave_human_fitting_data.mat', 'D_xb_data','ND_xb_data','freqs')

freqs=freqs(1:11); % didn't fit to 100 Hz
sd=[1 0 0 0 1]; % strain dependencies in human model
md=3; % metabolite dependence in human model
c=colororder;
[blue, red, green, purple]=deal(c(1,:),c(2,:),c(5,:),c(4,:));
    
figure('Position',[104,137,872,642])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
colour={[0 0 0],blue,red,green,purple};

for d=[0 1]
    if d
        data=D_xb_data;
        ps=load('D_xb_fit','x_p').x_p;
    else
        data=ND_xb_data;
        ps=load('ND_xb_fit','x_p').x_p;
    end
for i=1:width(data)-1 % plotting each metabolite
        [~,Yam]=XBmodel_2024_linear_perms(ps,sd,md,data{2,i+1});
    
        nexttile(1+d*2)
        errorbar(freqs,real(data{3,i+1}(1:11)),data{4,i+1}(1:11),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
        hold on
        ax=gca;
        ax.XScale='log';
        ax.FontSize=12;
        semilogx(logspace(-1,2,100),real(Yam),'LineWidth',1,'Color',colour{i})
        xlim([0.1 100])
        xticks([0.1 1 10 100])
        xticklabels({'0.1' '1' '10' '100'})
        if d; xlabel('Frequency (Hz)'); end
        ylabel('Elastic Modulus (MPa)')
        box off
    
        nexttile(2+d*2)
        errorbar(freqs,imag(data{3,i+1}(1:11)),(data{5,i+1}(1:11)),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
        hold on
        ax=gca;
        ax.XScale='log';
        ax.FontSize=12;
        semilogx(logspace(-1,2,100),imag(Yam),'LineWidth',1,'Color',colour{i})
        xlim([0.1 100])
        xticks([0.1 1 10 100])
        xticklabels({'0.1' '1' '10' '100'})
        if d; xlabel('Frequency (Hz)'); end
        ylabel({'','Viscous Modulus (MPa)'})
        box off
    
end
end
    
    % tidying figure
    nexttile(1)
    legend('','(B) Baseline','','(0.1A) 0.1 mM ATP','','(1A) 1 mM ATP','','(0P) 0 mM Pi','','(10P) 10 mM Pi','Location','northwest','FontSize',11)
    text(0.03,max(ylim),'A','FontSize',24,'FontWeight','bold')
    text(111,0.3,{'0.1A','1A','0P','B','10P'},'FontSize',12)
    nexttile(3)
    text(0.027,max(ylim),'B','FontSize',24,'FontWeight','bold')
    text(111,0.25,{'0.1A','1A','0P','B','10P'},'FontSize',12)


%% Figure 5: Isometric twitch - trabecula/muscle vs myocyte model

red=[0.85,0.325,0.098];

model=@Mmodel_2025_Human;
clear F_twitch

figure('Position',[456,441,614,257])
subplot('Position',[0.07 0.15 0.43 0.78])
hold on
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')
ylim([0 50])
text(-0.13,50,'A','FontSize',14)
title('Muscle')
subplot('Position',[0.56 0.15 0.43 0.78])
hold on
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')
text(-0.13,50,'B','FontSize',14)
title('Myocyte')

for diabetic=[0 1]

[params,ca]=load_inputs(diabetic,diabetic);

if diabetic
    col=red;
else
    col = [0 0 0];
end


params.met=[5 1];
params.mode='sarcomere';
L=2.2;

% ca transient input
dias=ca(1);
t=linspace(0,1,1001);
c_input = [t;ca];

% Run to steady state
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);

y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end

%calculating force trace
[~,F_twitch(:,diabetic+1),Fp]=model(t,y,L,c_input,params);


% ---MYO SIMULATION----
params.xb(9)=params.xb(9)/params.M_frac;

% Run to steady state
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);

y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end

%calculating force trace
[~,F,Fp]=model(t,y,L,c_input,params);
F_twitchmyo(:,diabetic+1)=F;

subplot('Position',[0.07 0.15 0.43 0.78])
plot(t,F_twitch(:,diabetic+1),"LineWidth",2,"Color",col)
subplot('Position',[0.56 0.15 0.43 0.78])
plot(t,F_twitchmyo(:,diabetic+1),"LineWidth",2,"Color",col)

end



legend('Non-diabetic','Diabetic')

% some quick analysis of the twitches
[d_amp,d_t]=twitch_analysis(F_twitch(:,2),t);
[n_amp,n_t]=twitch_analysis(F_twitch(:,1),t);
d_sys=max(F_twitch(:,2));
n_sys=max(F_twitch(:,1));
% d_amp/n_amp
% d_t/n_t
% d_sys/n_sys

%% Figure 6 & 7: Work loops - baseline w complete models

red=[0.85,0.325,0.098];

model=@Mmodel_2025_Human;
WLfig_gen=figure('Units','normalized','Position',[0.3 0.4 0.31 0.39]);
tiledlayout(3,2,'TileSpacing','tight','Padding','tight')

WLfig_out=figure('Position',[100 100 758 430]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact')

for diabetic=[0 1]

[params,ca]=load_inputs(diabetic,diabetic);

if diabetic; col=red;
else; col=[0 0 0];  
end


params.met=[5 1];

% ca transient 
dias=ca(1);
t=linspace(0,1,1001);
c_input = [t;ca];

% run simulation
outputs=WL_function(model,params,c_input);

% plotting the full workloops
figure(WLfig_gen)

set(WLfig_gen,'defaultAxesColorOrder',[[0 0 0];[0 0 0]]);
for i=1:10

ax=nexttile(diabetic+1);
ax.FontSize=10;
hold on
plot(outputs.time{i},outputs.stress{i},'Color',col)
if i==10
    plot(outputs.time{i},outputs.stress{i},'Color',col,'LineWidth',2)
    ylim([0 21])
    xlabel('Time (s)')
    yyaxis right
    ylim([0 21/max(outputs.stress{10})])
    yticks([0 0.5 1])
    ylabel('Relative Afterload')
    yyaxis left
end

ax=nexttile(diabetic+3);
ax.FontSize=10;
hold on
plot(outputs.time{i},outputs.length{i}/2.2,'Color',col)
if i==10
    plot(outputs.time{i},outputs.length{i}/2.2,'Color',col,'LineWidth',2)
    ylim([0.85 1.01])
    xlabel('Time (s)')
end

ax=nexttile(diabetic+5);
ax.FontSize=10;
hold on
plot(outputs.length{i}/2.2,outputs.stress{i},'Color',col)
if i==10
    plot(outputs.length{i}/2.2,outputs.stress{i},'Color',col,'LineWidth',2)
    ylim([0 21])
    xlabel('L/L_o')
    xlim([0.85 1.01])
    yyaxis right
    ylim([0 21/max(outputs.stress{10})])
    yticks([0 0.5 1])
    ylabel('Relative Afterload')
    yyaxis left
end

end
nexttile(1); title('Non-diabetic','FontSize',12,'FontWeight','bold')
ylabel('Stress (kPa)')
nexttile(2); title('Diabetic','FontSize',12,'FontWeight','bold')
ylabel('Stress (kPa)')
nexttile(3); ylabel('L/L_o')
nexttile(4); ylabel('L/L_o')
nexttile(5); ylabel('Stress (kPa)')
nexttile(6); ylabel('Stress (kPa)')


% plotting the WL output metrics
six_panel_WL(WLfig_out,outputs,':.',18,col,col)
if diabetic
    nexttile(2); 
    legend('Non-diabetic','Diabetic','Location','southwest','FontSize',10); 

end
end
nexttile(1)
text(-0.28,max(ylim),'A','FontSize',18,'FontWeight','bold')
nexttile(2)
text(-0.25,max(ylim),'B','FontSize',18,'FontWeight','bold')
nexttile(3)
text(-0.25,max(ylim),'C','FontSize',18,'FontWeight','bold')
nexttile(4)
text(-0.25,max(ylim),'D','FontSize',18,'FontWeight','bold')
nexttile(5)
text(-0.25,max(ylim),'E','FontSize',18,'FontWeight','bold')
nexttile(6)
text(-0.25,max(ylim),'F','FontSize',18,'FontWeight','bold')




%% Figure 8: Raised Pi
clear F_twitch max_eff max_energy max_pow max_short max_vel max_work
red=[0.85,0.325,0.098];

model=@Mmodel_2025_Human;

pis=linspace(0,10,11);
pis(1)=1e-6;
for i=1:11
p=pis(i);
for diabetic=[0 1]

[params, ca]=load_inputs(diabetic,diabetic);

params.mode='sarcomere';
params.met=[5 p];

% ca transient 
dias=min(ca);

t=linspace(0,1,1001);
c_input=[t; ca];
outputs=WL_function(model,params,c_input);
afters=outputs.after;
workint=interp1(afters,outputs.Work,linspace(min(afters),max(afters),1000));
powint=interp1(afters,outputs.velocity.*afters,linspace(min(afters),max(afters),1000));
max_work(i,diabetic+1)=max(workint);
max_short(i,diabetic+1)=max(outputs.shortening);
max_vel(i,diabetic+1)=max(outputs.velocity);
max_pow(i,diabetic+1)=max(powint);
max_energy(i,diabetic+1)=max(outputs.Energy);
max_eff(i,diabetic+1)=max(outputs.Work./outputs.Energy*100);

% Run isometric for two options
if p==1 || p==10
L=2.2;
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);

y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end
%calculating force trace
[~,F,Fp]=model(t,y,L,c_input,params);
F_twitch(:,diabetic+1+round(p/5))=F;
end

end
end
% plot setup
figure('Position',[445,230,790,278])
tiledlayout(2,5,'TileSpacing','tight','Padding','compact')
% isometric plot
nexttile(1,[2 2])
plot(t,F_twitch(:,1),'k',"LineWidth",2)
hold on
plot(t,F_twitch(:,2),"LineWidth",2,'Color',red)
plot(t,F_twitch(:,3),'k--',"LineWidth",2)
plot(t,F_twitch(:,4),'--',"LineWidth",2,'Color',red)
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')
ylim([0 25])
l=legend('Non-diabetic (1 mM P_i)','Diabetic (1 mM P_i)'...
    ,'Non-diabetic (10 mM P_i)','Diabetic (10 mM P_i)');
l.ItemTokenSize=[20 18];
text(-0.12,max(ylim),'A','FontSize',14,'FontWeight','bold')

nexttile(3)
plot(pis,max_work(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_work(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 0.65])
ylabel('Max Work (kJ\cdotm^-^3)      ')
text(-3,max(ylim),'B','FontSize',12,'FontWeight','bold')

nexttile(4)
plot(pis,max_short(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_short(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 15])
ylabel('Max Shortening (%)        ')
text(-3,max(ylim),'C','FontSize',12,'FontWeight','bold')

nexttile(8)
plot(pis,max_vel(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_vel(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 4])
ylabel('Max Velocity (s^-^1)      ')
xlabel('[P_i] (mM)')
text(-3,max(ylim),'E','FontSize',12,'FontWeight','bold')

nexttile(9)
plot(pis,max_pow(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_pow(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 30])
ylabel('Max Power (kW\cdotm^-^3)         ')
xlabel('[P_i] (mM)')
text(-3,max(ylim),'F','FontSize',12,'FontWeight','bold')

nexttile(5)
plot(pis,max_energy(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_energy(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 41])
ylabel('Max Energy (kJ\cdotm^-^3)      ')
text(-3,max(ylim),'D','FontSize',12,'FontWeight','bold')

nexttile(10)
plot(pis,max_eff(:,1),'k.','MarkerSize',15)
hold on
plot(pis,max_eff(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 13])
ylabel('Max Efficiency (%)       ')
text(-3,max(ylim),'G','FontSize',12,'FontWeight','bold')
xlabel('[P_i] (mM)')

% some twitch analysis
[d1_amp,d1_t]=twitch_analysis(F_twitch(:,2),t);
[n1_amp,n1_t]=twitch_analysis(F_twitch(:,1),t);
[d10_amp,d10_t]=twitch_analysis(F_twitch(:,4),t);
[n10_amp,n10_t]=twitch_analysis(F_twitch(:,3),t);
% d10_amp/d1_amp
% n10_amp/n1_amp

%% Figure 9: Lowered ATP

clear F_twitch max_eff max_energy max_pow max_short max_vel max_work
red=[0.85,0.325,0.098];

model=@Mmodel_2025_Human;
atps=linspace(1,10,10);
for i=1:10
atp=atps(i);
for diabetic=[0 1]

[params, ca]=load_inputs(diabetic,diabetic);

params.mode='sarcomere';
params.met=[atp 1];

% ca transient 
dias=min(ca);

t=linspace(0,1,1001);
c_input=[t; ca];
outputs=WL_function(model,params,c_input);
afters=outputs.after;
workint=interp1(afters,outputs.Work,linspace(min(afters),max(afters),1000));
powint=interp1(afters,outputs.velocity.*afters,linspace(min(afters),max(afters),1000));
max_work(i,diabetic+1)=max(workint);
max_short(i,diabetic+1)=max(outputs.shortening);
max_vel(i,diabetic+1)=max(outputs.velocity);
max_pow(i,diabetic+1)=max(powint);
max_energy(i,diabetic+1)=max(outputs.Energy);
max_eff(i,diabetic+1)=max(outputs.Work./outputs.Energy*100);

% Run isometric for two options
if atp==5 || atp==1
L=2.2;
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);

y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end
%calculating force trace
[~,F,Fp]=model(t,y,L,c_input,params);
F_twitch(:,diabetic+1+2*mod(round(atp),5))=F;
end

end
end
% plot setup
figure('Position',[145,230,790,278])
tiledlayout(2,5,'TileSpacing','tight','Padding','compact')
% isometric plot
nexttile(1,[2 2])
plot(t,F_twitch(:,1),'k',"LineWidth",2)
hold on
plot(t,F_twitch(:,2),"LineWidth",2,'Color',red)
plot(t,F_twitch(:,3),'k--',"LineWidth",2)
plot(t,F_twitch(:,4),'--',"LineWidth",2,'Color',red)
box off
ylabel('Stress (kPa)')
xlabel('Time (s)')
ylim([0 25])
l=legend('Non-diabetic (5 mM ATP)','Diabetic (5 mM ATP)'...
    ,'Non-diabetic (1 mM ATP)','Diabetic (1 mM ATP)');
l.ItemTokenSize=[20 18];
text(-0.12,max(ylim),'A','FontSize',14,'FontWeight','bold')

nexttile(3)
plot(atps,max_work(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_work(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 0.7])
ylabel('Max Work (kJ\cdotm^-^3)    ')
text(-3,max(ylim),'B','FontSize',12,'FontWeight','bold')

nexttile(4)
plot(atps,max_short(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_short(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 15])
ylabel('Max Shortening (%)      ')
text(-3,max(ylim),'C','FontSize',12,'FontWeight','bold')

nexttile(8)
plot(atps,max_vel(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_vel(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 4])
ylabel('Max Velocity (s^-^1)  ')
xlabel('[ATP] (mM)')
text(-3,max(ylim),'E','FontSize',12,'FontWeight','bold')

nexttile(9)
plot(atps,max_pow(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_pow(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 30])
ylabel('Max Power (kW\cdotm^-^3)      ')
xlabel('[ATP] (mM)')
text(-3,max(ylim),'F','FontSize',12,'FontWeight','bold')

nexttile(5)
plot(atps,max_energy(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_energy(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 27])
ylabel('Max Energy (kJ\cdotm^-^3)       ')
text(-3,max(ylim),'D','FontSize',12,'FontWeight','bold')

nexttile(10)
plot(atps,max_eff(:,1),'k.','MarkerSize',15)
hold on
plot(atps,max_eff(:,2),'.','Color',red,'MarkerSize',15)
ylim([0 35])
ylabel('Max Efficiency (%)    ')
xlabel('[ATP] (mM)')
text(-3,max(ylim),'G','FontSize',12,'FontWeight','bold')

% twitch analysis
[d5_amp,d5_t]=twitch_analysis(F_twitch(:,2),t);
[n5_amp,n5_t]=twitch_analysis(F_twitch(:,1),t);
[d1_amp,d1_t]=twitch_analysis(F_twitch(:,4),t);
[n1_amp,n1_t]=twitch_analysis(F_twitch(:,3),t);
% d1_amp/d5_amp
% n1_amp/n5_amp
% d1_t/d5_t
% n1_t/n5_t

%% Figure 10: XB/Ca model matrix
bf=colororder;
[blue, red, green, purple]=deal(bf(1,:),bf(2,:),bf(5,:),bf(4,:));

model=@Mmodel_2025_Human;
clear F_twitch

figure('Position',[445,230,790,278])
tiledlayout(2,5,'TileSpacing','tight','Padding','compact')


for ca_stuff=[0 1]
for xb=[0 1]

[params,ca]=load_inputs(xb,ca_stuff);
if xb
    % set non-diabetic passive for diabetic model
    params.passive=load('ND_pass_fit','xPFL').xPFL;
    if ca_stuff
        marker=':square';
        m_size=5;
        col=red;
        dash='-';
    else
        marker=':o';
        m_size=4;
        col=purple;
        dash=':';
    end
    fill='w';
else
    if ca_stuff
        marker=':.';
        m_size=15;
        col=green;
        fill=green;
        dash='--';
    else
        marker=':square';
        m_size=4;
        col=[0 0 0];
        fill=[0 0 0];
        dash='-';
    end  
end


params.met=[5 1];
params.mode='sarcomere';
L=2.2;

% ca transient input
dias=ca(1);
t=linspace(0,1,1001);
c_input = [t;ca];

% Run to steady state
y0=model();
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

tspan=[0 0.1];
y0(8)=L;
[~,y]=SSsim_par(model,tspan,y0,L,dias,params);
y_dias=y(end,:);

y_last(7)=1000;
y_curr=y_dias;
%run the sim until steady state on N
while abs(y_last(7)-y_curr(7))>1e-5
    y_last=y_curr;
    [t,y]=ode15s(@(t,y)model(t,y,L,c_input,params),t,y_last,options);
    y_curr(1:8)=y(end,1:8); % IntF must be forced to zero at start each time
end

%calculating force trace
[~,F_twitch(:,xb+1+2*ca_stuff),Fp]=model(t,y,L,c_input,params);

% Work loop simulations
outputs=WL_function(model,params,c_input); 

% isometric plot
nexttile(1,[2 2])
plot(t,F_twitch(:,xb+1+2*ca_stuff),dash,"LineWidth",2,"Color",col)
hold on
box off
ylabel('Total Stress (kPa)')
xlabel('Time (s)')
ylim([0 25])
text(-0.12,max(ylim),'A','FontSize',14,'FontWeight','bold')

afterloads=outputs.after;
norm_after=outputs.after/outputs.after(end);

nexttile(3)
plot(norm_after,outputs.Work,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])

nexttile(4)
plot(norm_after,outputs.shortening,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])
ylabel('Shortening Extent(%)')

nexttile(8)
plot(norm_after,outputs.velocity,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])
ylabel('Max shortening velocity (s^-^1)')
xlabel('Relative Afterload')

nexttile(9)
plot(norm_after,afterloads.*outputs.velocity,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])
ylabel('Shortening Power (kW\cdotm^-^3)')
xlabel('Relative Afterload')

nexttile(5)
plot(norm_after,outputs.Energy,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])
ylabel('XB Energy (kJ\cdotm^-^3)')

nexttile(10)
plot(norm_after,outputs.Work./outputs.Energy*100,marker,'MarkerSize',m_size,"MarkerFaceColor",fill,"Color",col,"LineWidth",1.5)
hold on
xlim([0 1])
ylabel('XB Efficiency (%)')
xlabel('Relative Afterload')

end
end

nexttile(1,[2 2])
legend('ND XBs, ND Ca','D XBs, ND Ca','ND XBs, D Ca','D XBs, D Ca')
text(-0.25,max(ylim),'A','FontSize',14,'FontWeight','bold')

nexttile(3);
ylabel('Work (kJ\cdotm^-^3)')
text(-0.29,max(ylim),'B','FontSize',12,'FontWeight','bold')

nexttile(4); 
    l=legend('  ND/ND','  D/ND','  ND/D','  D/D','Position',[0.629,0.58,0.11,0.11],...
        'NumColumns',2,'Orientation','horizontal'); 
    legend('boxoff')
    l.ItemTokenSize=[0 18];
    annotation(gcf,'textbox',...
    [0.6504 0.6906 0.0648 0.0755],'String',{'XB/Ca'},...
    'FitBoxToText','off','EdgeColor','none');
    annotation(gcf,'line',[0.6544 0.7101],[0.6896 0.6871]);

ylabel('Shortening Extent(%)       ')
text(-0.25,max(ylim),'C','FontSize',12,'FontWeight','bold')

nexttile(8)
ylabel('Shortening Velocity (s^-^1)            ')
text(-0.25,max(ylim),'E','FontSize',12,'FontWeight','bold')

nexttile(9)
ylabel('Short. Power (kW\cdotm^-^3)          ')
text(-0.25,max(ylim),'F','FontSize',12,'FontWeight','bold')

nexttile(5)
text(-0.25,max(ylim),'D','FontSize',12,'FontWeight','bold')
ylabel('XB Energy (kJ\cdotm^-^3)    ')

nexttile(10)
ylabel('XB Efficiency (%)    ')
text(-0.25,max(ylim),'G','FontSize',12,'FontWeight','bold')

% twitch analysis
[d_amp,d_t]=twitch_analysis(F_twitch(:,4),t);
[n_amp,n_t]=twitch_analysis(F_twitch(:,1),t);
[x_amp,x_t]=twitch_analysis(F_twitch(:,2),t);
[c_amp,c_t]=twitch_analysis(F_twitch(:,3),t);
% d_amp/n_amp
% x_amp/n_amp
% c_amp/n_amp
% d_t/n_t
% x_t/n_t
% c_t/n_t

%% Figure 11: muscle sensitivity analysis

full_sensitivity_for_paper(false)

%% helper functions

% plotting WL outputs
function []=six_panel_WL(work_pl,outputs,marker,m_size,col,fill)
    % Work vs afterload
    afterloads=outputs.after;
    norm_after=outputs.after/outputs.after(end);
    
    figure(work_pl)
    nexttile(1)
    hold on
    plot(norm_after,outputs.Work,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('Work (kJ\cdotm^-^3)')
    xlim([0 1])
    box off
    set(gca,'FontSize',11)

    
    nexttile(2)
    hold on
    plot(norm_after,outputs.shortening,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('Shortening Extent (%)         ')
    xlim([0 1])
    box off
    set(gca,'FontSize',11)
    
    nexttile(4)
    hold on
    plot(norm_after,outputs.velocity,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('Shortening Velocity (s^-^1)           ')
    xlabel ('Relative Afterload')
    xlim([0 1])
    box off
    set(gca,'FontSize',11)
    
    nexttile(5)
    hold on
    plot(norm_after,outputs.velocity.*afterloads,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('Short. Power (kW\cdotm^-^3)        ')
    xlabel ('Relative Afterload')
    xlim([0 1])
    box off
    set(gca,'FontSize',11)
    
    nexttile(3)
    hold on
    plot(norm_after,outputs.Energy,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('XB Energy (kJ\cdotm^-^3)    ')
    box off
    xlim([0 1])
    set(gca,'FontSize',11)
    
    nexttile(6)
    hold on
    plot(norm_after,outputs.Work./outputs.Energy*100,marker,'Color',col,'MarkerSize',m_size,'MarkerFaceColor',fill,'LineWidth',1.5)
    ylabel('XB Efficiency (%)  ')
    xlabel ('Relative Afterload')
    box off
    xlim([0 1])
    set(gca,'FontSize',11)
end
 

%% function to set up params and Ca transient for the model 
% takes boolean input of diabetic XBs or Ca
function [params,ca_T]=load_inputs(XB, Ca)

if ~XB % ND XB model
    params.xb=load('ND_xb_fit','x_i').x_i; % intact model parameters
    params.passive=load('ND_pass_fit','xPFL').xPFL;
    params.M_frac=0.361; % directly from Musgrave et al.
else % D XB model
    params.passive=load('D_pass_fit','xPFL').xPFL;
    params.xb=load('D_xb_fit','x_i').x_i;
    params.M_frac=0.292; % directly from Musgrave et al.
end

params.ca=load('thin_fil_ps.mat','ca').ca; % Land Human parameters
if ~Ca % ND Ca model/input
    ca_T=load('Ca_transients_paper.mat','nd_Ca').nd_Ca;
    params.ca(1)=0.33; % from Jones et al.
else % D Ca model/input
    ca_T=load('Ca_transients_paper.mat','d_Ca').d_Ca;
    params.ca(1)=0.408; % from Jones et al.
end

end

%% passive model functions:

% basic model for steady state:

function [F] = passive_model_SS(L,params)


%----------------------------
% assigning the parameters 
%kPE1=params(1); % kPa
kPE2=params(2); % kPa
%eta=params(3); % s^-1
phi_e=params(4); % unitless exponential length dependence


Lr=2.2*.85;


F2= kPE2*(exp(phi_e*(L-Lr))-1);
F=F2;

end
% 
% function for the passive CM model:
function [Y] = passive_model_linear(params)

freqs = logspace(-1,2,100);

% model parameters
kPE1=params(1);
kPE2=params(2);
eta=params(3);
phi_e=params(4);

% frequency response parameters
om=freqs*2*pi;
omi=om*1i; %omega*i
L0= 2.2; % experiment SL (um)
Lr= L0*0.85; % rest SL (um)

% partial derivatives
d11=-kPE1/eta;
du1=kPE1;

HF1=du1*omi./(omi-d11);
dF2=kPE2*phi_e*exp(phi_e*(L0-Lr));

% passive complex modulus
Y=L0/1000*(HF1+dF2);
end
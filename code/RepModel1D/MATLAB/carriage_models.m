%% INITIAL
Qfs     = 10/60; % total fresh air flow rate saloon (single HVAC) [m^3/s]
Qcs     = 40/60; % total circulation rate saloon (single HVAC) [m^3/s]
phi     = 1/10; % fraction of Qd/Qcs [-]
theta   = (Qcs-Qfs)/(Qcs*(1-phi)); % fraction of Qcs/Qc [-]
zeta    = 0; % fraction of Ql/Qe [-]
qD      = 0; % downward flow: open->passenger [m^2/s]
Vs      = 100; % volume saloon [m^3]
Vv      = 10; % volume vestibule [m^3]
xd      = 20/2; % distance till vestibule door from centre [m]
x1      = 6; % distance till start HVAC from centre [m]
x2      = 8; % distance till end HVAC from centre [m]
A       = Vs/(xd*2); % cross-sectional area saloon [m^2]
H       = 2; % height saloon [m]
W       = A/H; % width saloon [m] 
p       = 0.49/3600; % breathing rate person [m^3/s]

xs      = [2 5]; % source location(s) from centre [m] 
ms      = 2/60; % source rate saloon [unit/s]
svl     = 1/60; % source rate vestibule left [unit/s]
svr     = 2/60; % source rate vestibule right [unit/s]
sigma   = 0.5; % sigma width of source profile (same for all sources) [m]

x_Nstep_half     = 200; % number of grid points for half of the carriage
x_step_st_coarse = 1; % # of times coarse of storage grid
[x,x_Nstep,x_step,X_L,X_R,X_HVAC_L,X_HVAC_R,x_st,x_idx_st] = x_grid(x1,x2,xd,x_Nstep_half,x_step_st_coarse); % make grid along x-axis

Diff    = 0.05; % turbulent diffusion coefficient
tracer  = "Nebu"; % "covid" or "CO2" or "Nebu"
[kd]    = decay_rate(tracer); % decay rate [1/s]
F       = 0; % ventilation efficiency (0 = nothing removed, 1 = all removed)

t_start     = 0; % time start (s)
t_step      = 0.01; % time step (s)
t_step_st   = 2*60; % time step of storage value (need to be a multiplication of t_step)
t_end       = 60*20; % time stop (s)
[t,t_Nstep,t_idx_st,t_idx_st_f,t_st] = t_grid(t_start,t_step,t_end,t_step_st); % make time series

cs0     = zeros(1,x_Nstep); % initial concentraiton saloon
cvl0    = 0; % intial concentraiton left vestibule
cvr0    = 0; % intial concentraiton right vestibule
ss      = source_1DM(x,x_step,x_Nstep,xs,sigma,ms); % source term saloon

[Qc,Qd,Qf,Qe,Qcv,Qr,Ql,Qx,DQx]       = flow_rates(Qcs,theta,phi,zeta,x,x1,x2,xd); % get flowrates

%% 1D-model
% allocate
[cs_ghost,cs_ghost_p1,cs_ghost_m1]  = deal(NaN(1,x_Nstep+2)); % ghost grid points for concentraiton 
cs_st                               = NaN(length(t_st),length(x_st)); % store saloon concentrations
[cvl_st,cvr_st]                     = deal(NaN(1,length(t_st))); % store vestibule concentrations

% initialize 
i   = 1; % iteration number
cs  = cs0; % saloon concentration
cvl = cvl0; % south vestibule concentration
cvr = cvr0; % north vestibule concentraiton

% store values
if ismember(i,x_idx_st)
    k           = t_idx_st_f(i);
    
    cs_st(k,:)  = cs(x_idx_st);
    cvl_st(k)   = cvl;
    cvr_st(k)   = cvr;
end

plot_ini_1D(x,xd,Qx,cs0,cvl0,cvr0,ss); % plot overview initial model

% show stored concentrations
figure(2); set(gcf, 'Position',[600 200 900 500])
plot(x,cs,'k')
hold on
plot(-xd,cvl,'bs','MarkerFaceColor','b')
plot(+xd,cvr,'bs','MarkerFaceColor','b')


i = i+1;

% run 1D model
while i<=t_Nstep  
    [cs,cvl,cvr] = iteration_1DM(cs,cvl,cvr,cs_ghost,cs_ghost_p1,cs_ghost_m1,ss,svl,svr,x1,x2,xd,X_L,X_R,X_HVAC_L,X_HVAC_R,x_step,t_step,A,Vv,F,kd,Diff,theta,Qe,Qcv,Qr,Qd,Qx,DQx); 

    % store values
    if ismember(i,t_idx_st)
        k = t_idx_st_f(i);
        
        cs_st(k,:)  = cs(x_idx_st);
        cvl_st(k)   = cvl;
        cvr_st(k)   = cvr;

        figure(2)
        plot(x,cs,'k')
        hold on
        plot(-xd,cvl,'bs','MarkerFaceColor','b')
        plot(+xd,cvr,'bs','MarkerFaceColor','b')
        hold off
        drawnow
    end
    
    i = i+1; % next iteration
end
plot_style("$x$ (m)","$c$ (unit m\textsuperscript{-3})",[]);

cs_st_1D  = cs_st;
cvl_st_1D = cvl_st;
cvr_st_1D = cvr_st;
cs_1D     = cs;

%% FMM (fully-mixed model)
ss_FMM  = length(xs)*ms;
cs0_FMM = sum(cs0)*x_Nstep;
cv0_FMM = (cvl0+cvr0)/2;
sv_FMM  = (svl+svr)/2;

% run FMM model
[cs_FFM,cv_FFM,cs1_FFM,cv1_FFM] = FMM(Qcs,Qc,Qd,Qf,Qe,Qcv,Qr,Ql,Vs,Vv,kd,F,cs0_FMM,cv0_FMM,ss_FMM,sv_FMM,t,t_end);

figure(1)
hold on
plot(t/60,cs_FFM,'r')
plot(t/60,cv_FFM,'k')
hold off
plot_style("$t$ (min)","$c$ (unit m\textsuperscript{-3})",[])

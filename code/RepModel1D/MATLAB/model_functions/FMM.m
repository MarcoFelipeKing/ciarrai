function [cs,cv,cs1,cv1,cs_int] = FMM(Qcs,Qc,Qd,Qf,Qe,Qcv,Qr,Ql,Vs,Vv,kd,F,cs0,cv0,ss,sv,t,t_end);

% essential terms
Qsv      = 2*(Qd+Qr*Qcv/Qc); % total flow from saloon into vestibule
Qss      = 2*Qr*Qcs/Qc; % saloon recirculation
QHsv     = Qr*Qcv/Qc; % saloon to vestibule through HVAC
Ss       = ss/Vs; % source term saloon
Sv       = sv/Vv; % source term vestibule
Ks       = (Qsv+F*Qss)/Vs + kd; % exponential decay saloon
Kv       = Qe/Vv + kd; % exponential decay vestibule
Gsv      = (Qd+(1-F)*QHsv)/Vv; % combination term saloon and vestibule
Sv_ef    = Sv + Gsv/Ks*Ss; % effective source term vestibule
Gsv_ef   = Gsv*(cs0-Ss/Ks); % effective combination term vestibule
cs_tilde = Ss/Ks; % steady-state concentration saloon
cv_tilde = Sv_ef/Kv; % steady-state concentration vestibule
t_Nstep  = length(t); % number of steps

on_vent  = Qcv ~= 0; % check if ventilation is turned on
on_dec   = kd  ~= 0; % check if natural decay is included


% allocate
cs      = NaN(1,t_Nstep); % concentration saloon
cv      = NaN(1,t_Nstep); % concentraiton vestibule
cs_int  = NaN(1,t_Nstep); 

cs1 = NaN; % end concentration saloon
cv1 = NaN; % end concentration vestibule

% calculate concentrations per scenario
if on_vent == true % if ventilation is turned on
    cs      = cs0*exp(-Ks*t)+Ss/Ks*(1-exp(-Ks*t));
    cv      = Sv_ef/Kv-Gsv_ef/(Ks-Kv)*exp(-Ks*t)+(cv0-Sv_ef/Kv+Gsv_ef/(Ks-Kv))*exp(-Kv*t);
    cs_int  = (cs0*Ks-Ss+exp(-Ks*t)*(-cs0*Ks+Ss)+Ks*Ss*t)/Ks^2; % concentration saloon integral over time

            
    cs1 = cs0*exp(-Ks*t_end)+Ss/Ks*(1-exp(-Ks*t_end));
    cv1 = Sv_ef/Kv-Gsv_ef/(Ks-Kv)*exp(-Ks*t_end)+(cv0-Sv_ef/Kv+Gsv_ef/(Ks-Kv))*exp(-Kv*t_end);
elseif on_vent == false && on_dec == true % if ventilation turned off, but with natural decay
    cs      = cs0*exp(-kd*t)+Ss/kd*(1-exp(-kd*t)); 
    cv      = cv0*exp(-kd*t)+Sv/kd*(1-exp(-kd*t));
    cs_int  = (cs0*kd-Ss+exp(-kd*t)*(-cs0*kd+Ss)+kd*Ss*t)/kd^2;
    
    cs1 = cs0*exp(-kd*t_end)+Ss/kd*(1-exp(-kd*t_end)); 
    cv1 = cv0*exp(-kd*t_end)+Sv/kd*(1-exp(-kd*t_end));
elseif on_vent == false && on_dec == false % if ventilation turned off and no natural decay
    cs      = Ss*t+cs0;
    cv      = Sv*t+cv0;
    cs_int  = Ss/2*t.^2+cs0*t;
    
    cs1 = Ss*t_end+cs0;
    cv1 = Sv*t_end+cv0;
end
end
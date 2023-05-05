function [Qc,Qd,Qf,Qe,Qcv,Qr,Ql,Qx,DQx,ACHs,ACHv] = flow_rates(Qcs,theta,phi,zeta,x,x1,x2,xd,Vs,Vv)
% flow rates (m^3/s)
Qc   = (1/theta)*Qcs; % HVAC outlet saloon and vestibule
Qd   = (phi)*Qcs; % from saloon into vestibule through door
Qf   = (1/theta-1+phi)*Qcs; % Fresh air
Qe   = ((1/theta-1+phi)/(1-zeta))*Qcs; % extraction vestibule
Qcv  = (1/theta-1)*Qcs; % HVAC outlet vestibule
Qr   = (1-phi)*Qcs; % HVAC inlet saloon
Ql   = (zeta*(1/theta-1+phi)/(1-zeta))*Qcs; % Leakage vestibule

if exist('Vs','var') && exist('Vv','var')
    % Air changes per hour
    ACHs = 3600/((Vs/2)/(Qcs/Qc*Qf)); % saloon
    ACHv = 3600/((Vv/2)/(Qcv/Qc*Qf)); % vestibule
end

if exist('x','var') && exist('x1','var') && exist('x2','var') && exist('xd','var')
    % different domains
    X1          = abs(x)<x1;
    X2          = (x>=-x2 & x<=-x1)| (x>=+x1 & x<=+x2);
    X3          = abs(x)>x2;

    x_Nstep     = length(x);
    Qx          = NaN(1,x_Nstep);
    DQx         = NaN(1,x_Nstep);

    % longitudinal flow rate    
    Qx(X1)      = Qcs*x(X1)/xd;
    Qx(X2)      = (Qcs/xd-Qr/(x2-x1))*x(X2)+sign(x(X2)).*(Qr/(x2-x1)).*x1;
    Qx(X3)      = Qcs*(x(X3)/xd)+sign(x(X3))*(Qd-Qcs);
    
    % differential of longitudinal flow rate  
    DQx(X1)     = Qcs/xd;
    DQx(X2)     = (Qcs/xd-Qr/(x2-x1));
    DQx(X3)     = Qcs/xd;    
end
end
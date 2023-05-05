function [cs,cvl,cvr,cs_int] = iteration_1DM(cs,cvl,cvr,cs_ghost,cs_ghost_p1,cs_ghost_m1,ss,svl,svr,x1,x2,xd,X_L,X_R,X_HVAC_L,X_HVAC_R,x_step,t_step,A,Vv,F,kd,Diff,theta,Qe,Qcv,Qr,Qd,Qx,DQx,cs_int) 
    
    % source
    cs  = cs + ss*t_step/A; 
    cvl = cvl + svl*t_step/Vv; 
    cvr = cvr + svr*t_step/Vv;

    % advection & diffusion 
    Ablocked = 0;
    [cs]=adv_diff(cs,cs_ghost,cs_ghost_m1,cs_ghost_p1,Diff,t_step,x_step,Qcv,Qx,DQx,A,Ablocked);
    
    % removal saloon 
    KRs = Qr/(A*(x2-x1))*t_step; % reduction factor HVAC saloon
    MRl = A*KRs*sum(cs(X_HVAC_L))*x_step; % mass removed at HVAC inlet
    MRr = A*KRs*sum(cs(X_HVAC_R))*x_step; 
    
    Rcsl = (1-F)*theta*MRl/(A*xd); % concentraiton increase by HVAC outlet
    Rcsr = (1-F)*theta*MRr/(A*xd);
    
    Rcvl = (1-F)*(1-theta)*MRl/Vv; % concentraiton increase by HVAC outlet vestibule
    Rcvr = (1-F)*(1-theta)*MRr/Vv; 
    
    Dcvl = Qd*cs(1)  /Vv*t_step; % concentration increase vestibule through door
    Dcvr = Qd*cs(end)/Vv*t_step;
    
    cvl = cvl+Rcvl+Dcvl; % add to vestibule
    cvr = cvr+Rcvr+Dcvr;
    
    cs(X_HVAC_L) = cs(X_HVAC_L)*(1-KRs); % remove at HVAC inlet
    cs(X_HVAC_R) = cs(X_HVAC_R)*(1-KRs);
    
    cs(X_L) = cs(X_L)+Rcsl; % add at HVAC outlet 
    cs(X_R) = cs(X_R)+Rcsr; 
   
    % removal and decay vestibule
    KRv     = (Qe/Vv+kd)*t_step; % reduction factor vestibule
    cvl     = cvl*(1-KRv);
    cvr     = cvr*(1-KRv);
    
    % decay saloon
    cs     = cs*(1-kd*t_step);  
    
    if exist('cs_int','var')
        % integral
        cs_int = cs_int + cs*t_step;
    end
    
end
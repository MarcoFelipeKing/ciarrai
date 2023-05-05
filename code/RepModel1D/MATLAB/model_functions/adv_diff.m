function [cs]=adv_diff(cs,cs_ghost,cs_ghost_m1,cs_ghost_p1,Diff,t_step,x_step,Qcv,Qx,DQx,A,Ablocked)
    cs_i              = cs;
    cs_ghost(2:end-1) = cs_i;
    cs_ghost(1)       = cs_ghost(3);
    cs_ghost(end)     = cs_ghost(end-2);

    cs_ghost_p1(1:end-1)   = cs_ghost(2:end);
    cs_ghost_i             = cs_ghost;
    cs_ghost_m1(2:end)     = cs_ghost(1:end-1);

    cs_p1 = cs_ghost_p1(2:end-1);
    cs_m1 = cs_ghost_m1(2:end-1);

    dcdx1 = (cs_p1-cs_m1)/(2*x_step);
    dcdx2 = (cs_p1-2*cs_i+cs_m1)/(x_step^2);

    u       = Qx/(A-Ablocked);
    dudx1   = DQx/(A-Ablocked);
    
    on_vent  = Qcv ~= 0; % check if ventilation is turned on
    if on_vent == true
        dcdt = Diff*dcdx2 - dudx1.*cs_i - u.*dcdx1;
    elseif on_vent == false
        dcdt = Diff*dcdx2;
    end
    cs = cs_i + dcdt*t_step;
end
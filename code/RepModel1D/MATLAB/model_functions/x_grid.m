function [x,x_Nstep,x_step,X_L,X_R,X_HVAC_L,X_HVAC_R,x_st,x_idx_st] = x_grid(x1,x2,xd,x_Nstep_half,x_step_st_coarse)
% calculation grid
x_half  = linspace(0,xd,x_Nstep_half); % make grid for half of saloon
x       = [-flip(x_half(2:end)) x_half]; % grid complete saloon
x_Nstep = length(x); % gridpoints
x_step  = mean(diff(x));  % grid step

% storage grid (to store data)
x_half_idx_st = (x_Nstep_half-1) + (1:x_step_st_coarse:x_Nstep_half); % make grid for half of saloon
x_idx_st      = [x_Nstep+1-flip(x_half_idx_st(2:end)) x_half_idx_st]; % idx caculation grid identical to storage grid
if x_idx_st(1) ~= 1 % include grid boundaries 
    x_idx_st = [1 x_idx_st];
end
if x_idx_st(end) ~= x_Nstep
    x_idx_st = [x_idx_st x_Nstep];
end
x_st = x(x_idx_st); % get grid point locations 

% grid domains
X_L         = x<0; % south half
X_R         = x>0; % north half
X_HVAC_L    = x>=-x2 & x<=-x1; % HVAC south
X_HVAC_R    = x>=+x1 & x<=+x2; % HVAC north
end
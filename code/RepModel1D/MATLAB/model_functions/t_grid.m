function [t,t_Nstep,t_idx_st,t_idx_st_f,t_st] = t_grid(t_start,t_step,t_end,t_step_st) % make time series
% calculation time series
t           = t_start:t_step:t_end; % time series
t_Nstep     = length(t); % total time timesteps 

% storage time series
if exist('t_step_st','var')
    t_idx_st     = 1:round(t_step_st/t_step):length(t); % storage index time series
    t_st         = t(t_idx_st); % time series storage
    t_idx_st_num = 1:length(t_idx_st); 
    t_idx_st_f   = NaN(size(t));
    t_idx_st_f(t_idx_st) = t_idx_st_num; % link calculation and storage idx
end
end
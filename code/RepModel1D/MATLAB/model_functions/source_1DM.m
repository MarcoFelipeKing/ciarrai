function [ss] = source_1DM(x,x_step,x_Nstep,xs,sigma,ms) % get source profile
ss      = zeros(1,x_Nstep); % set source to zero
for i = 1:length(xs) % include every source
    xs_i      = xs(i); % source location
    f_i       = 1/(sigma*sqrt(2*pi))*exp(-1/2*((x-xs_i)/sigma).^2); % normal distribution
    f_norm_i  = f_i/(sum(f_i)*x_step); % normalize
    ss_i      = f_norm_i*ms; % set source strength 
    ss        = ss + ss_i; % add source to source term
end
end
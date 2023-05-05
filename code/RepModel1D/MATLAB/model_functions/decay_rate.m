function [kd] =  decay_rate(tracer)

if tracer == "covid"
    % Derive decay constant (https://www.tandfonline.com/doi/full/10.1080/02786826.2020.1829536)
    S    = 0; % integrated UVB irradiance, in W/m^2
    RH   = 40; % relative Humidity in %
    T    = 20; % temperature in degC
    dCS  = 0.05; % width of HVAC outlet

    kd = 7.569 + 1.411*((T-20.54)/10.66)...
                + 0.022*((RH-45.235)/28.665)...
                + 7.553*((S-50)/50)...
                + 1.397*((T-20.54)/10.66)*((S-50)/50); % per minitue 
    kd = kd/60;
else
    kd = 0;
end


end
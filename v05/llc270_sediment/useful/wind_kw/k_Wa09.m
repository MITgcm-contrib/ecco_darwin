function k = k_Wa09(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Wanninkhof et al. (2009)
%k660 = 3. + (0.1 * U) + (0.064 * U^2) + (0.011 * U^3)
 
%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = (3.0 + 0.1 .* windSpeed + 0.064 .* windSpeed .^ 2 + 0.011 .* windSpeed .^ 3) ...
    .* (660 ./ Sc) .^ 0.5;
    
end
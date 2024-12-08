function k = k_Mc01(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of McGillis et al. (2001)
%k660 = 3.3 + (0.026 * U^3)

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = 3.3 + (0.026 .* windSpeed .^ 3) .* (660 ./ Sc) .^ 0.5;

end
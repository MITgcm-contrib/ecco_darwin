function k = k_Ho06(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Ho et al (2006)
%k600 = 0.266 * U^2

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = (0.266 .* windSpeed .^ 2) .* (600 ./ Sc) .^ 0.5;

end
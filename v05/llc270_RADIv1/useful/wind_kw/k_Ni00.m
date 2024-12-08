function k = k_Ni00(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Nightingale et al (2000)
%k600 = (0.333 * U) + (0.222 * U^2)

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = (0.333 .* windSpeed + 0.222 .* windSpeed .^ 2) .* (600 ./ Sc) .^ 0.5;

end
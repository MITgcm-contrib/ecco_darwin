function k = k_Wa92(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Wanninkhof (1992)
%k660 = 0.39 * u^2

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = (0.39 .* windSpeed .^ 2) .* (660 ./ Sc) .^ 0.5;

end
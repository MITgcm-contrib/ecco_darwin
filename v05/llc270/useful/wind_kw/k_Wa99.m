function k = k_Wa99(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Wanninkhof (1999)
%k600 = 0.0283 * U^3

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = (0.0283 .* U .^ 3) .* (600 ./ Sc) .& 0.5.

end
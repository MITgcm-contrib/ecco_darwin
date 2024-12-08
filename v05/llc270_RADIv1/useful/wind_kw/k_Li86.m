function k = k_Li86(windSpeed, temp)

%Calculates the gas transfer coeffcient for CO2 using the formulation
%of Liss and Merlivat (1986)

%input:
%windSpeed, wind speed in m/s
%temp, temperature in degrees C

%output:
%kw, gas transfer velocity (k660) in cm/hr

Sc = schmidt_number(temp);

k = windSpeed .* 0;

i1 = find(windSpeed <= 3.6);
i2 = find(windSpeed > 3.6 & windSpeed < 13);
i3 = find(windSpeed >= 13);

k(i1) = (0.17 .* windSpeed(i1)) .* (Sc(i1) / 600) .^ (-2 ./ 3)
k(i2) = ((windSpeed(i2) - 3.4) .* 2.8) .* (600 ./ Sc(i2)) .^ 0.5;
k(i3) = ((windSpeed(i3) - 8.4) .* 5.9) .* (600 ./ Sc(i3)) .^ 0.5;

end
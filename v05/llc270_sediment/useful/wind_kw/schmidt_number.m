function Sc = schmidt_number(T)

%Calculates the Schmidt number as defined by Jahne et al. (1987) and listed
%in Wanninkhof (2014) Table 1.

%input:
%temp, temperature in degrees C

%output:
%Sc, schmidt number (dimensionless)

a = 2116.8;
b = -136.25;
c = 4.7353;
d = -0.092307;
e = 0.0007555;

%T = T + 273.15;

Sc = a + b .* T + c .* T .^ 2 + d .* T .^ 3 + e .* T .^ 4;

end
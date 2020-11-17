function a = calculate_scaled_alpha(U2,Sc,ice,scaling)

mask = find(isnan(U2) | isnan(Sc));

weight = 1 - ice;
weight(mask) = 0;

a = round(scaling ./ (nansum(nansum(U2 .* Sc .* weight)) ./ nansum(nansum(weight))), 4);

end
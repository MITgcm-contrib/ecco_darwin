function [lat lon interpField] = interp_grid_lat_lon(field,grid)

interpWeightsDir = '/nobackup/dcarrol2/for_jorge/mat/interp_weights/';

if (strcmp(grid,'cs510'))

    load([interpWeightsDir 'interp_weights_cs510_025x025.mat']);
    
elseif (strcmp(grid,'llc270'))
    
    load([interpWeightsDir 'interp_weights_llc270_025x025_deg.mat']);
    
else
    
	return;
    
end

z = field(:).';
zi = sum(z(tri) .* w,2);

interpField = reshape(zi,siz);

end

function output = rapid_mCDR_model(input)

%integrate rapid-mCDR prognostic equations for dAlk and dDIC in time
%time step is defined in the input data

%1s (useful to convert data to seconds)
%dts = seconds(1);

time = input.time;
n_time = length(input.time);

%z_l, depth of the mid layers - first value will be zero
zl = abs(input.Zl);

%z, depth of full layers
z = 0.5 .* (zl(2:end) + zl(1:end-1));
dz = zl(2:end) - zl(1:end-1);
nz = length(z);

%vertical velocity and diffusivity on time*z_l
w = -input.w;
kDiff = input.k_diff;

%all units of dDIC and dAlk are in uM = meq m-3
%1 uM = 10^-3 mol m^-3
dDIC = nan(n_time, nz);
dAlk = nan(n_time, nz);

fCO2 = nan(n_time, 1);

dDIC(1, :) = 0;
dAlk(1, :) = 0;
fCO2(1) = 0;

for t = 1:n_time-1
    
    %integrate from time t --> t+1
    %dt in seconds
    dt = seconds(time(t+1) - time(t));
    
    forcCO2 = -(input.dpco2_ov_ddic(t) * dDIC(t, 1) + input.dpco2_ov_dalk(t) .* dAlk(t, 1)) ...
        .* input.k_surf(t) .* (1 - input.siarea(t));
    
    %alk_forcing is the total forcing, need to divide by area to get it in the per-area units
    forcAlk = input.alk_forcing(t) ./ input.area(t);
    
    dDIC(t+1,:) = solveProgEqn(dDIC(t, :),w(t, :),kDiff(t, :),z,zl,dt,forcCO2);
    dAlk(t+1,:) = solveProgEqn(dAlk(t, :),w(t, :),kDiff(t, :),z,zl,dt,forcAlk);
    
    fCO2(t) = forcCO2;
    
end

output = struct( ...
    'time', time', ...
    'z', z, ...
    'fCO2', fCO2, ...
    'dic', dDIC, ...
    'alk', dAlk, ...
    'dDIC_profile', (dDIC .* input.area2), ...
    'dALK_profile', (dAlk .* input.area2), ...
    'cdr_potential', (input.dpco2_ov_dalk ./ input.dpco2_ov_ddic), ...
    'dpco2_ov_ddic', input.dpco2_ov_ddic, ...
    'dpco2_ov_dalk', input.dpco2_ov_dalk, ...
    'k_surf', input.k_surf, ...
    'area', input.area, ...
    'alk_forcing', input.alk_forcing, ...
    'dz', dz);

end

%%

function varNew = solveProgEqn(var,w,k,z,zl,dt,forcSurf)

a = zeros(size(var));
b = zeros(size(var));
c = zeros(size(var));

dt = seconds(dt);

d = var;
d(1) = d(1) + forcSurf .* dt ./ (zl(2) - zl(1));

A = dt ./ (zl(3:end-1) - zl(2:end-2));
f1 = A .* w(2:end-2);
f2 = A .* w(3:end-1);

%factors for upwind approximation (alphaK,alphaKp1 = 0.5 --> central differencing)
alphaK = heaviside2(w(2:end-2),0);
alphaKp1 = heaviside2(w(3:end-1),0);

f3 = k(2:end-2) * dt ./ (z(2:end-1) - z(1:end-2)) ./ (zl(3:end-1) - zl(2:end-2));
f4 = k(3:end-1) * dt ./ (z(3:end) - z(2:end-1)) ./ (zl(3:end-1) - zl(2:end-2));

a(2:end-1) = -f1 .* alphaK - f3;
b(2:end-1) = 1 + f2 .* alphaKp1 - f1 .* (1 - alphaK) + f3 + f4;
c(2:end-1) = f2 .* (1 - alphaKp1) - f4;

%boundary conditions
%zero tendency at the bottom
b(end) = 1;

%zero flux at the surface
alpha_0 = heaviside2(w(2), 0);

f1s = dt ./ (zl(2) - zl(1)) * w(2);
f2s = k(2) * dt / (z(2) - z(1)) ./ (zl(2) - zl(1));
b(1) = 1 + f1s * alpha_0 + f2s;
c(1) = f1s * (1 - alpha_0) - f2s;

varNew = TDMASolve(a, b, c, d);

end

%%

function x = TDMASolve(a, b, c, d)

n = length(a);

ac = a;
bc = b;
cc = c;
dc = d;

for j = 2:n
    
    if bc(j - 1) == 0
        
        error('Error: Zero pivot element');
        
    end
    
    ac(j) = ac(j) ./ bc(j-1);
    bc(j) = bc(j) - ac(j) * cc(j-1);
    
end

for j = 2:n
    
    dc(j) = dc(j) - ac(j) * dc(j-1);
    
end

dc(n) = dc(n) ./ bc(n);

for j = n-1:-1:1
    
    dc(j) = (dc(j) - cc(j) * dc(j+1)) ./ bc(j);
    
end

x = dc;

end

%%

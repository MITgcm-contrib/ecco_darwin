function [kelp_fr, kelp_ar] = growth_v5(kelp_fr,kelp_ar,envt,farm,time,envt_step,growth_step)
% Calculate uptake and growth to derive delta_Ns and delta_Nt
% Input: kelp_fr, envt, farm, time, envt_step, growth_step
% Output: kelp_fr at t+1


global param
                
            
%% KELP - Known Variables

    Ns = kelp_fr.Ns;
    Nf = kelp_fr.Nf;
    Nf_capacity = kelp_fr.Nf_capacity;
    Age = kelp_fr.Age;
    kelploc = kelp_fr.kelploc;
                
                
%% UPTAKE
% Uptake already determined during uptake loop which happens faster than
% the growth loop in order to update transport function. The average
% rate iss tored in UptakeNavg 
% [mg N/g(dry)/h]: NO3+NH4+Urea
% Weighted since uptake not calculated in uniform increments

    Uptake = kelp_fr.UptakeNavg;
    
    
%% GROWTH
% Growth limited by internal nutrients, temperature, PAR, and height (stops
% growing once reaches max height); [h-1]

    Growth =  glim_v1(kelp_fr.Q,kelp_fr.Type,kelp_fr.Height_tot,envt,farm,envt_step,kelploc);
           
           
%% MORTALTIY                
% Mortality_wave = frond loss due to waves; dependent on Hs, significant
% wave height [m]; Rodrigues et al. 2018 demonstrates linear relationship
% between Hs and frond loss rate [h-1] (continuous)

    M_wave  = param.d_wave_m .* envt.Hs(1,envt_step);

% d_blade = blade erosion [h-1] (continuous); Multiplied by the
% frBlade -> fraction of total as blade
    
    M_blade = param.d_blade .* kelp_fr.frBlade;                 
           
%% Nf at t+1
            
% Nf(t+1) = Nf(t) + Growth - Mortality

     dNf1 = Growth .* Ns .* time.dt_Gr ...
          - M_blade .* Nf .* time.dt_Gr ...
          - M_wave .* Nf .* time.dt_Gr;

     Nf_new = sum(cat(3,Nf,dNf1),3);

                 
%% Frond Characteristics, t+1          

% Nf_Capacity is based on height at current time step but is important for
% allometric conversions at t+1 [mg N/g(dry)/m frond]

% Nf_capacity changes when frond transitions from subsurface to canopy
% type. Create a smoothing function based on height that slowly
% transitions from subsurface to watercolum allometry. Otherwise there
% will be an abrupt redistribution of Nf upwards.

    Nf_capacity_new = (param.Nf_capacity_subsurface-param.Nf_capacity_watercolumn)./(1+exp(1.5.*(kelp_fr.Height_tot-(farm.z_cult*1.2))))+ param.Nf_capacity_watercolumn;
  

% Calculate Age at t+1; 9999, indicator for "harvest" cut [hours]

    Age_new = Age + time.dt_Gr;
    Age_new(Age == 9999) = 9999;
                            
                
%% APICAL GROWTH
% Nf redistributed upwards if Nf > Nf_capacity. Evaluate each depth bin
% separately from the bottom towards the surface. The surface is left alone
% so that the canopy accumulates. The surface bin is z=1.

        
    for z = farm.z_cult:-1:2 

        % delNf is the amount of biomass greater than
        % Nf_capacity. Nf_capacity was previously derived and
        % is dependent on height as kelp transitions from
        % subsurface to canopy.
        delNf = NaN(size(Nf_new(:,z)));       
        delNf(Nf_new(:,z) > Nf_capacity) = Nf_new(Nf_new(:,z) > Nf_capacity,z) - Nf_capacity(Nf_new(:,z) > Nf_capacity); % then calculate the amount of Nf beyond capacity...
        delNf(Nf_new(:,z) <= Nf_capacity) = 0;

        % redistribute Nf based upon a carrying capacity
        Nf_new(delNf>0,z-1) = nansum([delNf(delNf>0) Nf_new(delNf>0,z-1)],2);
        Nf_new(delNf>0,z) = Nf_capacity_new(delNf>0);

    end
    clear z delNf
                
                
%% Ns at t+1
% Ns(t+1) = Ns(t) + Uptake - Growth - Mortality 
% For uptake, only biomass in blades (frBlade) contributes 

    dNs1 = Uptake .* kelp_fr.B .* kelp_fr.frBlade .* time.dt_Gr ...
         - Growth .* Ns .* time.dt_Gr ...
         - param.d_dissolved .* Ns .* time.dt_Gr ...
         - M_blade .* Ns .* time.dt_Gr ...
         - M_wave .* Ns .* time.dt_Gr;

    dNs1_sum = sum(cat(3,Ns,dNs1),3);

% redustribute (translocation) as a function of fractional Nf

    fNf = Nf_new ./ nansum(Nf_new,2);
    Ns_new = nansum(dNs1_sum,2) .* fNf;


%% WHOLE-FROND SENESCENCE
            
% If Type is senescing (Type == 3), Nf and Ns decrease at a rate of
% moartlity_frond;

    Ns_new(kelp_fr.Type==3,:) = Ns(kelp_fr.Type==3,:) - Ns(kelp_fr.Type==3,:) .* param.d_frond .* time.dt_Gr;
    Nf_new(kelp_fr.Type==3,:) = Nf(kelp_fr.Type==3,:) - Nf(kelp_fr.Type==3,:) .* param.d_frond .* time.dt_Gr;


% If Type == senescing && Nf and Ns below threshold set Nf and Ns to NaN
% Threshold set to be 10% of Nf_capacity

    Nf_new(kelp_fr.Type==3  & Nf(:,farm.z_cult) < 0.1 * param.Nf_capacity_subsurface,:) = NaN;
    Ns_new(kelp_fr.Type==3  & Nf(:,farm.z_cult) < 0.1 * param.Nf_capacity_subsurface,:) = NaN;
    Age_new(kelp_fr.Type==3 & Nf(:,farm.z_cult) < 0.1 * param.Nf_capacity_subsurface,:) = NaN;


%% DON, PON
% rDON, the rate of kelp contribution to DON; [mg N/h] -> [mmol N/m3/h]
%   There are four sources of DON contribution from the Ns pool
%   1. dissolved loss
%   2. blade erosion
%   3. wave-based mortality
%   2. senescence

    rDON_new = ...
       (  nansum(param.d_dissolved .* Ns) ...
       +  nansum(M_blade .* Ns) ...
       +  nansum(M_wave .* Ns) ...
       +  nansum(param.d_frond .* Ns(kelp_fr.Type==3,:))) ...
       ./ param.MW_N;
  
% rPON, the rate of kelp contribution to PON; [mg N/m3/h]
%   There are three sources of PON contribution from the Nf pool
%   1. blade erosion
%   2. wave-based mortality
%   3. senescence

    rPON_new = ...
       (  nansum(param.d_blade .* Nf) ...
       +  nansum(M_wave .* Nf) ...
       +  nansum(param.d_frond .* Nf(kelp_fr.Type==3,:)));

                          
%% NEW FROND
% A new frond is initiated at a rate of Frond_init, and happens as a
% discrete event every 1/Frond_init (hours), dependent on Q

% A new frond = Nf equivalent of 1 m (Nf_capacity) initiated at cultivation
% depth

    Si = 1 / (param.Frond_init(1) * nanmean(kelp_fr.Q) + param.Frond_init(2));
    
    % Just in case Q is > 40; but it shouldn't be ...
    if nanmean(kelp_fr.Q) > 40
        Si = 1/ (param.Frond_init(1) * 40 + param.Frond_init(2));
    end

% Evaluate whether or not it is time to start a new frond
% if Current hours is greater than initiation of the "lastFrond" -> YES

    if growth_step * time.dt_Gr >= kelp_fr.lastFrond + Si
        [Nf_new, Ns_new, Nf_capacity_new, Age_new, kelp_fr.ID] = Si_v1(kelp_fr.Q,kelp_fr.ID,Nf_new,Ns_new,Nf_capacity_new,Age_new,farm);
        kelp_fr.lastFrond = growth_step .* time.dt_Gr; % cumulative number of fronds
    end
    clear Si         

    
%% UPDATE STATE VARIABLES
% Note, these replace existing matrices (don't append)

kelp_fr.Nf = Nf_new;
kelp_fr.Ns = Ns_new;
kelp_fr.Nf_capacity = Nf_capacity_new;
kelp_fr.Age = Age_new;

kelp_ar.rDON(kelploc(1),kelploc(2),1:farm.z_cult) = rDON_new;
kelp_ar.rPON(kelploc(1),kelploc(2),1:farm.z_cult) = rPON_new;


end
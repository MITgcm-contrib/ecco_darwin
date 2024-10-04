#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 07:03:32 2021

@author: blsaenz
"""
import copy


# These parameters are for 24-hour timesteps
# In V10, to run sub-daily (1 hour and 3 hour was tested), make these transforms:
# p['mp_spp_kcap'] *= 0.75 #1600 #2300 #1600 #p['mp_spp_kcap']*0.466667
# p['mp_spp_Gmax_cap'] *= 2 (tropical), or *= 1.5 (temperate)

Eucheuma = {
    'spp': 'Eucheuma',         
    'mp_spp_Vmax': 9.2*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below - gMACMODS Jan 2022
    'mp_spp_Ks_NO3': 5600,        # [umol N/m3]  % L. Roberson  - gMACMODS Jan 2022
    #'mp_spp_Vmax': 11.9*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below
    #'mp_spp_Ks_NO3': 11250,        # [umol N/m3]  %
    'mp_spp_kcap': 2000.,         # [g(dry)/m] # 2023 paper: 2000
    'mp_spp_Gmax_cap': 0.2,       # [1/day]
    'mp_spp_PARs': 125.9/4.57,       # [W/m2]
    'mp_spp_PARc': 13.5/4.57,      # [W/m2]
    'mp_spp_Q0': 25.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 5.76,           # [mg N/g(dry)]
    'mp_spp_Qmax': 44.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]
    'mp_spp_Topt1': 22.5,         # [deg C]
    'mp_spp_K1': 0.09,            # temp func slope 1
    'mp_spp_Topt2': 27.5,         # [deg C]
    'mp_spp_K2': 0.09,             # temp func slope 2
    'mp_spp_CD':0.5,             # drag coefficient (unitless)
    'mp_spp_dry_sa': 94.8,       # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.01,            # [d-1] % No info specific for Eucheuma
    'mp_spp_seed': 200.0,          # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Eucheuma['mp_spp_Vmax'] = Eucheuma['mp_spp_Vmax'] * Eucheuma['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d


Agardhiella = {
    'spp': 'Agardhiella',         
    'mp_spp_Vmax': 17.9 * 24,     # [umol N/g-DW/d], expecting umol N/m2/d, conversion below - gMACMODS Jan 2022
    'mp_spp_Ks_NO3': 2000,        # [umol N/m3]  % tuned
    'mp_spp_kcap': 2000.,         # [g(dry)/m] # 2023 paper: 2000
    'mp_spp_Gmax_cap': 0.25,      # [1/day]
    'mp_spp_PARs': 125.9/4.57,    # [W/m2]
    'mp_spp_PARc': 13.5/4.57,     # [W/m2]
    'mp_spp_Q0': 25.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 5.0,           # [mg N/g(dry)]
    'mp_spp_Qmax': 35.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]
    'mp_spp_Topt1': 20.5,         # [deg C]
    'mp_spp_K1': 0.09,            # temp func slope 1
    'mp_spp_Topt2': 27.5,         # [deg C]
    'mp_spp_K2': 0.09,            # temp func slope 2
    'mp_spp_CD':0.5,              # drag coefficient (unitless)
    'mp_spp_dry_sa': 40,          # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.005,            # [d-1] % No info specific for Agardhiella
    'mp_spp_seed': 200.0,         # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Agardhiella['mp_spp_Vmax'] = Agardhiella['mp_spp_Vmax'] * Agardhiella['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d



Sargassum = {
    'spp': 'Sargassum',         
    'mp_spp_Vmax': 16.0*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below
    'mp_spp_Ks_NO3': 2950,        # [umol N/m3]  % 
    'mp_spp_kcap': 500.,         # [g(dry)/m]
    'mp_spp_Gmax_cap': 0.2,       # [1/day]
    'mp_spp_PARs': 303.9/4.57,       # [W/m2]
    'mp_spp_PARc': 26./4.57,      # [W/m2]
    'mp_spp_Q0': 25.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 5.76,           # [mg N/g(dry)]
    'mp_spp_Qmax': 44.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]

    #'mp_spp_Topt1': 20.5,         # [deg C]
    #'mp_spp_K1': 0.03,            # temp func slope 1
    #'mp_spp_Topt2': 25.5,         # [deg C]
    #'mp_spp_K2': 0.08,             # temp func slope 2
    'mp_spp_Topt1': 22.5,         # [deg C]
    'mp_spp_K1': 0.09,            # temp func slope 1
    'mp_spp_Topt2': 27.5,         # [deg C]
    'mp_spp_K2': 0.09,             # temp func slope 2

    'mp_spp_CD':0.5,             # drag coefficient (unitless)
    'mp_spp_dry_sa': 333.0,       # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.01,            # [d-1] % No info specific for Eucheuma
    'mp_spp_seed': 50.0,          # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Sargassum['mp_spp_Vmax'] = Sargassum['mp_spp_Vmax'] * Sargassum['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d

Saccharina = {
    'spp': 'Saccharina',         
    'mp_spp_Vmax': 11.8*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below
    'mp_spp_Ks_NO3': 2000.,        # [umol N/m3]  % 
    'mp_spp_kcap': 2000.,         # [g(dry)/m]
    'mp_spp_Gmax_cap': 0.2,       # [1/day]
    'mp_spp_PARs': 76.3/4.57,       # [W/m2]
    'mp_spp_PARc': 15.5/4.57,      # [W/m2]
    'mp_spp_Q0': 32.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 10.18,           # [mg N/g(dry)]
    'mp_spp_Qmax': 54.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]
    'mp_spp_Topt1': 10.0,         # [deg C]
    'mp_spp_K1': 0.03,            # temp func slope 1
    'mp_spp_Topt2': 15.0,         # [deg C]
    'mp_spp_K2': 0.1,             # temp func slope 2
    'mp_spp_CD':0.5,             # drag coefficient (unitless)
    'mp_spp_dry_sa': 58.0,         # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.01,            # [d-1] % No info specific for Eucheuma
    'mp_spp_seed': 50.0,          # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Saccharina['mp_spp_Vmax'] = Saccharina['mp_spp_Vmax'] * Saccharina['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d


Macrocystis = {
    'spp': 'Macrocystis',         
    'mp_spp_Vmax': 12.8*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below
    'mp_spp_Ks_NO3': 10130,        # [umol N/m3]  % 
    'mp_spp_kcap': 2000.,         # [g(dry)/m]
    'mp_spp_Gmax_cap': 0.2,       # [1/day]
    'mp_spp_PARs': 212.4/4.57,       # [W/m2]
    'mp_spp_PARc': 20.45/4.57,      # [W/m2]
    'mp_spp_Q0': 32.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 10.18,           # [mg N/g(dry)]
    'mp_spp_Qmax': 54.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]
    'mp_spp_Topt1': 13.0,         # [deg C]
    'mp_spp_K1': 0.04,            # temp func slope 1
    'mp_spp_Topt2': 18.0,         # [deg C]
    'mp_spp_K2': 0.05,             # temp func slope 2
    'mp_spp_CD':0.5,             # drag coefficient (unitless)
    'mp_spp_dry_sa': 58.0,         # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.01,            # [d-1] % No info specific for Eucheuma
    'mp_spp_seed': 50.0,          # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Macrocystis['mp_spp_Vmax'] = Macrocystis['mp_spp_Vmax'] * Macrocystis['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d


Porphyra = {
    'spp': 'Porphyra',         
    'mp_spp_Vmax': 52.2*24,         # [umol N/g-DW/d], expecting umol N/m2/d, conversion below
    'mp_spp_Ks_NO3': 5200,        # [umol N/m3]  % 
    'mp_spp_kcap': 125.,         # [g(dry)/m]
    'mp_spp_Gmax_cap': 0.2,       # [1/day]
    'mp_spp_PARs': 104./4.57,       # [W/m2]
    'mp_spp_PARc': 24.8/4.57,      # [W/m2]
    'mp_spp_Q0': 32.0,            # initial Q [mg N/g(dry]
    'mp_spp_Qmin': 10.18,           # [mg N/g(dry)]
    'mp_spp_Qmax': 54.0,          # [mg N/g(dry)]
    'mp_spp_BtoSA': 1.0,          # Hmm not used right now???
    'mp_spp_line_sep': 1.0,       # m
    'mp_spp_kcap_rate': 0.05 ,    # [1/day]
    'mp_spp_Topt1': 12.,         # [deg C]
    'mp_spp_K1': 0.03,            # temp func slope 1
    'mp_spp_Topt2': 17.,         # [deg C]
    'mp_spp_K2': 0.09,             # temp func slope 2
    'mp_spp_CD':0.5,             # drag coefficient (unitless)
    'mp_spp_dry_sa': 10.0,         # [g(dry)/m2]
    'mp_spp_dry_wet': 0.094,      # [g(dry)/g(wet)] % Not changed from the macrocystis values
    'mp_spp_E': 0.01,             # [d-1] % No info specific for Eucheuma
    'mp_spp_seed': 10.0,          # initial biomass [g(dry)/m]
    'mp_spp_death': 0.01,         # death rate [1/day]        
}
Porphyra['mp_spp_Vmax'] = Porphyra['mp_spp_Vmax'] * Porphyra['mp_spp_dry_sa'] # umol N/g-DW/d  * g-DW/m2 = umol N/m2/d


def spp_std_harvest_params(spp):
    my_params = {}
    if spp == 'Macrocystis':
        my_params['mp_harvest_type'] = 1    # w/ conditional schedule, harvest _f at end?
        my_params['mp_harvest_kg'] = 1.35 #0.9
        my_params['mp_harvest_freq'] = 220
        my_params['mp_harvest_span'] =  150
        my_params['mp_harvest_schedule'] = 1
        my_params['mp_harvest_nmax'] = 2
        my_params['mp_harvest_f'] = 0.8

    elif spp == 'Saccharina':
        my_params['mp_harvest_type'] = 0      # w/ conditional schedule, harvest 99% at end
        my_params['mp_harvest_kg'] = 1.35 #0.9
        my_params['mp_harvest_freq'] = 180
        my_params['mp_harvest_span'] =  120
        my_params['mp_harvest_schedule'] = 1
        my_params['mp_harvest_nmax'] = 2
        my_params['mp_harvest_f'] = 0.8
         
    elif spp == 'Eucheuma':
        my_params['mp_harvest_type'] = 1      # harvest fractionally
        my_params['mp_harvest_kg'] = 0.8
        my_params['mp_harvest_freq'] = 364
        my_params['mp_harvest_span'] =  320
        my_params['mp_harvest_schedule'] = 2
        my_params['mp_harvest_nmax'] = 8
        my_params['mp_harvest_f'] = 0.8

        # my_params['mp_harvest_type'] = 0      # harvest to seed weight
        # my_params['mp_harvest_kg'] = 0.8      # not used for fixed
        # my_params['mp_harvest_freq'] = 45
        # my_params['mp_harvest_span'] =  0     # not used for fixed
        # my_params['mp_harvest_schedule'] = 0  # fixed
        # my_params['mp_harvest_nmax'] = 8      # not used for fixed, but will be 8
        # my_params['mp_harvest_f'] = 0.8       # not used for harvest_type 0
        # my_params['mp_spp_kcap_rate'] = 0.05          
        # my_params['mp_spp_kcap'] = 2000.0

    elif spp == 'Agardhiella':
        my_params['mp_harvest_type'] = 1      # harvest fractionally
        my_params['mp_harvest_kg'] = 0.8
        my_params['mp_harvest_freq'] = 364
        my_params['mp_harvest_span'] =  320
        my_params['mp_harvest_schedule'] = 2
        my_params['mp_harvest_nmax'] = 8
        my_params['mp_harvest_f'] = 0.8

        # my_params['mp_harvest_type'] = 0      # harvest to seed weight
        # my_params['mp_harvest_kg'] = 0.8      # not used for fixed
        # my_params['mp_harvest_freq'] = 45
        # my_params['mp_harvest_span'] =  0     # not used for fixed
        # my_params['mp_harvest_schedule'] = 0  # fixed
        # my_params['mp_harvest_nmax'] = 8      # not used for fixed, but will be 8
        # my_params['mp_harvest_f'] = 0.8       # not used for harvest_type 0
        # my_params['mp_spp_kcap_rate'] = 0.05          
        # my_params['mp_spp_kcap'] = 2000.0
        
    elif spp == 'Porphyra':
        my_params['mp_harvest_type'] = 1    
        my_params['mp_harvest_kg'] = 0.08
        my_params['mp_harvest_freq'] = 150
        my_params['mp_harvest_span'] =  110
        my_params['mp_harvest_schedule'] = 1
        my_params['mp_harvest_nmax'] = 6
        my_params['mp_harvest_f'] = 0.8

        my_params['mp_spp_seed'] = 10.0

    elif spp == 'Sargassum':
        # my_params['mp_harvest_type'] = 0      # harvest to seed weight
        # my_params['mp_harvest_kg'] = 0.8      # not used for fixed
        # my_params['mp_harvest_freq'] = 45
        # my_params['mp_harvest_span'] =  0     # not used for fixed
        # my_params['mp_harvest_schedule'] = 0  # fixed
        # my_params['mp_harvest_nmax'] = 8      # not used for fixed, but will be 8
        # my_params['mp_harvest_f'] = 0.8       # not used for harvest_type 0
        # my_params['mp_spp_kcap_rate'] = 0.05          
        # my_params['mp_spp_kcap'] = 800.0

        my_params['mp_harvest_type'] = 1      
        my_params['mp_harvest_kg'] = 0.4
        my_params['mp_harvest_freq'] = 364
        my_params['mp_harvest_span'] =  320
        my_params['mp_harvest_schedule'] = 2
        my_params['mp_harvest_nmax'] = 8
        my_params['mp_harvest_f'] = 0.8

    return my_params


spp_p_dict = {'Eucheuma':Eucheuma,'Sargassum':Sargassum,
              'Saccharina':Saccharina,'Macrocystis':Macrocystis,
              'Porphyra':Porphyra,'Agardhiella':Agardhiella}


spp_p_dict_harvest = copy.deepcopy(spp_p_dict)
for spp,p_dict in spp_p_dict_harvest.items():
    spp_p_dict_harvest[spp].update(spp_std_harvest_params(spp))





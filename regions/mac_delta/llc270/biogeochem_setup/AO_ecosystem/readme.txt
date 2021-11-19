The new Darwin ecosystem has been created to simulate general plankton dynamics taking place in the Arctic Ocean.
The new ecosystem setup includes the RADTRANS package and allows for CDOM.

Use git clone -b cdom-carbon https://github.com/jahn/darwin3

This set up involves 34 passive tracers (Only 33 in version 0 excluding RDOC):

1. DIC     13. POC    25. c4
2. NO3     14. PON    26. c5
3. NO2     15. POP    27. c6
4. NH4     16. POFe   28. c7
5. PO4     17. POSi   29. c8
6. FeT     19. PIC    30. Chl1
7. SiO2    19. ALK    31. Chl2
8. DOC     20. O2     32. Chl3
9. RDOC    21. CDOM   33. Chl4
10. DON    22. c1     34. Chl5
11. DOP    23. c2
12. DOFe   24. c3

This ecosystem contains:
  
  > 5 phytoplankton types:
      - Pico-eukaryote (analog of Micromonas spp - low-light adapted).
        	Nominal size: 3um ESD, quantum yield = 6e-5 mmol C mg Chl^-1 s^-1.
      - Haptophyte (which could be Phaeocystis).
        	Nominal size: 4.5um ESD, quantum yield = 4e-5.
      - High-light adapted diatoms (analog of centric diatom – potentially Thalassioira or Chaetoceras spp).
        	Nominal size: 10um ESD, quantum yield = 4e-5
      - Low-light adapted diatom (analog of pennate diatom – potentially Fragilariopsis spp)
        	Nominal size: 10um ESD, quantum yield = 6e-5; growth rate 0.8 * diatom(3) as cost for LL adaption, and inhibition value = 1.2.
      - Dinoflagellate (which will be mixotrophic, feeding on pico-euk and bacteria).
        	Nominal size: 10um ESD, quantum yield=4e-5.
  
  > 2 zooplankton types:
      - Small zooplankton (will feed on pico-euk, and bacteria).
        	Nominal size 45um ESD.
      - Large zooplankton (will feed on diatoms, haptophyte and dinoflagellate).
        	Nominal size 100um.
        
  > 1 bacteria:
      - Free-living heterotrophic bacteria.
        	Nominal size 0.6um ESD.

First approach adopted for (unknown) initial/boundary conditions
  
  > Tracers 1 to 20: use same initial/boundary conditions as Carroll et al. (2020) with RDOC added.
  	
  	*For tracer 9 RDOC: use IC of 50 uM everywhere and BCs with the same value. The Carroll et al. (2020) ecosystem
  	is then run for 14 years (2003–2017). Then the last RDOC field from the simulation is used as the IC for a second
  	14-year simulation with same BCs. The last field from the second simulation is used as the final IC.  		
  	
  > Tracer 21 (CDOM): use Carroll et al. (2020) DOP initial/boundary conditions.
  	Units should be converted from mmol P m^-3 to mmol C m^-3 using R_PC ratio (120/1).
  
  > Tracers 22, 24, 25, and 26 (phytoplankton): use sum of 5 phytoplankton from Carroll et al. (2020) and divide by 4.
    Set negative values to 0.
  
  > Tracer 23: use initial/boundary conditions set to 0 (Haptophytes are not considered in the Mackenzie shelf)
  
  > Tracer 27 and 28 (zooplankton): use the sum of 2 zooplankton fields from Carroll et al. (2020) divided by 2.
  	Set negative values to 0.
  	
  > Tracer 29 (bacteria): field is set constant at 0.1 mmol N m^-3 from 0-100m and 0.001 mmol N m^-3 from 100m to the bottom. 
    See  Le Fouest et al. 2013 and Le Fouest et al. 2013
    Field in mmol N m^-3 is then converted into mmol C m^-3 using the R_CN ratio (120/16).
  
  > Tracers 30 to 34: use same method as tracers 22-26 but with Chl fields (reminder: Haptophytes not yet considered, therefore Tracer 31 is set to 0).
  
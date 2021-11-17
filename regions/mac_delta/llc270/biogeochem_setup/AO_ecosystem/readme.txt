The new Darwin Ecosystem has been created to simulate the general plankton dynamic taking place in the Arctic Ocean.

The new Ecosystem set up involve RADTRANS package and allow CDOM treatment (upcoming updates).

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

This Arctic Ecosystem contains:
  > 5 phytoplankton types:
      - Pico-eukaryote (analog of Micromonas spp - low light adapted)
        nominal size: 3um ESD, quantum yield=6e-5 mmolC/mgChl/s
      - Haptophyte (which could be Phaeocystis)
        nominal size: 4.5um ESD, quantum yield=4e-5
      - High-light adapted diatoms (analog of centric diatom – potentially Thalassioira or Chaetoceras spp)
        nominal size: 10um ESD, quantum yield=4e-5
      - Low-light adapted diatom (analog of pennate diatom – potentially Fragilariopsis spp)
        nominal size: 10um ESD, quantum yield=6e-5; growth rate 0.8*diatom(3) as cost for LL adaption, and inhibition value 1.2
      - Dinoflagellate (which will be mixotrophic, feeding on pico-euk and bacteria)
        nominal size: 10um ESD, quantum yield=4e-5
  > 2 Zooplankton types:
      - Small zooplankton (will feed on pico-euk, and bacteria)
        nominal size 45um ESD
      - Large zooplankton (will feed on diatoms, haptophyte and dinoflagellate)
        nominal size 100um
  > 1 Bacteria:
      - Free-living heterotrophic bacteria
        nominal size 0.6um ESD

First approach adopted for (unknown) Initial/Boundary Conditions
  > Tracers 1 to 20: Use same initial/Boundary conditions than Carroll et al. (2020) set up modified with tDOC
  > Tracer 21 (CDOM): Use DOP initial/Boundary conditions (the unit should be converted from mmolP/m3 to mmolC/m3)
  > Tracer 22, 24, 25 & 26 (Phytoplankton): Use sum of 5 PFTs of the Carroll et al. (2020) set up divide by 4.
                                                  Absolute value of the 5 PFTs fields is used to get rid of very small existing values in its.
  > Tracer 23 : Use initial/Boundary conditions set to 0 on the whole volume (Haptophytes are not considered in the Mackenzie shelf)
  > Tracer 27 & 28 (Zooplankton): Use sum of 2 zoo's (absolute value) fields of the Carroll et al. (2020) set up divide by 2.
  > Tracer 29 (Bacteria): Field is set constant at 0.1 mmolN/m3 [0-100m] & 0.001 mmolN/m3 [100m-bottom] (see Le Fouest et al. 2015, Le Fouest et al. 2013)
                          Field in mmolN/m3 is then converted into mmolC/m3 using R_CN ratio set to 120/16
  > Tracer 30 to 34: Use same method as PFTs concentration with Chlorophyll fields (Reminder : Haptophytes not considered yet Tracer 3à set to 0)

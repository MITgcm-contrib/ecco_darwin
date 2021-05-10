ARCTIC ECOSYSTEM README FILE:

The new Darwin Ecosystem has been created to simulate the general plankton dynamic taking place in the Arctic Ocean.
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

The new Ecosystem set up involve RADTRANS package and CDOM treatment (upcoming updates)

I. First approach adopted for (unknown) Initial/Boundary Conditions
  a. Initial conditions
      - Phytoplankton: Use the same initial condition fiels for 5 phyto types (Haptophytes are intially not considered s set to 0)
                       The field is set constant at 0.2 mmolN/m3 from 0 to 100m and 0.002 from 100m to bottom (see Le Fouest et al. 2015)
                       the field was converted into mmolC/m3 using R_CN ratio set to 120/16
      - Chlorophyll: convert Phytoplankton field from mmolC/m3 into mmolChl/m3 using R_ChlC ratio set to 16/120
      - Zooplankton & Bacteria: Use the same initial condition field for 2 zoo types and Bacteria
                                The field is set constant at 0.1 mmolN/m3 from 0 to 100m and 0.001 from 100m to bottom (see Le Fouest et al. 2015, Le Fouest et al. 2013)
                                the field was converted into mmolC/m3 using R_CN ratio set to 120/16
      - CDOM: Use DOP initial conditions (the unit is the same mmolP/m3)

  b. Boundary conditions
      - Phytoplankton: Sum the concentration of the 5 Phytoplankton fields of old ecosystem
                       Divide the sum by 4 to distribute it on the 4 new phytoplankton types (Reminder: Haptophytes not considered yet)
      - Chlorophyll: Sum the concentration of the 5 Chlorophyll fields of old ecosystem
                     Divide the sum by 4 to distribute it on the 4 new phytoplankton types (Reminder: Haptophytes not considered yet)
      - Zooplankton: Sum the concentration of the 2 zooplankton fields of old ecosystem
                     Divide the sum by 2 to distribute it on the 2 new zooplankton types
      - Bacteria: Set same Boundary conditions as initial field
      - CDOM: Use DOP boundary conditions

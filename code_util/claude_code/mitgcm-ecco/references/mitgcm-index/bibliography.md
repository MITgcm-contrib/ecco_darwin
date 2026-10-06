# Bibliography: papers behind MITgcm / darwin3 schemes

Generated from the manuals' `manual_references.bib` (MITgcm, plus darwin3's extra ecosystem entries) and
every `:cite:` in the manual. **Part 1** lists, per manual page, which papers it cites, where (`L<line>`), under which
section, and the parameter on that line. **Part 2** gives the full references. `grep -n <param or author> bibliography.md`.
Papers the manuals don't cite (ECCO, ECCO-Darwin) are in `references/literature.md`.

## Part 1: citations by manual page

### darwin3 `doc/phys_pkgs/darwin_airsea.rst`
- L44 `wannink:92` Wanninkhof 1992 — Air-sea exchanges
- L45 `weiss:80` Weiss & Price 1980 — Air-sea exchanges
- L71 `garcia:92` Garcia & Gordon 1992 — Air-sea exchanges
- L72 `keeling:98` Keeling et al. 1998 — Air-sea exchanges

### darwin3 `doc/phys_pkgs/darwin_carbon.rst`
- L30 `munhoven:13` Munhoven 2013 — Carbon chemistry options
- L31 `follows:06` Follows et al. 2006 — Carbon chemistry options
- L34 `sulpis:22` Sulpis et al. 2022 — Carbon chemistry options
- L62 `uppstrom:74` Uppström 1974 — Carbon chemistry options
- L63 `lee:10` Lee et al. 2010 — Carbon chemistry options
- L70 `riley:65` Riley 1965 — Carbon chemistry options
- L71 `culkin:65` Culkin 1965 — Carbon chemistry options
- L78 `dickson:79` Dickson & Riley 1979 — Carbon chemistry options
- L79 `perez:87` Perez & Fraga 1987 — Carbon chemistry options
- L86 `mehrbach:73` Mehrbach et al. 1973 — Carbon chemistry options
- L86 `millero:95` Millero 1995 — Carbon chemistry options
- L87 `roy:93` Roy et al. 1993 — Carbon chemistry options
- L88 `millero:95` Millero 1995 — Carbon chemistry options
- L89 `lueker:00` Lueker et al. 2000 — Carbon chemistry options
- L90 `millero:10b` Millero 2010 — Carbon chemistry options
- L91 `waters:13` Waters & Millero 2013 — Carbon chemistry options
- L91 `waters:14` Waters et al. 2014 — Carbon chemistry options
- L145 `keir:80` Keir 1980 — Calcite dissolution
- L152 `naviaux:19` Naviaux et al. 2019 — Calcite dissolution

### darwin3 `doc/phys_pkgs/darwin_iron.rst`
- L98 `parekh:2005` Parekh et al. 2005 — Scavenging — `DARWIN_PART_SCAV`
- L125 `parekh:2004` Parekh et al. 2004 — Scavenging
- L126 `dutkiewicz:2005` Dutkiewicz et al. 2005 — Scavenging

### darwin3 `doc/phys_pkgs/darwin_spectral.rst`
- L9 `dutkiewicz:2015` Dutkiewicz et al. 2015 — Spectral Light
- L163 `dutkiewicz:2020` Dutkiewicz et al. 2020 — Allometric scaling of absorption and scattering spectra
- L196 `montagnes:1994` Montagnes et al. 1994 — Total scattering

### `doc/algorithm/adv-schemes.rst`
- L52 `adcroft:95` Adcroft 1995 — Centered second order advection-diffusion
- L52 `adcroft:97` Adcroft et al. 1997 — Centered second order advection-diffusion
- L237 `roe:85` Roe 1985 — Second order flux limiters

### `doc/algorithm/algorithm.rst`
- L428 `durran:91` Durrann 1991 — Adams-Bashforth III
- L1534 `adcroft:98` Adcroft & Marshall 1998 — Lateral dissipation
- L1598 `wajsowicz:93` Wajsowicz 1993 — Vertical dissipation
- L1599 `griffies:00` Griffies & Hallberg 2000 — Vertical dissipation
- L2101 `cam:04` Campin et al. 2004 — Time-stepping of tracers: ABII
- L2101 `griffies:00` Griffies & Hallberg 2000 — Time-stepping of tracers: ABII
- L2176 `shapiro:70` Shapiro 1970 — Shapiro Filter
- L2344 `bryan:75` Bryan et al. 1975 — Reynolds-Number Limited Eddy Viscosity
- L2356 `bryan:75` Bryan et al. 1975 — Reynolds-Number Limited Eddy Viscosity — `viscAhReMax`
- L2377 `smag:93` Smagorinsky 1993 — Vertical Eddy Viscosities
- L2384 `smag:63` Smagorinsky 1963 — Smagorinsky Viscosity
- L2384 `smag:93` Smagorinsky 1993 — Smagorinsky Viscosity
- L2417 `smag:93` Smagorinsky 1993 — Smagorinsky Viscosity
- L2450 `griffies:00` Griffies & Hallberg 2000 — Smagorinsky Viscosity — `viscC2Smag`
- L2451 `smag:93` Smagorinsky 1993 — Smagorinsky Viscosity
- L2454 `smag:93` Smagorinsky 1993 — Smagorinsky Viscosity
- L2470 `leith:68` Leith 1968 — Leith Viscosity
- L2470 `leith:96` Leith 1996 — Leith Viscosity
- L2560 `bachman:17` Bachman et al. 2017 — Quasi-Geostrophic Leith Viscosity
- L2584 `bachman:17` Bachman et al. 2017 — Quasi-Geostrophic Leith Viscosity
- L2616 `bachman:17` Bachman et al. 2017 — Quasi-Geostrophic Leith Viscosity — `ALLOW_LEITH_QG`
- L2643 `griffies:00` Griffies & Hallberg 2000 — Courant–Freidrichs–Lewy Constraint on Viscosity
- L2652 `holland:78` Holland 1978 — Biharmonic Viscosity
- L2683 `griffies:00` Griffies & Hallberg 2000 — Biharmonic Viscosity
- L2738 `griffies:00` Griffies & Hallberg 2000 — Biharmonic Viscosity
- L2746 `griffies:00` Griffies & Hallberg 2000 — Biharmonic Viscosity
- L2769 `gill:82` Gill 1982 — Mercator, Nondimensional Equations

### `doc/algorithm/c-grid.rst`
- L7 `arakawa:77` Arakawa & Lamb 1977 — C grid staggering of variables

### `doc/algorithm/finitevol-meth.rst`
- L44 `adcroft:97` Adcroft et al. 1997 — The finite volume method: finite volumes versus finite difference

### `doc/algorithm/nonlinear-freesurf.rst`
- L108 `cam:04` Campin et al. 2004 — Free surface effect on column total thickness (Non-linear free-surface)
- L250 `cam:04` Campin et al. 2004 — Tracer conservation with non-linear free-surface
- L423 `adcroft:04a` Adcroft & Campin 2004 — Non-linear free-surface and vertical resolution

### `doc/algorithm/vert-grid.rst`
- L62 `adcroft:97` Adcroft et al. 1997 — Topography: partially filled cells

### `doc/autodiff/autodiff.rst`
- L11 `griewank:08` Griewank & Walther 2008 — Automatic Differentiation
- L26 `giering:98` Giering & Kaminski 1998 — Automatic Differentiation
- L27 `giering:00` Giering 2000 — Automatic Differentiation
- L29 `maro-eta:99` Marotzke et al. 1999 — Automatic Differentiation
- L30 `stammer:02` Stammer et al. 2002 — Automatic Differentiation
- L30 `stammer:97` Stammer et al. 1997 — Automatic Differentiation
- L35 `naumann:06` Naumann et al. 2006 — Automatic Differentiation
- L35 `utke:08` Utke et al. 2008 — Automatic Differentiation
- L40 `gaikwad:24` Gaikwad et al. 2025 — Automatic Differentiation
- L548 `griewank:92` Griewank 1992 — Storing vs. recomputation in reverse mode
- L549 `restrepo:98` Restrepo et al. 1998 — Storing vs. recomputation in reverse mode
- L1346 `gil-lem:89` Gilbert & Lemarechal 1989 — Control variable handling for optimization applications
- L1755 `gaikwad:24` Gaikwad et al. 2025 — Adjoint code generation using Tapenade
- L1755 `hascoet:24` Hascoet et al. 2024 — Adjoint code generation using Tapenade

### `doc/contributing/contributing.rst`
- L1186 `«BIB_REFERENCE»` (not in bib) — Citations
- L1197 `bryan:79` Bryan & Lewis 1979 — Citations
- L1198 `bryan:79` Bryan & Lewis 1979 — Citations

### `doc/examples/advection_in_gyre/advection_in_gyre.rst`
- L18 `dutay:02` Dutkiewicz et al. 2005 — Ocean Gyre Advection Schemes
- L23 `marshall:06` Marshall et al. 2006 — Ocean Gyre Advection Schemes

### `doc/examples/baroclinic_gyre/baroclinic_gyre.rst`
- L14 `cox:84` Cox & Bryan 1984 — Baroclinic Ocean Gyre
- L86 `marshall:97a` Marshall et al. 1997 — Equations solved
- L229 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L229 `gill:82` Gill 1982 — Numerical Stability Criteria
- L446 `munk:50` Munk 1950 — PARM01 - Continuous equation parameters — `no_slip_sides`
- L448 `stommel:48` Stommel 1948 — PARM01 - Continuous equation parameters
- L1238 `pedlosky:96` Pedlosky 1996 — Model solution
- L1239 `vallis:17` Vallis 17 — Model solution
- L1245 `cushmanroisin:11` Cushman-Roisin & Beckers 2011 — Model solution
- L1245 `vallis:17` Vallis 17 — Model solution

### `doc/examples/barotropic_gyre/barotropic_gyre.rst`
- L10 `bryan:63` Bryan 1963 — Barotropic Ocean Gyre
- L10 `munk:50` Munk 1950 — Barotropic Ocean Gyre
- L10 `stommel:48` Stommel 1948 — Barotropic Ocean Gyre
- L13 `cushmanroisin:11` Cushman-Roisin & Beckers 2011 — Barotropic Ocean Gyre
- L13 `vallis:17` Vallis 17 — Barotropic Ocean Gyre
- L48 `marshall:97a` Marshall et al. 1997 — Equations Solved
- L94 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L111 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L127 `pedlosky:87` Pedlosky 1987 — Numerical Stability Criteria
- L140 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L693 `pedlosky:87` Pedlosky 1987 — Model Solution

### `doc/examples/deep_convection/deep_convection.rst`
- L102 `marshall:97a` Marshall et al. 1997 — Equations solved

### `doc/examples/examples.rst`
- L299 `neale:01` Neale & Hoskins 2001 — Additional Example Experiments: Forward Model Setups
- L333 `stevens:90` Stevens 1990 — Additional Example Experiments: Forward Model Setups
- L347 `neale:01` Neale & Hoskins 2001 — Additional Example Experiments: Forward Model Setups
- L350 `held-suar:94` Held & Suarez 1994 — Additional Example Experiments: Forward Model Setups
- L359 `ferrari:10` Ferrari et al. 2010 — Additional Example Experiments: Forward Model Setups
- L363 `ferrari:08` Ferrari et al. 2008 — Additional Example Experiments: Forward Model Setups
- L378 `gas-eta:90` Gaspar et al. 1990 — Additional Example Experiments: Forward Model Setups
- L394 `gas-eta:90` Gaspar et al. 1990 — Additional Example Experiments: Forward Model Setups
- L419 `held-suar:94` Held & Suarez 1994 — Additional Example Experiments: Forward Model Setups
- L422 `held-suar:94` Held & Suarez 1994 — Additional Example Experiments: Forward Model Setups
- L425 `held-suar:94` Held & Suarez 1994 — Additional Example Experiments: Forward Model Setups
- L435 `klymaklegg10` Klymak & Legg 2010 — Additional Example Experiments: Forward Model Setups
- L443 `hellmer:89` Hellmer & Olbers 1989 — Additional Example Experiments: Forward Model Setups
- L456 `hibler:87` Hibler & Bryan 1987 — Additional Example Experiments: Forward Model Setups
- L559 `munhoven:13` Munhoven 2013 — Additional Example Experiments: Forward Model Setups
- L572 `smag:63` Smagorinsky 1963 — Additional Example Experiments: Forward Model Setups
- L582 `lar-eta:94` Large et al. 1994 — Additional Example Experiments: Forward Model Setups
- L586 `gas-eta:90` Gaspar et al. 1990 — Additional Example Experiments: Forward Model Setups
- L593 `mellor:82` Mellor & Yamada 1982 — Additional Example Experiments: Forward Model Setups
- L596 `pal-rom:97` Paluszkiewicz & Romea 1997 — Additional Example Experiments: Forward Model Setups
- L599 `pacanowski:81` Pacanowski & Philander 1981 — Additional Example Experiments: Forward Model Setups
- L674 `hellmer:89` Hellmer & Olbers 1989 — Additional Example Experiments: Adjoint Model Setups

### `doc/examples/global_oce_biogeo/global_oce_biogeo.rst`
- L14 `trenberth:89` Trenberth et al. 1989 — Overview
- L15 `jiang:99` Jiang et al. 1999 — Overview
- L17 `levitus:94a` Levitus & Boyer 1994 — Overview
- L17 `levitus:94b` Levitus & Boyer 1994 — Overview
- L18 `gen-mcw:90` Gent & McWilliams 1990 — Overview
- L27 `yamanaka:97` Y. & Tajika 1997 — Overview
- L29 `yamanaka:97` Y. & Tajika 1997 — Overview
- L30 `martin:87` Martin et al. 1987 — Overview
- L33 `follows:06` Follows et al. 2006 — Overview
- L35 `wannink:92` Wanninkhof 1992 — Overview
- L37 `wannink:92` Wanninkhof 1992 — Overview
- L38 `dutkiewicz:05` Dutkiewicz et al. 2005 — Overview
- L126 `dutkiewicz:05` Dutkiewicz et al. 2005 — Equations Solved
- L127 `mckinley:04` McKinley et al. 2004 — Equations Solved

### `doc/examples/global_oce_in_p/global_oce_in_p.rst`
- L27 `trenberth:90` Trenberth et al. 1990 — Overview
- L28 `jiang:99` Jiang et al. 1999 — Overview
- L28 `levitus:94a` Levitus & Boyer 1994 — Overview
- L28 `levitus:94b` Levitus & Boyer 1994 — Overview
- L34 `haney:71` Haney 1971 — Overview
- L84 `deszoeke:02` de Szoeke & Samelson 2002 — Discrete Numerical Configuration
- L123 `marshall:97a` Marshall et al. 1997 — Discrete Numerical Configuration
- L124 `cam:04` Campin et al. 2004 — Discrete Numerical Configuration
- L426 `jackett:95` Jackett & McDougall 1995 — File :filelink:`input/data <verification/tutorial_global_oce_in_p/input/data>`
- L427 `mcdougall:03` McDougall et al. 2003 — File :filelink:`input/data <verification/tutorial_global_oce_in_p/input/data>`
- L673 `levitus:94a` Levitus & Boyer 1994 — Files ``input/lev_t.bin`` and ``input/lev_s.bin``
- L673 `levitus:94b` Levitus & Boyer 1994 — Files ``input/lev_t.bin`` and ``input/lev_s.bin``
- L684 `trenberth:90` Trenberth et al. 1990 — Files ``input/trenberth_taux.bin`` and ``input/trenberth_tauy.bin``
- L690 `levitus:94b` Levitus & Boyer 1994 — File ``input/lev_sst.bin``
- L697 `jiang:99` Jiang et al. 1999 — Files ``input/shi_qnet.bin`` and ``input/shi_empmr.bin``

### `doc/examples/global_oce_latlon/global_oce_latlon.rst`
- L14 `bryan:84` Bryan 1984 — Global Ocean Simulation
- L23 `trenberth:90` Trenberth et al. 1990 — Overview
- L24 `kalnay:96` Kalnay et al. 1996 — Overview
- L24 `levitus:94a` Levitus & Boyer 1994 — Overview
- L24 `levitus:94b` Levitus & Boyer 1994 — Overview
- L30 `haney:71` Haney 1971 — Overview
- L117 `marshall:97a` Marshall et al. 1997 — Discrete Numerical Configuration
- L197 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L213 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L239 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L249 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L262 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L440 `jackett:95` Jackett & McDougall 1995 — File :filelink:`input/data <verification/tutorial_global_oce_latlon/input/data>`
- L586 `trenberth:90` Trenberth et al. 1990 — Files ``input/trenberth_taux.bin`` and ``input/trenberth_tauy.bin``

### `doc/examples/global_oce_optim/global_oce_optim.rst`
- L20 `levitus:94a` Levitus & Boyer 1994 — Overview
- L20 `levitus:94b` Levitus & Boyer 1994 — Overview
- L26 `ferriera:05` Ferreira et al. 2005 — Overview
- L26 `stammer:02` Stammer et al. 2002 — Overview
- L41 `levitus:94a` Levitus & Boyer 1994 — Overview
- L63 `stammer:02` Stammer et al. 2002 — Overview
- L72 `gil-lem:89` Gilbert & Lemarechal 1989 — Overview

### `doc/examples/held_suarez_cs/held_suarez_cs.rst`
- L9 `held-suar:94` Held & Suarez 1994 — Held-Suarez Atmosphere
- L11 `adcroft:04a` Adcroft & Campin 2004 — Held-Suarez Atmosphere
- L13 `adcroft:04b` Adcroft et al. 2004 — Held-Suarez Atmosphere
- L25 `held-suar:94` Held & Suarez 1994 — Overview
- L39 `held-suar:94` Held & Suarez 1994 — Overview
- L42 `shapiro:70` Shapiro 1970 — Overview
- L46 `held-suar:94` Held & Suarez 1994 — Overview
- L55 `held-suar:94` Held & Suarez 1994 — Forcing
- L59 `held-suar:94` Held & Suarez 1994 — Forcing
- L123 `adcroft:04b` Adcroft et al. 2004 — Set-up description
- L179 `marshall:97a` Marshall et al. 1997 — Set-up description
- L190 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L199 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L209 `adcroft:95` Adcroft 1995 — Numerical Stability Criteria
- L560 `shapiro:70` Shapiro 1970 — File :filelink:`input/data.pkg <verification/tutorial_held_suarez_cs/input/data.pkg>`
- L589 `shapiro:70` Shapiro 1970 — File :filelink:`input/data.shap <verification/tutorial_held_suarez_cs/input/data.shap>`

### `doc/examples/reentrant_channel/reentrant_channel.rst`
- L19 `sverdrup:33` Sverdrup 1933 — Southern Ocean Reentrant Channel Example
- L22 `marshall:03` Marshall & Radko 2003 — Southern Ocean Reentrant Channel Example
- L23 `marshall:12` Marshall & Speer 2012 — Southern Ocean Reentrant Channel Example
- L23 `olbers:04` Olbers & Visbeck 2004 — Southern Ocean Reentrant Channel Example
- L24 `nikurashin:12` Nikurashin & Vallis 2012 — Southern Ocean Reentrant Channel Example
- L25 `armour:16` Armour et al. 2016 — Southern Ocean Reentrant Channel Example
- L25 `sallee:18` Sallée 2018 — Southern Ocean Reentrant Channel Example
- L28 `abernathy:11` Abernathey et al. 2011 — Southern Ocean Reentrant Channel Example
- L40 `gen-mcw:90` Gent & McWilliams 1990 — Southern Ocean Reentrant Channel Example
- L41 `danabasoglu:94` Danabasoglu et al. 1994 — Southern Ocean Reentrant Channel Example
- L43 `gent:11` Gent 2011 — Southern Ocean Reentrant Channel Example
- L144 `stewart:17` Stewart et al. 2017 — Discrete Numerical Configuration
- L247 `roach:15` Roach et al. 2015 — Numerical Stability Criteria
- L292 `gen-mcw:90` Gent & McWilliams 1990 — File :filelink:`code/packages.conf <verification/tutorial_reentrant_channel/code/packages.conf>`
- L395 `jamart:86` Jamart & Ozer 1986 — PARM01 - Continuous equation parameters
- L599 `danabasoglu:95` Danabasoglu & J.C. McWilliams 1995 — File :filelink:`input/data.gmredi <verification/tutorial_reentrant_channel/input/data.gmredi>`
- L611 `gr:98` Griffies 1998 — File :filelink:`input/data.gmredi <verification/tutorial_reentrant_channel/input/data.gmredi>`
- L776 `adcroft:97` Adcroft et al. 1997 — File ``input/bathy.50km.bin``
- L958 `gent:11` Gent 2011 — Coarse Resolution Solution
- L1008 `bryan:91` Burridge & Haseler 1977 — Coarse Resolution Solution
- L1008 `deacon:37` Deacon 1937 — Coarse Resolution Solution
- L1009 `doos:94` Döös & Webb 1994 — Coarse Resolution Solution
- L1009 `speer:00` Speer & Sloyan 2000 — Coarse Resolution Solution
- L1030 `ferrari:03` Ferrari & Pumb 2003 — Coarse Resolution Solution
- L1031 `wolfe:14` Wolfe 2014 — Coarse Resolution Solution
- L1047 `danabasoglu:94` Danabasoglu et al. 1994 — Coarse Resolution Solution
- L1060 `abernathy:11` Abernathey et al. 2011 — Coarse Resolution Solution
- L1068 `veronis:75` Veronis 1975 — Coarse Resolution Solution
- L1163 `leith:68` Leith 1968 — Eddy Permitting Solution
- L1163 `leith:96` Leith 1996 — Eddy Permitting Solution
- L1234 `dufour:12` Dufour et al. 2012 — Eddy Permitting Solution
- L1234 `viebahn:12` Viebahn & Eden 2012 — Eddy Permitting Solution
- L1253 `ferriera:05` Ferreira et al. 2005 — Eddy Permitting Solution

### `doc/examples/tracer_adjsens/tracer_adjsens.rst`
- L21 `hill:04` Hill et al. 2004 — Overview of the experiment
- L41 `gen-eta:95` Gent et al. 1995 — Passive tracer equation
- L41 `gen-mcw:90` Gent & McWilliams 1990 — Passive tracer equation
- L340 `giering:99` Giering 1999 — File ``makefile``

### `doc/getting_started/getting_started.rst`
- L1492 `bryan:79` Bryan & Lewis 1979 — C Preprocessor Options — `ALLOW_BL79_LAT_VARY`
- L2180 `bryan:72` Bryan & Cox 1972 — Parameters: Equation of State
- L2190 `fofonoff:83` Fofonoff & R. Millard 1983 — Parameters: Equation of State
- L2196 `jackett:95` Jackett & McDougall 1995 — Parameters: Equation of State
- L2203 `jackett:95` Jackett & McDougall 1995 — Parameters: Equation of State
- L2211 `mcdougall:03` McDougall et al. 2003 — Parameters: Equation of State
- L2220 `ioc:10` IOC et al. 2010 — Parameters: Equation of State
- L2221 `mcdougall:11` McDougall & Barker 2011 — Parameters: Equation of State
- L2222 `roquet:15` Roquet et al. 2015 — Parameters: Equation of State
- L2240 `millero:10` Millero 2010 — Parameters: Equation of State
- L2240 `pawlowicz:13` Pawlowicz 2013 — Parameters: Equation of State
- L2385 `burridge:77` (not in bib) — Configuration
- L2385 `sadourny:75` Sadourny 1975 — Configuration
- L2723 `bryan:79` Bryan & Lewis 1979 — Tracer Diffusivities — `diffKrBL79surf`
- L2725 `bryan:79` Bryan & Lewis 1979 — Tracer Diffusivities — `diffKrBL79deep`
- L2727 `bryan:79` Bryan & Lewis 1979 — Tracer Diffusivities — `diffKrBL79scl`
- L2729 `bryan:79` Bryan & Lewis 1979 — Tracer Diffusivities — `diffKrBL79Ho`

### `doc/ocean_state_est/ocean_state_est.rst`
- L24 `for-eta:15` Forget et al. 2015 — ECCO: model-data comparisons using gridded data sets
- L442 `fuku-etal:14` Fukumori et al. 2015 — Generic Integral Function
- L442 `heim-eta:11` Heimbach et al. 2011 — Generic Integral Function
- L442 `maro-eta:99` Marotzke et al. 1999 — Generic Integral Function
- L541 `smith:19` Smith & Heimbach 2019 — Custom Cost Functions
- L1324 `weaver:01` Weaver & Courtier 2001 — Generic Control Processing Options
- L1617 `gil-lem:89` Gilbert & Lemarechal 1989 — General features
- L1983 `gil-lem:89` Gilbert & Lemarechal 1989 — Alternative code to :filelink:`optim` and :filelink:`lsopt`
- L2078 `for-eta:15` Forget et al. 2015 — Test Cases For Estimation Package Capabilities

### `doc/overview/finding_pressure.rst`
- L93 `harlow:65` Harlow & Welch 1965 — Non-hydrostatic pressure
- L93 `potter:73` Potter 1973 — Non-hydrostatic pressure
- L93 `williams:69` Williams 1969 — Non-hydrostatic pressure
- L125 `williams:69` Williams 1969 — Boundary Conditions
- L146 `marshall:97a` Marshall et al. 1997 — Boundary Conditions
- L146 `marshall:97b` Marshall et al. 1997 — Boundary Conditions

### `doc/overview/global_atmos_hs.rst`
- L16 `held-suar:94` Held & Suarez 1994 — Global atmosphere: ‘Held-Suarez’ benchmark
- L29 `adcroft:04b` Adcroft et al. 2004 — Global atmosphere: ‘Held-Suarez’ benchmark

### `doc/overview/hydrostatic.rst`
- L35 `marshall:97a` Marshall et al. 1997 — Hydrostatic, Quasi-hydrostatic, Quasi-nonhydrostatic and Non-hydrostatic forms
- L95 `marshall:97a` Marshall et al. 1997 — Hydrostatic and quasi-hydrostatic forms
- L123 `marshall:97a` Marshall et al. 1997 — Hydrostatic and quasi-hydrostatic forms
- L143 `marshall:97a` Marshall et al. 1997 — Non-hydrostatic Ocean
- L143 `white:95` (not in bib) — Non-hydrostatic Ocean

### `doc/overview/overview.rst`
- L64 `hill:95` Hill & Marshall 1995 — Introduction
- L68 `marshall:97a` Marshall et al. 1997 — Introduction
- L72 `marshall:97b` Marshall et al. 1997 — Introduction
- L76 `adcroft:97` Adcroft et al. 1997 — Introduction
- L80 `mars-eta:98` Marshall et al. 1998 — Introduction
- L85 `adcroft:99` Adcroft et al. 1999 — Introduction
- L91 `hill:99` Hill et al. 1999 — Introduction
- L96 `maro-eta:99` Marotzke et al. 1999 — Introduction
- L100 `adcroft:04a` Adcroft & Campin 2004 — Introduction
- L105 `adcroft:04b` Adcroft et al. 2004 — Introduction
- L108 `marshall:04` Marshall et al. 2004 — Introduction
- L111 `adcroft:04c` Adcroft et al. 2004 — Introduction

### `doc/overview/soln_strategy.rst`
- L22 `marshall:97a` Marshall et al. 1997 — Solution strategy

### `doc/phys_pkgs/fizhi.rst`
- L32 `moorsz:92` Moorthi & Suarez 1992 — Sub-grid and Large-scale Convection
- L60 `moorsz:92` Moorthi & Suarez 1992 — Sub-grid and Large-scale Convection
- L122 `sudm:88` Sud & Molod 1988 — Sub-grid and Large-scale Convection
- L243 `rosen:87` Rosenfield et al. 1987 — Cloud Formation
- L253 `chou:90` Chou 1990 — Shortwave Radiation
- L253 `chou:92` Chou 1992 — Shortwave Radiation
- L257 `lhans:74` Lacis & Hansen 1974 — Shortwave Radiation
- L343 `chsz:94` Chou & Suarez 1994 — Longwave Radiation
- L498 `yam:77` (not in bib) — Turbulence
- L537 `helflab:88` Helfand & Labraga 1988 — Turbulence
- L619 `helfschu:95` Helfand & Schubert 1995 — Turbulence
- L624 `yagkad:74` Yamada 1977 — Turbulence
- L642 `kondo:75` Kondo 1975 — Turbulence
- L642 `larpond:81` Large & Pond 1981 — Turbulence
- L644 `dorsell:89` Dorman & Sellers 1989 — Turbulence
- L648 `pano:73` Panofsky 1973 — Turbulence
- L660 `clarke:70` Clarke 1970 — Turbulence
- L752 `ks:91` Koster & Suarez 1991 — Surface Type
- L752 `ks:92` Koster & Suarez 1992 — Surface Type
- L756 `deftow:94` Defries & Townshend 1994 — Surface Type
- L757 `dorsell:89` Dorman & Sellers 1989 — Surface Type
- L813 `helfschu:95` Helfand & Schubert 1995 — Surface Roughness
- L814 `kondo:75` Kondo 1975 — Surface Roughness
- L814 `larpond:81` Large & Pond 1981 — Surface Roughness
- L834 `zhouetal:95` Zhou et al. 1995 — Gravity Wave Drag
- L854 `taksz:96` Takacs & Suarez 1996 — Gravity Wave Drag
- L1405 `helflab:88` Helfand & Labraga 1988 — ET - Diffusivity Coefficient for Temperature and Moisture (m^2/sec)
- L1426 `helflab:88` Helfand & Labraga 1988 — ET - Diffusivity Coefficient for Temperature and Moisture (m^2/sec)
- L1446 `helflab:88` Helfand & Labraga 1988 — EU - Diffusivity Coefficient for Momentum (m^2/sec)
- L1468 `helflab:88` Helfand & Labraga 1988 — EU - Diffusivity Coefficient for Momentum (m^2/sec)

### `doc/phys_pkgs/generic_advdiff.rst`
- L54 `prather:86` Prather 1986 — Key subroutines, parameters and files — `GAD_ALLOW_TS_SOM_ADV`
- L57 `smolark:89` Smolarkiewicz 1989 — Key subroutines, parameters and files — `GAD_SMOLARKIEWICZ_HACK`

### `doc/phys_pkgs/ggl90.rst`
- L14 `gas-eta:90` Gaspar et al. 1990 — Key subroutines, parameters and files

### `doc/phys_pkgs/gmredi.rst`
- L14 `redi1982` Redi 1982 — Introduction
- L17 `gen-eta:95` Gent et al. 1995 — Introduction
- L17 `gen-mcw:90` Gent & McWilliams 1990 — Introduction
- L34 `cox87` Cox 1987 — Description
- L48 `gretal:98` Griffies et al. 1998 — Description
- L155 `danabasoglu:95` Danabasoglu & J.C. McWilliams 1995 — GM parameterization
- L189 `gr:98` Griffies 1998 — Griffies Skew Flux
- L350 `visbeck:97` Visbeck et al. 1997 — Visbeck et al. 1997 GM diffusivity :math:`\kappa_{GM}(x,y)`
- L379 `marshall:12b` Marshall et al. 2012 — Marshall et al. 2012 GM diffusivity :math:`\kappa_{GM}(x,y)`
- L387 `mak:18` Mak et al. 2018 — Marshall et al. 2012 GM diffusivity :math:`\kappa_{GM}(x,y)`
- L387 `mak:22` Mak et al. 2022 — Marshall et al. 2012 GM diffusivity :math:`\kappa_{GM}(x,y)`
- L397 `ferriera:05` Ferreira et al. 2005 — Marshall et al. 2012 GM diffusivity :math:`\kappa_{GM}(x,y)`
- L432 `lar-eta:94` Large et al. 1994 — Tapering and stability
- L443 `cox87` Cox 1987 — Slope clipping
- L504 `gkw:91` Gerdes et al. 1991 — Tapering: Gerdes, Koberle and Willebrand, 1991 (GKW91)
- L545 `danabasoglu:95` Danabasoglu & J.C. McWilliams 1995 — Tapering: Danabasoglu and McWilliams, 1995 (DM95)
- L562 `lar-eta:97` Large et al. 1997 — Tapering: Large, Danabasoglu and Doney, 1997 (LDD97)

### `doc/phys_pkgs/gridalt.rst`
- L10 `mol:09` Molod 2009 — Introduction

### `doc/phys_pkgs/kl10.rst`
- L16 `klymaklegg10` Klymak & Legg 2010 — Introduction
- L48 `thorpe77` Thorpe 1977 — Introduction
- L53 `moum96` Moum 1996 — Introduction
- L53 `seimgregg94` Seim & Gregg 1994 — Introduction
- L53 `wesson94` Wesson & Gregg 1994 — Introduction
- L75 `klymaklegg10` Klymak & Legg 2010 — Introduction

### `doc/phys_pkgs/kpp.rst`
- L16 `lar-eta:94` Large et al. 1994 — Introduction
- L50 `lar-eta:97` Large et al. 1997 — Introduction
- L257 `lar-eta:94` Large et al. 1994 — Equations and key routines
- L295 `lar-eta:94` Large et al. 1994 — BLMIX: Mixing in the boundary layer

### `doc/phys_pkgs/obcs.rst`
- L348 `orl:76` Orlanski 1976 — Equations and key routines
- L408 `stevens:90` Stevens 1990 — Equations and key routines
- L468 `stevens:90` Stevens 1990 — Equations and key routines

### `doc/phys_pkgs/opps.rst`
- L15 `pal-rom:97` Paluszkiewicz & Romea 1997 — Key subroutines, parameters and files

### `doc/phys_pkgs/rbcs.rst`
- L203 `adcroft:97` Adcroft et al. 1997 — Experiments and tutorials that use rbcs

### `doc/phys_pkgs/remesh.rst`
- L80 `jordan:18` Jordan et al. 2018 — Description
- L130 `losch:08` Losch 2008 — Alternate boundary layer formulation

### `doc/phys_pkgs/seaice.rst`
- L177 `lemieux:12` Lemieux et al. 2012 — General flags and parameters — `SEAICEuseEVPstar`
- L179 `bouillon:13` Bouillon et al. 2013 — General flags and parameters — `SEAICEuseEVPrev`
- L191 `kimmritz:16` Kimmritz et al. 2016 — General flags and parameters — `SEAICEaEVPcStar`
- L194 `kimmritz:16` Kimmritz et al. 2016 — General flags and parameters — `SEAICEaEVPalphaMin`
- L339 `lemieux:15` Lemieux et al. 2015 — General flags and parameters — `SEAICEbasalDragK1`
- L347 `liu:22` Liu et al. 2022 — General flags and parameters — `SEAICESideDrag`
- L368 `hibler:79` Hibler 1979 — General flags and parameters — `useHibler79IceStrength`
- L368 `rothrock:75` Rothrock 1975 — General flags and parameters — `useHibler79IceStrength`
- L371 `hibler:79` Hibler 1979 — General flags and parameters — `SEAICEsimpleRidging`
- L373 `rothrock:75` Rothrock 1975 — General flags and parameters — `SEAICE_cf`
- L375 `thorndike:75` Thorndike et al. 1975 — General flags and parameters — `SEAICEpartFunc`
- L377 `hibler:80` Hibler 1980 — General flags and parameters — `SEAICEredistFunc`
- L383 `thorndike:75` Thorndike et al. 1975 — General flags and parameters — `SEAICEgStar`
- L385 `lipscomb:07` Lipscomb et al. 2007 — General flags and parameters — `SEAICEhStar`
- L385 `thorndike:75` Thorndike et al. 1975 — General flags and parameters — `SEAICEhStar`
- L388 `lipscomb:07` Lipscomb et al. 2007 — General flags and parameters
- L391 `lipscomb:07` Lipscomb et al. 2007 — General flags and parameters
- L397 `lipscomb:01` Lipscomb 2001 — General flags and parameters — `SEAICEuseLinRemapITD`
- L403 `lipscomb:01` Lipscomb 2001 — General flags and parameters — `Hlimit_c3`
- L413 `zhang:97` Zhang & Hibler 1997 — Description
- L415 `hibler:79` Hibler 1979 — Description
- L415 `hibler:80` Hibler 1980 — Description
- L417 `losch:10` Losch et al. 2010 — Description
- L424 `zhang:97` Zhang & Hibler 1997 — Description
- L425 `hunke:97` Hunke & Dukowicz 1997 — Description
- L426 `lemieux:10` Lemieux et al. 2010 — Description
- L427 `losch:14` Losch et al. 2014 — Description
- L430 `campin:08` Campin et al. 2008 — Description
- L430 `hibler:87` Hibler & Bryan 1987 — Description
- L444 `hibler:79` Hibler 1979 — Description
- L445 `flato:92` Flato & W. D. Hibler 1992 — Description
- L446 `hunke:97` Hunke & Dukowicz 1997 — Description
- L454 `zhang:97` Zhang & Hibler 1997 — Description
- L459 `semtner:76` Semtner 1976 — Description
- L463 `fenty:13` Fenty & Heimbach 2013 — Description
- L480 `winton:00` Winton 2000 — Compatibility with ice-thermodynamics package :filelink:`pkg/thsice`
- L626 `hibler:79` Hibler 1979 — Viscous-Plastic (VP) Rheology — `SEAICE_cStar`
- L651 `konig:10` König Beatty & Holland 2010 — Viscous-Plastic (VP) Rheology — `SEAICE_tensilFac`
- L711 `hibler:79` Hibler 1979 — Elliptical yield curve with normal flow rule
- L788 `ringeisen:20` Ringeisen et al. 2020 — Elliptical yield curve with non-normal flow rule
- L832 `hibler:00` Hibler & Schulson 2000 — Truncated ellipse method (TEM) for elliptical yield curve
- L833 `ringeisen:19` Ringeisen et al. 2019 — Truncated ellipse method (TEM) for elliptical yield curve
- L863 `ip:91` Ip et al. 1991 — Mohr-Coulomb yield curve with shear flow rule — `SEAICE_ALLOW_MCS`
- L882 `zha:05` Zhang & Rothrock 2005 — Teardrop yield curve with normal flow rule
- L901 `zha:05` Zhang & Rothrock 2005 — Parabolic lens yield curve with normal flow rule
- L940 `zhang:97` Zhang & Hibler 1997 — LSR and JFNK solver
- L943 `lemieux:10` Lemieux et al. 2010 — LSR and JFNK solver
- L945 `losch:14` Losch et al. 2014 — LSR and JFNK solver
- L1031 `losch:14` Losch et al. 2014 — LSR and JFNK solver
- L1073 `hutchings:04` Hutchings et al. 2004 — LSR and JFNK solver
- L1085 `hunke:97` Hunke & Dukowicz 1997 — Elastic-Viscous-Plastic (EVP) Dynamics
- L1099 `hunke:97` Hunke & Dukowicz 1997 — Elastic-Viscous-Plastic (EVP) Dynamics
- L1137 `hunke:97` Hunke & Dukowicz 1997 — Elastic-Viscous-Plastic (EVP) Dynamics
- L1149 `hunke:97` Hunke & Dukowicz 1997 — Elastic-Viscous-Plastic (EVP) Dynamics
- L1166 `bouillon:13` Bouillon et al. 2013 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1166 `hunke:01` Hunke 2001 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1166 `lemieux:12` Lemieux et al. 2012 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1168 `bouillon:13` Bouillon et al. 2013 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1168 `kimmritz:15` Kimmritz et al. 2015 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1168 `lemieux:12` Lemieux et al. 2012 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1187 `kimmritz:15` Kimmritz et al. 2015 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1213 `kimmritz:15` Kimmritz et al. 2015 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1225 `bouillon:13` Bouillon et al. 2013 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1228 `kimmritz:16` Kimmritz et al. 2016 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1240 `kimmritz:16` Kimmritz et al. 2016 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1245 `kimmritz:16` Kimmritz et al. 2016 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1258 `kimmritz:15` Kimmritz et al. 2015 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1258 `kimmritz:16` Kimmritz et al. 2016 — More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP
- L1270 `hibler:87` Hibler & Bryan 1987 — Ice-Ocean stress
- L1276 `hibler:87` Hibler & Bryan 1987 — Ice-Ocean stress
- L1281 `hibler:87` Hibler & Bryan 1987 — Ice-Ocean stress
- L1581 `semtner:76` Semtner 1976 — Zero-layer thermodynamics
- L1593 `parkinson:79` Parkinson & Washington 1979 — Zero-layer thermodynamics
- L1594 `manabe:79` Manabe et al. 1979 — Zero-layer thermodynamics
- L1647 `semtner:76` Semtner 1976 — Zero-layer thermodynamics
- L1656 `hibler:84` Hibler 1984 — Zero-layer thermodynamics
- L1665 `castro-morales:14` Castro-Morales et al. 2014 — Zero-layer thermodynamics
- L1669 `castro-morales:14` Castro-Morales et al. 2014 — Zero-layer thermodynamics
- L1679 `hibler:79` Hibler 1979 — Zero-layer thermodynamics
- L1679 `hibler:80` Hibler 1980 — Zero-layer thermodynamics
- L1684 `zha:98` Zhang et al. 1998 — Zero-layer thermodynamics
- L1694 `leppaeranta:83` Lepparanta 1983 — Zero-layer thermodynamics
- L1725 `winton:00` Winton 2000 — Advection of thermodynamic variables
- L1726 `semtner:76` Semtner 1976 — Advection of thermodynamic variables
- L1734 `roe:85` Roe 1985 — Advection of thermodynamic variables
- L1735 `winton:00` Winton 2000 — Advection of thermodynamic variables
- L1757 `thorndike:75` Thorndike et al. 1975 — Dynamical Ice Thickness Distribution (ITD)
- L1758 `rothrock:75` Rothrock 1975 — Dynamical Ice Thickness Distribution (ITD)
- L1760 `ungermann:17` Ungermann et al. 2017 — Dynamical Ice Thickness Distribution (ITD)
- L1769 `thorndike:75` Thorndike et al. 1975 — Distribution, participation and redistribution functions in ridging
- L1787 `hibler:79` Hibler 1979 — Distribution, participation and redistribution functions in ridging
- L1802 `thorndike:75` Thorndike et al. 1975 — Distribution, participation and redistribution functions in ridging
- L1803 `lipscomb:07` Lipscomb et al. 2007 — Distribution, participation and redistribution functions in ridging — `SEAICEpartFunc`
- L1805 `hibler:80` Hibler 1980 — Distribution, participation and redistribution functions in ridging
- L1806 `lipscomb:07` Lipscomb et al. 2007 — Distribution, participation and redistribution functions in ridging — `SEAICEredistFunc`
- L1808 `lipscomb:07` Lipscomb et al. 2007 — Distribution, participation and redistribution functions in ridging
- L1811 `lipscomb:07` Lipscomb et al. 2007 — Distribution, participation and redistribution functions in ridging
- L1850 `bitz:01` Bitz et al. 2001 — Distribution, participation and redistribution functions in ridging
- L1861 `lipscomb:01` Lipscomb 2001 — Distribution, participation and redistribution functions in ridging
- L1878 `hibler:79` Hibler 1979 — Ice strength parameterization
- L1883 `rothrock:75` Rothrock 1975 — Ice strength parameterization
- L1904 `ungermann:17` Ungermann et al. 2017 — Ice strength parameterization

### `doc/phys_pkgs/shelfice.rst`
- L54 `holland:99` Holland & Jenkins 1999 — SHELFICE configuration — `SHI_ALLOW_GAMMAFRICT`
- L142 `holland:99` Holland & Jenkins 1999 — SHELFICE run-time parameters — `SHELFICEuseGammaFrict`
- L158 `marshall:04` Marshall et al. 2004 — SHELFICE description
- L183 `beckmann:99` Beckmann et al. 1999 — SHELFICE description
- L311 `jenkins:01` Jenkins et al. 2001 — Three-equations thermodynamics
- L328 `holland:99` Holland & Jenkins 1999 — Three-equations thermodynamics
- L329 `jenkins:01` Jenkins et al. 2001 — Three-equations thermodynamics
- L376 `holland:99` Holland & Jenkins 1999 — Three-equations thermodynamics
- L405 `holland:99` Holland & Jenkins 1999 — Three-equations thermodynamics
- L424 `hellmer:89` Hellmer & Olbers 1989 — Three-equations thermodynamics
- L425 `jenkins:01` Jenkins et al. 2001 — Three-equations thermodynamics
- L535 `holland:99` Holland & Jenkins 1999 — Solving the three-equations system
- L601 `grosfeld:97` Grosfeld et al. 1997 — ISOMIP thermodynamics
- L629 `holland:99` Holland & Jenkins 1999 — Exchange coefficients
- L638 `losch:08` Losch 2008 — Remark

### `doc/phys_pkgs/streamice.rst`
- L355 `Macayeal:89` MacAyeal 1989 — Equations Solved
- L460 `asay-davis:16` Asay-Davis et al. 2016 — Equations Solved
- L499 `goldberg:2011` Goldberg 2011 — Hybrid SIA-SSA stress balance
- L515 `goldberg:2011` Goldberg 2011 — Hybrid SIA-SSA stress balance
- L525 `Albrecht:2011` Albrecht et al. 2011 — Ice front advance
- L579 `goldberg:2011` Goldberg 2011 — Numerical Details
- L663 `Goldberg:2015` Goldberg et al. 2015 — Boundary Stresses
- L697 `goldberg_heimbach:2013` Goldberg & Heimbach 2013 — Adjoint
- L698 `goldberg_openad_fixed:2016` Goldberg et al. 2016 — Adjoint
- L699 `christianson:94` Christianson 1994 — Adjoint

### `doc/software_arch/software_arch.rst`
- L125 `hoe:99` Hoe et al. 1999 — Target hardware
- L383 `hoe:99` Hoe et al. 1999 — Distributed memory communication

## Part 2: full references (`key`: citation; cited N times)

- `abernathy:11`: Abernathey et al. (2011). The Dependence of Southern Ocean Meridional Overturning on Wind Stress. J. Phys. Oceanogr. 41, 2261–2278. doi:10.1175/JPO-D-11-023.1 (2×)
- `adcroft:04a`: Adcroft & Campin (2004). Re-scaled height coordinates for accurate representation of free-surface flows in ocean circulation models. Ocean Modelling 7, 269-284. doi:10.1016/j.ocemod.2003.09.003 (3×)
- `adcroft:04b`: Adcroft et al. (2004). Implementation of an Atmosphere-Ocean General Circulation Model on the Expanded Spherical Cube. Mon. Wea. Rev. 132, 2845-2863. doi:10.1175/MWR2823.1 (4×)
- `adcroft:04c`: Adcroft et al. (2004). Overview of the Formulation and Numerics of the MITGCM. Proceedings of the ECMWF seminar series on Numerical Methods, Recent developments in numerical methods for atmosphere and ocean modelling, 139-149. http://mitgcm.org/pdfs/ECMWF2004-Adcroft.pdf (1×)
- `adcroft:95`: Adcroft (1995). Numerical Algorithms for use in a Dynamical Model of the Ocean. Imperial College, London. https://extranet.gfdl.noaa.gov/~aja/papers/adcroft_PhD_1995.pdf (13×)
- `adcroft:97`: Adcroft et al. (1997). Representation of topography by shaved cells in a height coordinate ocean model. Mon. Wea. Rev. 125, 2293-2315. doi:10.1175/1520-0493\%281997\%29125<2293:ROTBSC>2.0.CO;2 (6×)
- `adcroft:98`: Adcroft & Marshall (1998). How slippery are piecewise-constant coastlines in numerical ocean models?. Tellus 50, 95-108 (1×)
- `adcroft:99`: Adcroft et al. (1999). A new treatment of the Coriolis terms in C-grid models at both high and low resolutions. Mon. Wea. Rev. 127, 1928-1936. doi:10.1175/1520-0493\%281999\%29127<1928:ANTOTC>2.0.CO;2 (1×)
- `Albrecht:2011`: Albrecht et al. (2011). Parameterization for subgrid-scale motion of ice-shelf calving fronts. The Cryosphere 5, 35–44. doi:10.5194/tc-5-35-2011 (1×)
- `arakawa:77`: Arakawa & Lamb (1977). Computational design of the basic dynamical processes of the UCLA general circulation model. Meth. Comput. Phys. 17, 174-267 (1×)
- `armour:16`: Armour et al. (2016). Southern Ocean warming delayed by circumpolar upwelling and equatorward transport. Nature Geosci. 9, 549–554. doi:10.1038/ngeo2731 (1×)
- `asay-davis:16`: Asay-Davis et al. (2016). Experimental design for three interrelated marine ice sheet and ocean model intercomparison projects: MISMIP v. 3 (MISMIP +), ISOMIP v. 2 (ISOMIP +) and MISOMIP v. 1 (MISOMIP1). Geosci. model dev. 9, 2471–2497. doi:10.3929/ethz-b-000119139 (1×)
- `bachman:17`: Bachman et al. (2017). A scale-aware subgrid model for quasi-geostrophic turbulence. J. Geophys. Res. Ocean. 122, 1529–1554. doi:10.1002/2016JC012265 (3×)
- `beckmann:99`: Beckmann et al. (1999). A numerical model of the Weddell Sea: Large-scale circulation and water mass distribution. J. Geophys. Res. Oceans 104, 23375–23391. doi:10.1029/1999JC900194 (1×)
- `bitz:01`: Bitz et al. (2001). Simulating the ice-thickness distribution in a coupled climate model. J. Geophys. Res. 106, 2441. doi:10.1029/1999JC000113 (1×)
- `bouillon:13`: Bouillon et al. (2013). The Elastic-Viscous-Plastic Method Revisited. Ocean Modelling 71, 2–12. doi:10.1016/j.ocemod.2013.05.013 (4×)
- `bryan:63`: Bryan (1963). A numerical investigation of a nonlinear model of a wind-driven ocean. J. Atmos. Sci. 20, 594-606 (1×)
- `bryan:72`: Bryan & Cox (1972). An approximate equation of state for numerical models of ocean circulation. J. Phys. Oceanogr. 2, 510–514. doi:10.1175/1520-0485(1972)002<0510:AAEOSF>2.0.CO;2 (1×)
- `bryan:75`: Bryan et al. (1975). A global ocean-atmosphere climate model. Part II. The oceanic circulation. J. Phys. Oceanogr. 5, 30–46 (2×)
- `bryan:79`: Bryan & Lewis (1979). A water mass model of the world ocean. J. Geophys. Res. 84, 2503–2517. doi:10.1029/JC084iC05p02503 (7×)
- `bryan:84`: Bryan (1984). Accelerating the convergence to equilibrium of ocean-climate models. J. Phys. Oceanogr. 14, 666-673. doi:10.1175/1520-0485(1984)014<0666:ATCTEO>2.0.CO;2 (1×)
- `bryan:91`: Burridge & Haseler (1977). A Model for medium range weather forecasting: Adiabatic Formulation. Strategies for Future Climate Research, 196 pp.. https://www.ecmwf.int/sites/default/files/elibrary/1977/8495-model-medium-range-weather-forecasts-adiabatic-formulation.pdf (1×)
- `cam:04`: Campin et al. (2004). Conservation of properties in a free-surface model. Ocean Modelling 6, 221–244. doi:10.1016/s1463-5003(03)00009-x (4×)
- `campin:08`: Campin et al. (2008). Sea ice–ocean coupling using a rescaled vertical coordinate z*. Ocean Modelling 24, 1–14. doi:10.1016/j.ocemod.2008.05.005 (1×)
- `castro-morales:14`: Castro-Morales et al. (2014). Sensitivity of Simulated Arctic Sea Ice to Realistic Ice Thickness Distributions and Snow Parameterizations. J. Geophys. Res. Oceans 119, 559–571. doi:10.1002/2013JC009342 (2×)
- `chou:90`: Chou (1990). Parameterizations for the absorption of solar radiation by O_2 and CO_2 with applications to climate studies.. J. Clim. 3, 209-217 (1×)
- `chou:92`: Chou (1992). A solar radiation model for use in climate studies.. J. Atmos. Sci. 49, 762-772 (1×)
- `christianson:94`: Christianson (1994). Reverse accumulation and attractive fixed points. Optim. Method. Softw. 9, 307-322. doi:10.1080/10556789408805572 (1×)
- `chsz:94`: Chou & Suarez (1994). An efficient thermal infrared radiation parameterization for use in general circulation models. National Aeronautics and Space Administration (1×)
- `clarke:70`: Clarke (1970). Observational studies in the atmospheric boundary layer.. Q. J. R. Meteorol. Soc. 96, 91-114 (1×)
- `cox87`: Cox (1987). Isopycnal diffusion in a z-coordinate ocean model. Ocean modelling (unpublished manuscripts) 74, 1-5 (2×)
- `cox:84`: Cox & Bryan (1984). A Numerical Model of the Ventilated Thermocline. J. Phys. Oceanogr. 14, 674-687. doi:10.1175/1520-0485(1984)014<0674:ANMOTV>2.0.CO;2 (1×)
- `culkin:65`: Culkin (1965). The major constituents of seawater. Chemical Oceanography 1, 121–161 (1×)
- `cushmanroisin:11`: Cushman-Roisin & Beckers (2011). Introduction to Geophysical Fluid Dynamics, 2nd Edition. Academic Press, 875 pp. (2×)
- `danabasoglu:94`: Danabasoglu et al. (1994). The Role of Mesoscale Tracer Transports in the Global Ocean Circulation. Science 264, 1123–1126. doi:10.1126/science.264.5162.1123 (2×)
- `danabasoglu:95`: Danabasoglu & J.C. McWilliams (1995). Sensitivity of the Global Ocean Circulation to Parameterizations of Mesoscale Tracer Transports. J. Clim. 8, 2967–2987. doi:10.1175/1520-0442(1995)008<2967:SOTGOC>2.0.CO;2 (3×)
- `deacon:37`: Deacon (1937). The hydrology of the southern ocean. Discovery Rept. 15, 1-124 (1×)
- `deftow:94`: Defries & Townshend (1994). NDVI-derived Land Cover Classification at Global Scales.. Int'l J. Rem. Sens. 15, 3567-3586 (1×)
- `deszoeke:02`: de Szoeke & Samelson (2002). The duality between the Boussinesq and Non-Boussinesq hydrostatic equations of motion. J. Phys. Oceanogr. 32, 2194-2203. doi:10.1175/1520-0485(2002)032<2194:TDBTBA>2.0.CO;2 (1×)
- `dickson:79`: Dickson & Riley (1979). The estimation of acid dissociation constants in seawater media from potentionmetric titrations with strong base. I. The ionic product of water — Kw. Marine Chemistry 7, 89–99. doi:10.1016/0304-4203(79)90001-X (1×)
- `doos:94`: Döös & Webb (1994). The Deacon Cell and the Other Meridional Cells of the Southern Ocean. J. Phys. Oceanogr. 24, 429–442. doi:10.1175/1520-0485(1994)024<0429:TDCATO>2.0.CO;2 (1×)
- `dorsell:89`: Dorman & Sellers (1989). A global climatology of albedo, roughness length and stomatal resistance for atmospheric general circulation models as represented by the Simple Biosphere model (SiB).. J. Appl. Meteor. 28, 833-855 (2×)
- `dufour:12`: Dufour et al. (2012). Standing and Transient Eddies in the Response of the Southern Ocean Meridional Overturning to the Southern Annular Mode. J. Clim. 25, 6958 - 6974. doi:10.1175/JCLI-D-11-00309.1 (1×)
- `durran:91`: Durrann (1991). The Third-Order Adams-Bashforth Method: An Attractive Alternative to Leapfrog Time Differencing. Mon. Wea. Rev. 119, 702–720. doi:10.1175/1520-0493(1991)119<0702:TTOABM>2.0.CO;2 (1×)
- `dutay:02`: Dutkiewicz et al. (2005). A three-dimensional ocean-seaice-carbon cycle model and its coupling to a two-dimensional atmospheric model: Uses in climate change studies. Ocean Modelling 4, 47 pp.. doi:10.1016/S1463-5003(01)00013-0 (1×)
- `dutkiewicz:05`: Dutkiewicz et al. (2005). A three-dimensional ocean-seaice-carbon cycle model and its coupling to a two-dimensional atmospheric model: Uses in climate change studies. MIT Joint Program of the Science and Policy of Global Change, 47 pp.. http://web.mit.edu/globalchange/www/MITJPSPGC_Rpt122.pdf (2×)
- `dutkiewicz:2005`: Dutkiewicz et al. (2005). Interactions of the iron and phosphorus cycles: A three-dimensional model study. Global Biogeochemical Cycles 19. doi:10.1029/2004GB002342 (1×)
- `dutkiewicz:2015`: Dutkiewicz et al. (2015). Capturing optically important constituents and properties in a marine biogeochemical and ecosystem model. Biogeosciences 12, 4447–4481. doi:10.5194/bg-12-4447-2015 (1×)
- `dutkiewicz:2020`: Dutkiewicz et al. (2020). Dimensions of marine phytoplankton diversity. Biogeosciences 17, 609–634. doi:10.5194/bg-17-609-2020 (1×)
- `fenty:13`: Fenty & Heimbach (2013). Coupled sea icetextendashocean-state estimate in the Labrador Sea and Baffin Bay. J. Phys. Oceanogr. 43, 884-904. doi:10.1175/JPO-D-12-065.1 (1×)
- `ferrari:03`: Ferrari & Pumb (2003). Residual circulation in the ocean. Proceedings of the 13th 'Aha Huliko'a Hawaiian Winter Workshop 13, 219-228. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.518.57&rep=rep1&type=pdf (1×)
- `ferrari:08`: Ferrari et al. (2008). Parameterization of Eddy Fluxes near Oceanic Boundaries. J. Clim. 21, 2770–2789. doi:10.1175/2007JCLI1510.1 (1×)
- `ferrari:10`: Ferrari et al. (2010). A boundary-value problem for the parameterized mesoscale eddy transport. Ocean Modelling 32, 143-156. doi:10.1016/j.ocemod.2010.01.004 (1×)
- `ferriera:05`: Ferreira et al. (2005). Estimating eddy stresses by fitting dynamics to observations using a residual-mean ocean circulation model and its adjoint. J. Phys. Oceanogr. 35, 1891-1910. doi:10.1175/JPO2785.1 (3×)
- `flato:92`: Flato & W. D. Hibler (1992). Modeling pack ice as a cavitating fluid. J. Phys. Oceanogr. 22, 626-651 (1×)
- `fofonoff:83`: Fofonoff & R. Millard (1983). Algorithms for computation of fundamental properties of seawater. UNESCO (1×)
- `follows:06`: Follows et al. (2006). On the solution of the carbonate chemistry system in ocean biogeochemistry models. Ocean Modelling 12, 290-301. doi:10.1016/j.ocemod.2005.05.004 (2×)
- `for-eta:15`: Forget et al. (2015). ECCO version 4: an integrated framework for non-linear inverse modeling and global ocean state estimation. Geoscientific Model Development 8, 3071–3104. doi:10.5194/gmd-8-3071-2015 (2×)
- `fuku-etal:14`: Fukumori et al. (2015). A near-uniform fluctuation of ocean bottom pressure and sea level across the deep ocean basins of the Arctic Ocean and the Nordic Seas. Progress in Oceanography 134, 152 - 172. doi:10.1016/j.pocean.2015.01.013 (1×)
- `gaikwad:24`: Gaikwad et al. (2025). MITgcm-AD v2: Open source tangent linear and adjoint modeling framework for the oceans and atmosphere enabled by the Automatic Differentiation tool Tapenade. Future Generation Computer Systems 163, 107512. doi:10.1016/j.future.2024.107512 (2×)
- `garcia:92`: Garcia & Gordon (1992). Oxygen solubility in seawater: Better fitting equations. Limnology and Oceanography 37, 1307–1312. doi:10.4319/lo.1992.37.6.1307 (1×)
- `gas-eta:90`: Gaspar et al. (1990). A Simple Eddy Kinetic Energy Model for Simulations of the Oceanic Vertical Mixing: Tests at Station Papa and Long-Term Upper Ocean Study Site. J. Geophys. Res. 95, 16,179–16,193. doi:10.1029/JC095iC09p16179 (4×)
- `gen-eta:95`: Gent et al. (1995). Parameterizing eddy-induced tracer transports in ocean circulation models. J. Phys. Oceanogr. 25, 463-474. doi:10.1175/1520-0485(1995)025<0463:PEITTI>2.0.CO;2 (2×)
- `gen-mcw:90`: Gent & McWilliams (1990). Isopycnal mixing in ocean circulation models. J. Phys. Oceanogr. 20, 150-155. doi:10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2 (5×)
- `gent:11`: Gent (2011). The Gent–McWilliams parameterization: 20/20 hindsight. Ocean Modelling 39, 2-9. doi:10.1016/j.ocemod.2010.08.002 (2×)
- `giering:00`: Giering (2000). Tangent linear and adjoint biogeochemical models. Inverse Methods in Global Biogeochemical Cycles, 33-48. doi:10.1029/GM114p0033 (1×)
- `giering:98`: Giering & Kaminski (1998). Recipes for adjoint code construction. ACM Transactions on Mathematical Software 24, 437-474. doi:10.1145/293686.293695 (1×)
- `giering:99`: Giering (1999). Tangent linear and adjoint model compiler. users manual 1.4 (tamc version 5.2). Massachusetts Institute of Technology. http:autodiff.com/tamc/tamc_manual.ps.gz (1×)
- `gil-lem:89`: Gilbert & Lemarechal (1989). Some numerical experiments with variable-storage quasi-Newton algorithms. Math. Programming 45, 407–435. doi:10.1007/BF01589113 (4×)
- `gill:82`: Gill (1982). Atmosphere-Ocean Dynamics. Academic Press, 662 pp. (2×)
- `gkw:91`: Gerdes et al. (1991). The influence of numerical advection schemes on the results of ocean general circulation models. Clim. Dynamics 5, 211-226. doi:10.1007/BF00210006 (1×)
- `goldberg:2011`: Goldberg (2011). A variationally-derived, depth-integrated approximation to a higher-order glaciologial flow model. J. of Glaciology 57, 157–170 (3×)
- `Goldberg:2015`: Goldberg et al. (2015). Committed retreat of Smith, Pope, and Kohler Glaciers over the next 30 years inferred by transient model calibration. The Cryosphere 9, 2429–2446 (1×)
- `goldberg_heimbach:2013`: Goldberg & Heimbach (2013). Parameter and state estimation with a time-dependent adjoint marine ice sheet model. The Cryosphere 7, 1659–1678 (1×)
- `goldberg_openad_fixed:2016`: Goldberg et al. (2016). An optimized treatment for algorithmic differentiation of an important glaciological fixed-point problem. Geoscientific Model Development 9, 1891–1904 (1×)
- `gr:98`: Griffies (1998). The Gent-McWilliams Skew Flux. J. Phys. Oceanogr. 28, 831-841. doi:10.1175/1520-0485(1998)028<0831:TGMSF>2.0.CO;2 (2×)
- `gretal:98`: Griffies et al. (1998). Isoneutral diffusion in a z-coordinate ocean model. J. Phys. Oceanogr. 28, 805-830. doi:10.1175/1520-0485(1998)028<0805:IDIAZC>2.0.CO;2 (1×)
- `griewank:08`: Griewank & Walther (2008). Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation, Second Edition. SIAM, 426 pp. (1×)
- `griewank:92`: Griewank (1992). Achieving logarithmic growth of temporal and spatial complexity in reverse Automatic Differentiation. Optimization Methods and Software 1, 35-54 (1×)
- `griffies:00`: Griffies & Hallberg (2000). Biharmonic friction with a Smagorinsky-like viscosity for use in large-scale eddy-permitting ocean models. Mon. Wea. Rev. 128, 2935-2946 (7×)
- `grosfeld:97`: Grosfeld et al. (1997). Thermohaline circulation and interaction between ice shelf cavities and the adjacent open water. J. Geophys. Res. Oceans 102, 15595-15610. doi:10.1029/97JC00891 (1×)
- `haney:71`: Haney (1971). Surface thermal boundary conditions for ocean circulation models. J. Phys. Oceanogr. 1, 241-248. doi:10.1175/1520-0485(1971)001<0241:STBCFO>2.0.CO;2 (2×)
- `harlow:65`: Harlow & Welch (1965). Numerical Calculation of Time-Dependent Viscous Incompressible Flow of Fluid with Free Surface. Physics of Fluids 8, 2182-2189 (1×)
- `hascoet:24`: Hascoet et al. (2024). Profiling checkpointing schedules in adjoint ST-AD.. doi:10.48550/arXiv.2405.15590 (1×)
- `heim-eta:11`: Heimbach et al. (2011). Timescales and regions of the sensitivity of Atlantic meridional volume and heat transport: Toward observing system design. Deep Sea Research Part II: Topical Studies in Oceanography 58, 1858–1879 (1×)
- `held-suar:94`: Held & Suarez (1994). A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models. Bulletin of the American Meteorological Society 75(10), 1825-1830 (11×)
- `helflab:88`: Helfand & Labraga (1988). Design of a non-singular level 2.5 second-order closure model for the prediction of atmospheric turbulence.. J. Atmos. Sci. 45, 113-132 (5×)
- `helfschu:95`: Helfand & Schubert (1995). Climatology of the Simulated Great Plains Low-Level Jet and Its contribution to the Continental Moisture Budget of the United States.. J. Clim. 8, 784-806 (2×)
- `hellmer:89`: Hellmer & Olbers (1989). A two-dimensional model of the thermohaline circulation under an ice shelf. Antarct. Sci. 1, 325-336. doi:10.1017/S0954102089000490 (3×)
- `hibler:00`: Hibler & Schulson (2000). On Modeling the Anisotropic Failure and Flow of Flawed Sea Ice. J. Geophys. Res. Oceans 105, 17105–17120. doi:10.1029/2000JC900045 (1×)
- `hibler:79`: Hibler (1979). A Dynamic Thermodynamic Sea Ice Model. J. Phys. Oceanogr. 9, 815-846 (9×)
- `hibler:80`: Hibler (1980). Modeling a Variable Thickness Sea Ice Cover. Mon. Wea. Rev. 1, 1943–1973 (4×)
- `hibler:84`: Hibler (1984). The role of sea ice dynamics in modeling CO2 increases. Climate processes and climate sensitivity 29, 238–253 (1×)
- `hibler:87`: Hibler & Bryan (1987). A Diagnostic Ice-Ocean Model. J. Phys. Oceanogr. 17, 987–1015 (5×)
- `hill:04`: Hill et al. (2004). Evaluating carbon sequestration efficiency in an ocean circulation model by adjoint sensitivity analysis. J. Geophys. Res. Oceans 109. doi:10.1029/2002JC001598 (1×)
- `hill:95`: Hill & Marshall (1995). Application of a Parallel Navier-Stokes Model to Ocean Circulation in Parallel Computational Fluid Dynamics. Implementations and Results Using Parallel Computers, 545-552 (1×)
- `hill:99`: Hill et al. (1999). A Strategy for Terascale Climate Modeling. In Proceedings of the Eighth ECMWF Workshop on the Use of Parallel Processors in Meteorology, 406-425 (1×)
- `hoe:99`: Hoe et al. (1999). A personal supercomputer for climate research. SC’99: Proceedings of the 1999 ACM/IEEE Conference on Supercomputing, 59. doi:10.1109/SC.1999.10009 (2×)
- `holland:78`: Holland (1978). The role of mesoscale eddies in the general circulation of the ocean-numerical experiments using a wind-driven quasi-geostrophic model. J. Phys. Oceanogr. 8, 363-392 (1×)
- `holland:99`: Holland & Jenkins (1999). Modeling Thermodynamic Ice–Ocean Interactions at the Base of an Ice Shelf. J. Phys. Oceanogr. 29, 1787-1800. doi:10.1175/1520-0485(1999)029<1787:MTIOIA>2.0.CO;2 (7×)
- `hunke:01`: Hunke (2001). Viscous-Plastic Sea Ice Dynamics with the EVP Model: Linearization Issues. J. Comput. Phys. 170, 18–38. doi:10.1006/jcph.2001.6710 (1×)
- `hunke:97`: Hunke & Dukowicz (1997). An elastic-viscous-plastic model for sea ice dynamics. J. Phys. Oceanogr. 27, 1849–1867. doi:10.1175/1520-0485(1997)027<1849:AEVPMF>2.0.CO;2 (6×)
- `hutchings:04`: Hutchings et al. (2004). A Strength Implicit Correction Scheme for the Viscous-Plastic Sea Ice Model. Ocean Modelling 7, 111-133. doi:10.1016/S1463-5003(03)00040-4 (1×)
- `ioc:10`: IOC et al. (2010). The International Thermodynamic Equation of Seawater 2010 (TEOS-10): Calculation and Use of Thermodynamic Properties. UNESCO (English), 196 pp.. http://www.teos-10.org/pubs/TEOS-10_Manual.pdf (1×)
- `ip:91`: Ip et al. (1991). On the Effect of Rheology on Seasonal Sea-Ice Simulations. Ann. Glaciol. 15, 17–25 (1×)
- `jackett:95`: Jackett & McDougall (1995). Minimal adjustment of hydrographic profiles to achieve static stability. J. Atmos. Ocean. Technol. 12, 381-389. doi:10.1175/1520-0426(1995)012<0381:MAOHPT>2.0.CO;2 (4×)
- `jamart:86`: Jamart & Ozer (1986). Numerical boundary layers and spurious residual flows. J. Geophys. Res. 91, 10621– 10631. doi:10.1029/JC091iC09p10621 (1×)
- `jenkins:01`: Jenkins et al. (2001). The role of meltwater advection in the formulation of conservative boundary conditions at an ice-ocean interface. J. Phys. Oceanogr. 31, 285-296. doi:10.1175/1520-0485(2001)031<0285:TROMAI>2.0.CO;2 (3×)
- `jiang:99`: Jiang et al. (1999). An assessment of the Geophysical Fluid Dynamics Laboratory ocean model with coarse resolution: Annual-mean climatology. J. Geophys. Res. 104, 25623-25645. doi:10.1029/1999JC900095 (3×)
- `jordan:18`: Jordan et al. (2018). Ocean-forced ice-shelf thinning in asynchronously coupled ice-ocean model. J. Geophys. Res. Oceans 123, 864-882. doi:10.1002/2017JC013251 (1×)
- `kalnay:96`: Kalnay et al. (1996). The NMC/NCAR 40-year reanalysis project. Bull. Am. Met. Soc. 77, 437-471. doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2 (1×)
- `keeling:98`: Keeling et al. (1998). Seasonal variations in the atmospheric O2/N2 ratio in relation to the kinetics of air-sea gas exchange. Global Biogeochemical Cycles 12, 141–163. doi:10.1029/97GB02339 (1×)
- `keir:80`: Keir (1980). The dissolution kinetics of biogenic calcium carbonates in seawater. Geochimica et Cosmochimica Acta 44, 241–252. doi:10.1016/0016-7037(80)90135-0 (1×)
- `kimmritz:15`: Kimmritz et al. (2015). On the Convergence of the Modified Elastic-Viscous-Plastic Method of Solving for Sea-Ice Dynamics. J. Comput. Phys. 296, 90–100. doi:10.1016/j.jcp.2015.04.051 (4×)
- `kimmritz:16`: Kimmritz et al. (2016). The Adaptive EVP Method for Solving the Sea Ice Momentum Equation. Ocean Modelling 101, 59–67. doi:10.1016/j.ocemod.2016.03.004 (6×)
- `klymaklegg10`: Klymak & Legg (2010). A simple mixing scheme for models that resolve breaking internal waves. Ocean Modelling 33, 224–234. doi:10.1016/j.ocemod.2010.02.005 (3×)
- `kondo:75`: Kondo (1975). Air-sea bulk transfer coefficients in diabatic conditions.. Bound. Layer Meteorol. 9, 91-112 (2×)
- `konig:10`: König Beatty & Holland (2010). Modeling Landfast Sea Ice by Adding Tensile Strength. J. Phys. Oceanogr. 40, 185–198. doi:10.1175/2009JPO4105.1 (1×)
- `ks:91`: Koster & Suarez (1991). A Simplified Treatment of SiB's Land Surface Albedo Parameterization.. National Aeronautics and Space Administration (1×)
- `ks:92`: Koster & Suarez (1992). Modeling the Land Surface Boundary in Climate Models as a Composite of Independent Vegetation Stands.. J. Geophys. Res. 97, 2697-2715. doi:10.1029/91JD01696 (1×)
- `lar-eta:94`: Large et al. (1994). Oceanic vertical mixing: A review and a model with nonlocal boundary layer parameterization. Rev. Geophys. 32, 363–403. doi:10.1029/94RG01872 (5×)
- `lar-eta:97`: Large et al. (1997). Sensitivity to surface forcing and boundary layer mixing in a global ocean model: Annual-mean climatology. J. Phys. Oceanogr. 27, 2418–2447. doi:10.1175/1520-0485(1997)027<2418:STSFAB>2.0.CO;2 (2×)
- `larpond:81`: Large & Pond (1981). Open ocean momentum flux measurements in moderate to strong winds. J. Phys. Oceanogr. 11, 324-336. doi:10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2 (2×)
- `lee:10`: Lee et al. (2010). The universal ratio of boron to chlorinity for the North Pacific and North Atlantic oceans. Geochimica et Cosmochimica Acta 74, 1801–1811. doi:10.1016/j.gca.2009.12.027 (1×)
- `leith:68`: Leith (1968). Large eddy simulation of complex engineering and geophysical flows. Physics of Fluids 10, 1409-1416 (2×)
- `leith:96`: Leith (1996). Stochastic models of chaotic systems. Physica D. 98, 481-491 (2×)
- `lemieux:10`: Lemieux et al. (2010). Improving the Numerical Convergence of Viscous-Plastic Sea Ice Models with the Jacobian-Free Newton-Krylov Method. J. Comput. Phys. 229, 2840–2852. doi:10.1016/j.jcp.2009.12.011c (2×)
- `lemieux:12`: Lemieux et al. (2012). A Comparison of the Jacobian-free Newton-Krylov Method and the EVP Model for Solving the Sea Ice Momentum Equation with a Viscous-Plastic Formulation: a Serial Algorithm Study. J. Comput. Phys. 231, 5926–5944. doi:10.1016/j.jcp.2012.05.024 (3×)
- `lemieux:15`: Lemieux et al. (2015). A Basal Stress Parameterization for Modeling Land-Fast Ice. J. Geophys. Res. 120, 3157–3173. doi:10.1002/2014JC010678 (1×)
- `leppaeranta:83`: Lepparanta (1983). A Growth Model for Black Ice, Snow Ican and Snow Thickness in Subarctic Basins. Nordic Hydrology 14, 59–70 (1×)
- `levitus:94a`: Levitus & Boyer (1994). World Ocean Atlas 1994 Volume 3: Salinity. National Oceanic and Atmospheric Administration (6×)
- `levitus:94b`: Levitus & Boyer (1994). World Ocean Atlas 1994 Volume 3: Temperature. National Oceanic and Atmospheric Administration. ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/PUBLICATIONS/WOA94_vol4a.pdf (6×)
- `lhans:74`: Lacis & Hansen (1974). A parameterization for the absorption of solar radiation in the Earth's atmosphere.. J. Atmos. Sci. 31, 118-133 (1×)
- `lipscomb:01`: Lipscomb (2001). Remapping the thickness distribution in sea ice models. J. Geophys. Res. 106, 13989–14000. doi:10.1029/2000JC000518 (3×)
- `lipscomb:07`: Lipscomb et al. (2007). Ridging, strength, and stability in high-resolution sea ice models. J. Geophys. Res. 112, 1–18. doi:10.1029/2005JC003355 (7×)
- `liu:22`: Liu et al. (2022). A New Parameterization of Coastal Drag to Simulate Landfast Ice in Deep Marginal Seas in the Arctic.. J. Geophys. Res. 127, e2022JC018413. doi:10.1029/2022JC018413 (1×)
- `losch:08`: Losch (2008). Modeling ice shelf cavities in a z-coordinate ocean general circulation model. J. Geophys. Res. Oceans 113, 129–144. doi:10.1029/2007JC004368 (2×)
- `losch:10`: Losch et al. (2010). On the formulation of sea-ice models. Part 1: Effects of different solver implementations and parameterizations. Ocean Modelling 33, 129–144. doi:10.1016/j.ocemod.2009.12.008 (1×)
- `losch:14`: Losch et al. (2014). A Parallel JAcobian-Free Newton-Krylov Solver for a Coupled Sea Ice-Ocean Model. J. Comput. Phys. 257, 901–910. doi:10.1016/j.jcp.2013.09.026 (3×)
- `lueker:00`: Lueker et al. (2000). Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. Marine Chemistry 70, 105–119. doi:10.1016/S0304-4203(00)00022-0 (1×)
- `Macayeal:89`: MacAyeal (1989). Large-scale ice flow over a viscous basal sediment: Theory and application to Ice Stream B, Antarctica. Journal of Geophysical Research – Solid Earth 94, 4071–4087 (1×)
- `mak:18`: Mak et al. (2018). Implementation of a geometrically informed and energetically constrained mesoscale eddy parameterization in an ocean circulation model. J. Phys. Oceanogr. 48, 2363–2382. doi:10.1175/JPO-D-18-0017.1 (1×)
- `mak:22`: Mak et al. (2022). Acute sensitivity of global ocean circulation and heat content to eddy energy dissipation time-scale. Geophys. Res. Lett. 49, e2021GL097259. doi:10.1029/2021GL097259 (1×)
- `manabe:79`: Manabe et al. (1979). A Global Ocean-Atmosphere Climate Model with Seasonal Variation for Future Studies of Climate Sensitivity. Dyn. Atmos. Oceans 3, 393–426. doi:10.1016/0377-0265(79)90021-6 (1×)
- `maro-eta:99`: Marotzke et al. (1999). Construction of the adjoint MIT ocean general circulation model and application to Atlantic heat transport variability. J. Geophys. Res. 104, 29,529–29,547. doi:10.1029/1999JC900236 (3×)
- `mars-eta:98`: Marshall et al. (1998). Efficient ocean modeling using non-hydrostatic algorithms. J. Mar. Sys. 18, 115–134. doi:10.1016/S0924-7963\%2898\%2900008-6 (1×)
- `marshall:03`: Marshall & Radko (2003). Residual-Mean Solutions for the Antarctic Circumpolar Current and Its Associated Overturning Circulation. J. Phys. Oceanogr. 33, 2341–2354. doi:10.1175/1520-0485(2003)033<2341:RSFTAC>2.0.CO;2 (1×)
- `marshall:04`: Marshall et al. (2004). Atmosphere-Ocean Modeling Exploiting Fluid Isomorphisms. Mon. Wea. Rev. 132, 2882-2894. doi:10.1175/MWR2835.1 (2×)
- `marshall:06`: Marshall et al. (2006). Estimates and implications of surface eddy diffusivity in the southern ocean derived from tracer transport. J. Phys. Oceanogr. 36, 1806-1821. doi:10.1175/JPO2949.1 (1×)
- `marshall:12`: Marshall & Speer (2012). Closure of the meridional overturning circulation through Southern Ocean upwelling. Nature Geosci. 5, 171–180. doi:10.1038/ngeo1391 (1×)
- `marshall:12b`: Marshall et al. (2012). A framework for parameterizing eddy potential vorticity fluxes. J. Phys. Oceanogr. 42, 539-557. doi:10.1175/JPO-D-11-048.1 (1×)
- `marshall:97a`: Marshall et al. (1997). Hydrostatic, quasi-hydrostatic, and nonhydrostatic ocean modeling. J. Geophys. Res. 102, 5733-5752. doi:10.1029/96JC02776 (13×)
- `marshall:97b`: Marshall et al. (1997). A finite-volume, incompressible Navier Stokes model for studies of the ocean on parallel computers. J. Geophys. Res. 102, 5753-5766. doi:10.1029/96JC02775 (2×)
- `martin:87`: Martin et al. (1987). VERTEX: carbon cycling in the northeast Pacific. Deep Sea Res. Part A. Oceanogr. Res. Papers 34, 267-285. doi:10.1016/0198-0149(87)90086-0 (1×)
- `mcdougall:03`: McDougall et al. (2003). Accurate and computationally efficient algorithms for potential temperature and density of seawater. J. Atmos. Ocean. Technol. 20, 730-741. doi:10.1175/1520-0426(2003)20<730:AACEAF>2.0.CO;2 (2×)
- `mcdougall:11`: McDougall & Barker (2011). Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox. SCOR/IAPSO WG127, 28 pp.. http://www.teos-10.org/pubs/gsw/pdf/Getting_Started.pdf (1×)
- `mckinley:04`: McKinley et al. (2004). Mechanisms of air-sea co2 flux variability in the equatorial pacific and the north atlantic. Global Biogeochem. Cycles 18. doi:10.1029/2003GB002179 (1×)
- `mehrbach:73`: Mehrbach et al. (1973). Measurement of the Apparent Dissociation Constants of Carbonic Acid in Seawater at Atmospheric Pressure. Limnology and Oceanography 18, 897–907. doi:10.4319/lo.1973.18.6.0897 (1×)
- `mellor:82`: Mellor & Yamada (1982). Development of a turbulence closure model for geophysical fluid problems. Rev. Geophys. 20, 851– 875. doi:10.1029/RG020i004p00851 (1×)
- `millero:10`: Millero (2010). History of the equation of state of seawater. Oceanography 23, 18-33. doi:10.5670/oceanog.2010.21 (1×)
- `millero:10b`: Millero (2010). Carbonate constants for estuarine waters. Marine and Freshwater Research 61, 139–142. doi:10.1071/MF09254 (1×)
- `millero:83`: Millero (1983). The estimation of the pK*HA of acids in seawater using the Pitzer equations. Geochimica et Cosmochimica Acta 47, 2121–2129. doi:10.1016/0016-7037(83)90037-6 (not cited in manual text)
- `millero:95`: Millero (1995). Thermodynamics of the carbon dioxide system in the oceans. Geochimica et Cosmochimica Acta 59, 661–677. doi:10.1016/0016-7037(94)00354-O (2×)
- `mol:09`: Molod (2009). Running GCM physics and dynamics on different grids: algorithm and tests. Tellus 61A, 381–393 (1×)
- `montagnes:1994`: Montagnes et al. (1994). Estimating carbon, nitrogen, protein, and chlorophyll a from volume in marine phytoplankton. Limnology and Oceanography 39, 1044–1060. doi:10.4319/lo.1994.39.5.1044 (1×)
- `moorsz:92`: Moorthi & Suarez (1992). Relaxed Arakawa Schubert: A parameterization of moist convection for general circulation models.. Mon. Wea. Rev. 120, 978-1002 (2×)
- `moum96`: Moum (1996). Energy-containing scales of turbulence in the ocean thermocline. J. Geophys. Res. 101, 14095–14109. doi:10.1029/96JC00507 (1×)
- `munhoven:13`: Munhoven (2013). Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388. doi:10.5194/gmd-6-1367-2013 (2×)
- `munk:50`: Munk (1950). On the wind-driven ocean circulation. J. Meteor. 7, 79-932 (2×)
- `naumann:06`: Naumann et al. (2006). Adjoint code by source transformation with OpenAD/F. European Conference on Computational Fluid Dynamics (ECCOMAS CFD 2006) (1×)
- `naviaux:19`: Naviaux et al. (2019). Calcite dissolution rates in seawater: Lab vs. in-situ measurements and inhibition by organic matter. Marine Chemistry 215, 103684. doi:10.1016/j.marchem.2019.103684 (1×)
- `neale:01`: Neale & Hoskins (2001). A standard test for AGCMs including their physical parametrizations: I: The proposal. Atmos. Sci. Letters. doi:10.1006/asle.2000.0022 (2×)
- `nikurashin:12`: Nikurashin & Vallis (2012). A Theory of the Interhemispheric Meridional Overturning Circulation and Associated Stratification. J. Phys. Oceanogr. 42, 1652–1667. doi:10.1175/JPO-D-11-0189.1 (1×)
- `olbers:04`: Olbers & Visbeck (2004). A Model of the Zonally Averaged Stratification and Overturning in the Southern Ocean. J. Phys. Oceanogr. 35, 1190–1205. doi:10.1175/JPO2750.1 (1×)
- `orl:76`: Orlanski (1976). A simple boundary condition for unbounded hyperbolic flows. J. Comput. Phys. 21, 251–269 (1×)
- `pacanowski:81`: Pacanowski & Philander (1981). Parameterization of Vertical Mixing in Numerical Models of Tropical Oceans. J. Phys. Oceanogr. 11, 1443–1451. doi:10.1175/1520-0485(1981)011<1443:POVMIN>2.0.CO;2 (1×)
- `pal-rom:97`: Paluszkiewicz & Romea (1997). A one-dimensional model for the parameterization of deep convection in the ocean. Dyn. Atmos. Oceans 26, 95–130 (2×)
- `pano:73`: Panofsky (1973). Tower Micrometeorology.. Workshop on Micrometeorology (1×)
- `parekh:2004`: Parekh et al. (2004). Modeling the global ocean iron cycle. Global Biogeochemical Cycles 18. doi:10.1029/2003GB002061 (1×)
- `parekh:2005`: Parekh et al. (2005). Decoupling of iron and phosphate in the global ocean. Global Biogeochemical Cycles 19. doi:10.1029/2004GB002280 (1×)
- `parkinson:79`: Parkinson & Washington (1979). A Large-Scale Numerical Model of Sea Ice. J. Geophys. Res. 84, 311–337. doi:10.1029/JC084iC01p00311 (1×)
- `patt:1994`: Patt & Gregg (1994). Exact closed-form geolocation algorithm for Earth survey sensors. International Journal of Remote Sensing 15, 3719–3734. doi:10.1080/01431169408954354 (not cited in manual text)
- `pawlowicz:13`: Pawlowicz (2013). Key Physical Variables in the Ocean: Temperature, Salinity, and Density. Nature Education Knowledge 4, 13. https://www.nature.com/scitable/knowledge/library/key-physical-variables-in-the-ocean-temperature-102805293/ (1×)
- `pedlosky:87`: Pedlosky (1987). Geophysical Fluid Dynamics, Second Edition. Spring-Verlag, 710 pp. (2×)
- `pedlosky:96`: Pedlosky (1996). Ocean Circulation Theory. Spring-Verlag, 456 pp. (1×)
- `perez:87`: Perez & Fraga (1987). Association constant of fluoride and hydrogen ions in seawater. Marine Chemistry 21, 161–168. doi:10.1016/0304-4203(87)90036-3 (1×)
- `potter:73`: Potter (1973). Computational Physics. John Wiley, 304 pp. (1×)
- `prather:86`: Prather (1986). Numerical advection by conservation of second-order moments. J. Geophys. Res. 91, 6671-6681. doi:10.1029/JD091iD06p06671 (1×)
- `redi1982`: Redi (1982). Oceanic Isopycnal Mixing by Coordinate Rotation. J. Phys. Oceanogr. 12, 1154–1158. doi:10.1175/1520-0485(1982)012<1154:OIMBCR>2.0.CO;2 (1×)
- `restrepo:98`: Restrepo et al. (1998). Circumventing storage limitations in variational data assimilation studies. SIAM J. Sci. Comput. 19, 1586-1605 (1×)
- `riley:65`: Riley (1965). The occurrence of anomalously high fluoride concentrations in the North Atlantic. Deep Sea Research and Oceanographic Abstracts 12, 219–220. doi:10.1016/0011-7471(65)90027-6 (1×)
- `ringeisen:19`: Ringeisen et al. (2019). Simulating Intersection Angles Between Conjugate Faults in Sea Ice with Different Viscous–Plastic Rheologies. The Cryosphere 13, 1167–1186. doi:10.5194/tc-13-1167-2019 (1×)
- `ringeisen:20`: Ringeisen et al. (2020). Non-Normal Flow Rules Affect Fracture Angles in Sea Ice Viscous-Plastic Rheologies. The Cryosphere Discussions 2020, 1–24. doi:10.5194/tc-2020-153 (1×)
- `rios:1998`: Ríos et al. (1998). Chemical composition of phytoplankton and Particulate Organic Matter in the Ría de Vigo (NW Spain). Sciencia Marina 6 2, 257–271. doi:10.3989/scimar.1998.62n3257 (not cited in manual text)
- `roach:15`: Roach et al. (2015). Detecting and Characterizing Ekman Currents in the Southern Ocean.. J. Phys. Oceanogr. 45, 1205–1223. doi:10.1175/JPO-D-14-0115.1 (1×)
- `roe:85`: Roe (1985). Some contributions to the modelling of discontinuous flows. Large-Scale Computations in Fluid Mechanics 22, 163–193 (2×)
- `roquet:15`: Roquet et al. (2015). Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling 90, 29-43. doi:10.1016/j.ocemod.2015.04.002 (1×)
- `rosen:87`: Rosenfield et al. (1987). A computation of the stratospheric diabatic circulation using an accurate radiative transfer model. J. Atmos. Sci. 44, 859-876 (1×)
- `rothrock:75`: Rothrock (1975). The energetics of the plastic deformation of pack ice by ridging. J. Geophys. Res. 80, 4514–4519. doi:10.1029/JC080i033p04514 (4×)
- `roy:93`: Roy et al. (1993). The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45°C. Marine Chemistry 44, 249–267. doi:10.1016/0304-4203(93)90207-5 (1×)
- `sadourny:75`: Sadourny (1975). The Dynamics of Finite-Difference Models of the Shallow-Water Equations. J. Atmos. Sci. 32, 680-689. doi:10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2 (1×)
- `sallee:18`: Sallée (2018). Southern Ocean Warming. Oceanography 31, 52 - 62. doi:10.5670/oceanog.2018.215 (1×)
- `seimgregg94`: Seim & Gregg (1994). Detailed observations of a naturally occurring shear instability. J. Geophys. Res. 99, 10049–10073. doi:10.1029/94JC00168 (1×)
- `semtner:76`: Semtner (1976). A model for the thermodynamic growth of sea ice in numerical investigations of climate. J. Phys. Oceanogr. 6, 379–389 (4×)
- `shapiro:70`: Shapiro (1970). Smoothing, filtering, and boundary effects. Rev. Geophys. Space Phys. 8, 359–387 (4×)
- `smag:63`: Smagorinsky (1963). General circulation experiments with the primitive equations I: The basic experiment. Mon. Wea. Rev. 91, 99–164 (2×)
- `smag:93`: Smagorinsky (1993). Large eddy simulation of complex engineering and geophysical flows. Evolution of Physical Oceanography, 3-36 (5×)
- `smith:19`: Smith & Heimbach (2019). Atmospheric origins of Variability in the South Atlantic Meridional Overturning Circulation. J. Clim. 32, 1483–1500. doi:10.1175/JCLI-D-18-0311.1 (1×)
- `smolark:89`: Smolarkiewicz (1989). Comment on "A positive definite advection scheme obtained by nonlinear renormalization of the advective fluxes". Mon. Wea. Rev. 117, 2626–2632. doi:10.1175/1520-0493(1989)117<2626:COPDAS>2.0.CO;2 (1×)
- `speer:00`: Speer & Sloyan (2000). The Diabatic Deacon Cell. J. Phys. Oceanogr. 30, 3212–3222. doi:10.1175/1520-0485(2000)030<3212:TDDC>2.0.CO;2 (1×)
- `stammer:02`: Stammer et al. (2002). The global ocean circulation and transports during 1992 - 1997, estimated from ocean observations and a general circulation model. J. Geophys. Res. 107, 3118. doi:10.1029/2001JC000888 (3×)
- `stammer:97`: Stammer et al. (1997). The global ocean circulation estimated from TOPEX/POSEIDON altimetry and a general circulation model. Massachusetts Institute of Technology. https://cgcs.mit.edu/publications/cgcs-report/global-ocean-circulation-estimated-topexposeidon-altimetry-and-mit-general (1×)
- `stevens:90`: Stevens (1990). On Open Boundary Conditions for Three Dimensional Primitive Equation Ocean Circulation Models. Geophys. Astrophys. Fl. Dyn. 51, 103–133 (3×)
- `stewart:17`: Stewart et al. (2017). Vertical resolution of baroclinic modes in global ocean models. Ocean Modelling 113, 50-65. doi:10.1016/j.ocemod.2017.03.012 (1×)
- `stommel:48`: Stommel (1948). The western intensification of wind-driven ocean currents. Trans. Am. Geophys. Union 29, 206 (2×)
- `sudm:88`: Sud & Molod (1988). The roles of dry convection, cloud-radiation feedback processes and the influence of recent improvements in the parameterization of convection in the GLA GCM.. Mon. Wea. Rev. 116, 2366-2387 (1×)
- `sulpis:22`: Sulpis et al. (2022). RADIv1: a non-steady-state early diagenetic model for ocean sediments in Julia and MATLAB/GNU Octave. Geoscientific Model Development 15, 2105–2131. doi:10.5194/gmd-15-2105-2022 (1×)
- `sverdrup:33`: Sverdrup (1933). On vertical circulation in the ocean due to the action of the wind with application to conditions within the Antarctic Circumpolar Current. Discovery Rept. 7, 141-169 (1×)
- `taksz:96`: Takacs & Suarez (1996). Dynamical aspects of climate simulations Using the GEOS General Circulation Model. National Aeronautics and Space Administration (1×)
- `thorndike:75`: Thorndike et al. (1975). The thickness distribution of sea ice. J. Geophys. Res. 80, 4501–4513 (6×)
- `thorpe77`: Thorpe (1977). Turbulence and mixing in a Scottish loch. Phil. Trans. R. Soc. Lond. 286, 125–181 (1×)
- `trenberth:89`: Trenberth et al. (1989). A global wind stress climatology based on ecmwf analyses. National Center for Atmospheric Research. doi:10.5065/D6ST7MR9 (1×)
- `trenberth:90`: Trenberth et al. (1990). The mean annual cycle in Global Ocean wind stress. J. Phys. Oceanogr. 20, 1742-1760. doi:10.1175/1520-0485(1990)020<1742:TMACIG>2.0.CO;2 (4×)
- `ungermann:17`: Ungermann et al. (2017). Impact of the Ice Strength Formulation on the Performance of a Sea Ice Thickness Distribution Model in the Arctic. J. Geophys. Res. 122, 2090–2107. doi:10.1002/2016JC012128 (2×)
- `uppstrom:74`: Uppström (1974). The boron/chlorinity ratio of deep-sea water from the Pacific Ocean. Deep Sea Research and Oceanographic Abstracts 21, 161–162. doi:10.1016/0011-7471(74)90074-6 (1×)
- `utke:08`: Utke et al. (2008). OpenAD/F: A modular open-source tool for automatic differentiation of Fortran codes. ACM Transactions on Mathematical Software (TOMS) 34, 18. doi:10.1145/1377596.1377598 (1×)
- `vallis:17`: Vallis (17). Atmospheric and Oceanic Fluid Dynamics: Fundamentals and Large-Scale Circulation, 2nd Edition. Cambridge University Press, 964 pp.. doi:10.1017/9781107588417 (3×)
- `veronis:75`: Veronis (1975). The role of models in tracer studies. Numerical Models of the Ocean Circulation, 133–146. https://books.google.com/books?hl=en&lr=&id=9S8rAAAAYAAJ&oi=fnd&pg=PA133&dq=Veronis,+G.,+1975:+The+role+of+models+in+tracer+studies.+Numerical+Models+of+the+Ocean+Circulation,+Natl.+Acad.+Sci.,+133%E2%80%93146.+&ots=xitpIWzXX3&sig=fXxEToFbFCutn1-7ZbMoY4oDEFE#v=onepage&q&f=false (1×)
- `viebahn:12`: Viebahn & Eden (2012). Standing Eddies in the Meridional Overturning Circulation. J. Phys. Oceanogr. 42, 1486 - 1508. doi:10.1175/JPO-D-11-087.1 (1×)
- `visbeck:97`: Visbeck et al. (1997). Specification of eddy transfer coefficients in coarse-resolution ocean circulation models. J. Phys. Oceanogr. 27, 381-402. doi:10.1175/1520-0485(1997)027<0381:SOETCI>2.0.CO;2 (1×)
- `wajsowicz:93`: Wajsowicz (1993). A consistent formulation of the anisotropic stress tensor for use in models of the large-scale ocean circulation. J. Comput. Phys. 105, 333-338 (1×)
- `wannink:92`: Wanninkhof (1992). Relationship between wind speed and gas exchange over the ocean. J. Geophys. Res. 97, 7373-7382. doi:10.1029/92JC00188 (3×)
- `waters:13`: Waters & Millero (2013). The free proton concentration scale for seawater pH. Marine Chemistry 149, 8–22. doi:10.1016/j.marchem.2012.11.003 (1×)
- `waters:14`: Waters et al. (2014). Corrigendum to “The free proton concentration scale for seawater pH”, [MARCHE: 149 (2013) 8–22]. Marine Chemistry 165, 66–67. doi:10.1016/j.marchem.2014.07.004 (1×)
- `weaver:01`: Weaver & Courtier (2001). Correlation modelling on the sphere using a generalized diffusion equation. Q. J. R. Meteorol. Soc. 127, 1815-1846. doi:10.1002/qj.49712757518 (1×)
- `weiss:80`: Weiss & Price (1980). Nitrous oxide solubility in water and seawater. Marine Chemistry 8, 347–359. doi:10.1016/0304-4203(80)90024-9 (1×)
- `wesson94`: Wesson & Gregg (1994). Mixing at Camarinal Sill in the Strait of Gibraltar. Q. J. R. Meteorol. Soc. 99 (C5), 9847–9878 (1×)
- `williams:69`: Williams (1969). Numerical integration of the three-dimensional Navier Stokes equations for incompressible flow. J. Fluid Mech. 37, 727–750 (2×)
- `winton:00`: Winton (2000). A reformulated three-layer sea ice model. J. Atmos. Ocean. Technol. 17, 525–531. doi:10.1175/1520-0426(2000)017<0525:ARTLSI>2.0.CO;2 (3×)
- `wolfe:14`: Wolfe (2014). Approximations to the ocean’s residual circulation in arbitrary tracer coordinates. Ocean Modelling 75, 20 - 35. doi:10.1016/j.ocemod.2013.12.004 (1×)
- `yagkad:74`: Yamada (1977). A numerical experiment on pollutant dispersion in a horizontally-homogeneous atmospheric boundary layer.. Atmos. Environ. 11, 1015-1024 (1×)
- `yamanaka:97`: Y. & Tajika (1997). Role of dissolved organic matter in the marine biogeochemical cycle: Studies using an ocean biogeochemical general circulation model. Global Biogeochem. Cycles 11, 599-612. doi:10.1029/97GB02301 (2×)
- `zha:05`: Zhang & Rothrock (2005). Effect of Sea Ice Rheology in Numerical Investigations of Climate. J. Geophys. Res. Oceans 110, C08014. doi:10.1029/2004JC002599 (2×)
- `zha:98`: Zhang et al. (1998). Arctic ice-ocean modeling with and without climate restoring. J. Phys. Oceanogr. 28, 191–217 (1×)
- `zhang:97`: Zhang & Hibler (1997). On an Efficient Numerical Method for Modeling Sea Ice Dynamics. J. Geophys. Res. 102, 8691–8702. doi:10.1029/96JC03744 (4×)
- `zhouetal:95`: Zhou et al. (1995). Impact of Orographically Induced Gravity Wave Drag in the GLA GCM. Q. J. R. Meteorol. Soc. 122, 903-927 (1×)

# Literature behind ECCO, ECCO-Darwin and Dustin's MITgcm work

Papers that the MITgcm/darwin3 manuals don't cite. For a scheme or package's own paper (seaice, KPP,
GGL90, GMRedi, z*, cubed sphere, TAF, Darwin ecosystem, radtrans…), grep `references/mitgcm-index/bibliography.md`,
which is generated from the manuals' bib files and says which manual section and parameter cites each paper.

Every DOI/handle below was resolved against Crossref, DataCite, or hdl.handle.net on 2026-10-05; Dustin's
papers come from his ORCID record (0000-0003-1686-5255). Notes say what each paper is used for here; they
aren't summaries. When citing, re-check details (pages, final versions) against the DOI.

## ECCO state estimation (framework, releases)

- Wunsch & Heimbach (2007). Practical global oceanic state estimation. Physica D 230, 197–208. doi:10.1016/j.physd.2006.09.040
  — the adjoint/least-squares philosophy behind ECCO (state estimate = free forward run of the model).
- Wunsch et al. (2009). The Global General Circulation of the Ocean Estimated by the ECCO-Consortium. Oceanography 22, 88–103. doi:10.5670/oceanog.2009.41
- Heimbach et al. (2005). An efficient exact adjoint of the parallel MIT General Circulation Model, generated via automatic differentiation. Future Generation Computer Systems 21, 1356–1371. doi:10.1016/j.future.2004.11.010
  — TAF adjoint of MITgcm: checkpointing, store directives, parallel adjoint (pairs with `adjoint.md` troubleshooting).
- Forget et al. (2015). ECCO version 4: an integrated framework for non-linear inverse modeling and global ocean state estimation. Geosci. Model Dev. 8, 3071–3104. doi:10.5194/gmd-8-3071-2015
  — **the** ECCO v4 reference: LLC90 grid, model settings, ctrl/cost set-up. Cite for v4 configurations.
- Fukumori et al. (2017). ECCO Version 4 Release 3. MIT DSpace technical note. hdl:1721.1/110380 (https://hdl.handle.net/1721.1/110380)
- ECCO Consortium et al. (2021). Synopsis of the ECCO Central Production Global Ocean and Sea-Ice State Estimate, Version 4 Release 4. Zenodo. doi:10.5281/zenodo.4533349
  — v4r4 (1992–2017, checkpoint66g); the physics behind ECCO-Darwin v05/v06 1° and LLC90 runs is v4r5, built on this lineage.
- Fenty & Heimbach (2013). Coupled Sea Ice–Ocean-State Estimation in the Labrador Sea and Baffin Bay. J. Phys. Oceanogr. 43, 884–904. doi:10.1175/jpo-d-12-065.1
  — sea-ice adjoint/state estimation (seaice ctrl variables).
- Nguyen et al. (2021). The Arctic Subpolar Gyre sTate Estimate (ASTE): Description and Assessment of a Data-Constrained, Dynamically Consistent Ocean-Sea Ice Estimate for 2002–2017. JAMES 13. doi:10.1029/2020ms002398
  — regional (Arctic/subpolar) ECCO-style estimate; physics for ASTE-BGC.

## ECCO-Darwin and BGC state estimation

- Brix et al. (2015). Using Green's Functions to initialize and adjust a global, eddying ocean biogeochemistry general circulation model. Ocean Modelling 95, 1–14. doi:10.1016/j.ocemod.2015.07.008
  — Green's-function optimisation of BGC parameters/ICs used for ECCO-Darwin.
- Follows et al. (2007). Emergent Biogeography of Microbial Communities in a Model Ocean. Science 315, 1843–1846. doi:10.1126/science.1138544
  — origin of the Darwin ecosystem approach (Darwin's equations themselves: darwin3 manual, `bibliography.md`).
- Carroll et al. (2020). The ECCO-Darwin Data-Assimilative Global Ocean Biogeochemistry Model: Estimates of Seasonal to Multidecadal Surface Ocean pCO2 and Air-Sea CO2 Flux. JAMES 12. doi:10.1029/2019ms001888
  — **the** ECCO-Darwin reference (configuration in `ecco_darwin/v04/llc270_JAMES_paper`).
- Carroll et al. (2022). Attribution of Space-Time Variability in Global-Ocean Dissolved Inorganic Carbon. Global Biogeochem. Cycles 36. doi:10.1029/2021gb007162
  — DIC budget closure in ECCO-Darwin (`ecco_darwin/v04/llc270_JAMES_budget`; budget recipes in `troubleshooting/diagnostics.md`).
- Savelli et al. (2026). Implementing riverine biogeochemical inputs in ECCO-Darwin: a sensitivity analysis of terrestrial fluxes in a data-assimilative global ocean biogeochemistry model. Geosci. Model Dev. 19, 867–885. doi:10.5194/gmd-19-867-2026
  — river carbon/nutrient forcing (v06 "+ rivers").
- Suselj et al. (2025). Quantifying Marine Carbon Dioxide Removal via Alkalinity Enhancement Across Circulation Regimes Using ECCO-Darwin and 1D Models. JAMES 17. doi:10.1029/2024ms004847
- Tyka et al. (2026). Substantial inter-model variation in OAE efficiency between the CESM2/MARBL and ECCO-Darwin ocean biogeochemistry models. Biogeosciences 23, 4943–4966. doi:10.5194/bg-23-4943-2026
  — OAE experiments (`ecco_darwin/v05/1deg_oaemip`, `3deg_CDR`).
- van der Zant et al. (2026). RADIv2: an adaptable and versatile diagenetic model for coastal and open-ocean sediments. Geosci. Model Dev. 19, 1965–1989. doi:10.5194/gmd-19-1965-2026
  — sediment model coupled in `ecco_darwin/v05/1deg_RADIv2` and `llc270_RADIv1`.
- Moseley et al. (2026). The ASTE-BGC Data-Assimilative Regional Ocean Biogeochemical Model. JAMES 18. doi:10.1029/2025ms004976
- Choi et al. (2024). A New Ecosystem Model for Arctic Phytoplankton Phenology From Ice-Covered to Open-Water Periods. Geophys. Res. Lett. 51. doi:10.1029/2024gl110155
- Yasunaka et al. (2023). An Assessment of CO2 Uptake in the Arctic Ocean From 1985 to 2018. Global Biogeochem. Cycles 37. doi:10.1029/2023gb007806

## Regional / river-plume configurations (Mackenzie)

- Bertin et al. (2023). Biogeochemical River Runoff Drives Intense Coastal Arctic Ocean CO2 Outgassing. Geophys. Res. Lett. 50. doi:10.1029/2022gl102377
  — Mackenzie shelf set-up (`ecco_darwin/regions/mac_delta`, ED-SBS ecosystem).
- Bertin et al. (2025). Paving the Way for Improved Representation of Coupled Physical and Biogeochemical Processes in Arctic River Plumes—A Case Study of the Mackenzie Shelf. Permafrost Periglac. Process. 36, 363–377. doi:10.1002/ppp.2271
- Bertin et al. (2025). Colored dissolved organic matter (CDOM) alters the seasonal physics and biogeochemistry of the Arctic Mackenzie River plume. Biogeosciences 22, 6607–6629. doi:10.5194/bg-22-6607-2025
  — CDOM/radtrans coupling in the plume.

## Glacier fjords and subglacial plumes (MITgcm, pre-ECCO work)

- Carroll et al. (2015). Modeling Turbulent Subglacial Meltwater Plumes: Implications for Fjord-Scale Buoyancy-Driven Circulation. J. Phys. Oceanogr. 45, 2169–2185. doi:10.1175/jpo-d-15-0033.1
- Carroll et al. (2016). The impact of glacier geometry on meltwater plume structure and submarine melt in Greenland fjords. Geophys. Res. Lett. 43, 9739–9748. doi:10.1002/2016gl070170
- Carroll et al. (2017). Subglacial discharge-driven renewal of tidewater glacier fjords. J. Geophys. Res. Oceans 122, 6611–6629. doi:10.1002/2017jc012962
- Jackson et al. (2017). Near-glacier surveying of a subglacial discharge plume: Implications for plume parameterizations. Geophys. Res. Lett. 44, 6886–6894. doi:10.1002/2017gl073602
- Amundson & Carroll (2018). Effect of Topography on Subglacial Discharge and Submarine Melting During Tidewater Glacier Retreat. J. Geophys. Res. Earth Surf. 123, 66–79. doi:10.1002/2017jf004376
- Slater et al. (2022). Characteristic Depths, Fluxes, and Timescales for Greenland's Tidewater Glacier Fjords From Subglacial Discharge-Driven Upwelling During Summer. Geophys. Res. Lett. 49. doi:10.1029/2021gl097081

## Adjoint tools

- Giering & Kaminski (1998), TAF recipes, and the other TAF papers: in `bibliography.md` (`giering:98`, `giering:99`, `giering:00`).
- Hascoët & Pascual (2013). The Tapenade automatic differentiation tool. ACM Trans. Math. Softw. 39, 1–43. doi:10.1145/2450153.2450158
  — Tapenade (MITgcm `-tap` builds; see `build.md` and `troubleshooting/adjoint.md`).

## Updating

Add a paper only after resolving its DOI (`curl -s https://api.crossref.org/works/<doi>`; DataCite/Zenodo DOIs via
`curl -sL -H 'Accept: application/vnd.citationstyles.csl+json' https://doi.org/<doi>`). Dustin's new papers:
`curl -s -H 'Accept: application/json' https://pub.orcid.org/v3.0/0000-0003-1686-5255/works`.

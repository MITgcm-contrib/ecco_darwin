# Observations for 1-D ECCO-Darwin columns (README_obs)

Assembled 2026-10-07 (gap fill for HOT, PAP and GEOTRACES Fe on 2026-10-08) for comparison against 1-D ECCO-Darwin columns (1992-2025, 50 levels, daily).
Sites (model column): HOT/ALOHA (22.75N 158.0W), BATS (31.67N 64.17W), Hydrostation S (32.17N 64.5W),
Papa / Line P P26 (50.1N 144.9W), PAP-SO (49.0N 16.5W).

## Layout

- `<SITE>/raw/<source>/` raw downloads, one dir per source (untrusted data; never executed)
- `<SITE>/<SITE>_obs.csv` tidy long-format file per site (SITE = HOT, BATS, HydroS, Papa, PAP)
- `scripts/standardize_obs.py` driver; per-source readers `std_hot.py`, `std_bios.py` (BATS + HydroS),
  `std_papa.py`, `std_pap.py`, `std_synth.py` (GLODAP/SOCAT/BGC-Argo/satellite for all sites),
  `std_tm.py` (dissolved trace metals, HOT + BATS; added 2026-10-08);
  `fetch_synth.sh`, `subset_glodap.py` (synthesis download/subset helpers); `summarize_core.py` (core summary)
- `obs_summary.csv` (site x variable x source) and `obs_summary_core.csv` (site x variable)

Regenerate: `python3 -I obs/scripts/standardize_obs.py --summary-csv obs/obs_summary.csv` then
`python3 -I obs/scripts/summarize_core.py obs/obs_summary.csv obs/obs_summary_core.csv`
(the driver re-adds the user site-packages dir, where pandas/numpy/netCDF4 live, since `-I` drops it).

## Tidy file columns

`time` (ISO 8601 UTC), `depth_m` (positive down; pressure in dbar converted or taken ~ m, see source notes),
`variable`, `value` (as reported, no unit conversion), `units` (original), `source` (= raw/<source> tag),
`qc` (original flag of the kept sample, '' if none), `lat`, `lon` (sample position; blank = at station).

Variable names: `temp` (in-situ T; all sources) and `theta` (potential T; HOT bottle file only), `salt`, `NO3`, `NO2`,
`NH4`, `PO4`, `SiO2`, `DIC`, `ALK`, `O2`, `Chl` (fluorometric extracted or sensor fluorescence), `Chl_HPLC`,
`HPLC_<pigment>`, `Chl_sat`, `PP`, `pCO2`, `fCO2`, `pCO2_air`, `pH` (scale/temperature in units column), `MLD`
(no site provides one), `POC`, `PON`, `POP`, `PIC`, `DOC`/`TOC`, `DON`, `DOP`, `TDN`, `TDP`, `bSi`, `N2O`,
`POC_flux`, `PON_flux`, `POP_flux`, `PIC_flux`, `bSi_flux`, `mass_flux`, `lith_flux`, `bbp`, cell counts
(`BactAbund`, `Prochlorococcus`/`ProAbund`, etc.), plus a few source-specific extras (`O2_CTD`, `NO3_LL`, `PO4_LL`,
`SRP_LL`, `Phaeo`, `d13C_DIC`, `PP_dark`, `lSi`), and dissolved trace metals `Fe`, `Mn`, `Co`, `Ni`, `Cu`, `Zn`, `Cd`,
`Pb`, `Al`, `Ti` (added 2026-10-08). See the per-source notes below.

## Flag convention (summary; details per source below)

Kept: HOT bottle 1 (unchecked) + 2 (good); BIOS 1 (unverified) + 2 (verified), drop bottle flag -3 and 3/4/9;
WOCE-style (GLODAP f, SOCAT, NCEI OCADS, Line P carbon) 2 (+6 for Line P replicate means); SOCAT dataset QC A-D;
OceanSITES 0 (no QC, range-checked), 1, 2 at PAP and 1, 2 at Papa; OOI QARTOD 1-2; BGC-Argo adjusted with QC 1/2/5/8.
Missing values (-9, -999, NaN, fill) dropped. IOS bottle ERDDAP (Papa) and PANGAEA have no flags (qc blank).
Added 2026-10-08: SPOTS (HOT) WOCE 2 + 6; HOT PP v5 bottle + light flags 1/2; HOT trap v4 no flags; trace metals
SeaDataNet/GEOTRACES 1 + 2 (HOT 962986, BAIT Fe 936824) and 1 only (BAIT bottle 937302, own scheme).
Newer versions of existing HOT sources are de-duplicated BY CRUISE: they only add cruises missing from the v1 file.

## Units and density conversion

Model tracers are mmol/m^3. Per-kg values need x rho/1000 (rho ~ 1.022-1.050 kg/L): nutrients, DIC, ALK, O2 at HOT,
BATS, HydroS, GLODAP, Line P carbon, PAP mooring O2, Argo O2/NO3 (umol/kg). Per-litre values (umol/L = mmol/m^3,
no conversion): Line P/IOS nutrients, PAP JC087 bottle nutrients/O2. Old IOS O2 in mL/L: x 44.66 to umol/L.
Chl is ug/L or ug/kg depending on source (BATS Chl ug/kg, Chl_HPLC ng/kg); POC/PON at BATS ug/kg. Always read
the `units` column.

## General caveats

- Synthesis subsets (GLODAP, SOCAT +/-1 deg; BGC-Argo within 100 km) are off-column; use lat/lon to filter.
  BATS and HydroS are ~60 km apart, so their synthesis subsets overlap.
- SOCAT contains the moorings (WHOTS, BTM, Papa, PAP), so it partly duplicates mooring pCO2/fCO2 records.
  Deduplicate by source before combining. PAP JC278 underway data are probably also in SOCAT.
- Several records pre-date 1992 (model start); kept for climatologies.
- Moored records are high frequency (3-hourly to 30-min); PAP_obs.csv is ~176 MB mostly from 2002-05 T/S chains.
  Average to daily before comparing.
- All times are UTC (HOT PP times corrected from local HST, +10 h; date-only records set to 00:00; trap records at
  deployment midpoint).

## Summary (site x variable, all sources combined; full breakdown in obs_summary.csv)

| site | variable | n | years | depth (m) | sources |
|---|---|---|---|---|---|
| BATS | ALK | 8508 | 1989-2025 | 1-4546 | bcodmo_bottle,glodap,scripps_co2 |
| BATS | Chl | 27817 | 1988-2025 | 0-2002 | bcodmo_pigments,bgcargo |
| BATS | Chl_HPLC | 6120 | 1989-2025 | 0-504 | bcodmo_pigments |
| BATS | Chl_sat | 346 | 1997-2026 | 0-0 | satchl |
| BATS | DIC | 9872 | 1988-2025 | 0-4546 | bcodmo_bottle,glodap,scripps_co2 |
| BATS | DOC | 140 | 2003-2021 | 3-4525 | glodap |
| BATS | Fe | 284 | 2019-2019 | 0-1678 | bcodmo_bait_fe,bcodmo_bait_tm_bottle |
| BATS | NO3 | 25176 | 1988-2025 | 0-4556 | bcodmo_bottle,bgcargo,glodap |
| BATS | O2 | 50165 | 1988-2025 | 0-4728 | bcodmo_bottle,bgcargo,glodap |
| BATS | PIC_flux | 1877 | 1978-2018 | 500-3200 | ncei_ofp |
| BATS | PO4 | 19313 | 1988-2025 | 0-4556 | bcodmo_bottle,glodap |
| BATS | POC | 9348 | 1988-2025 | 0-1011 | bcodmo_bottle |
| BATS | POC_flux | 3082 | 1978-2025 | 150-3200 | bcodmo_flux,ncei_ofp |
| BATS | PON | 9394 | 1988-2025 | 0-1011 | bcodmo_bottle |
| BATS | PP | 3684 | 1988-2025 | 0-161 | bcodmo_pp |
| BATS | SiO2 | 19592 | 1988-2025 | 0-4556 | bcodmo_bottle,glodap |
| BATS | bSi_flux | 675 | 2000-2015 | 500-3200 | ncei_ofp |
| BATS | bbp | 21961 | 2022-2025 | 1-2002 | bgcargo |
| BATS | fCO2 | 214945 | 1993-2025 | 0-5 | ncei_btm_mooring,socat |
| BATS | pH | 12702 | 2012-2024 | 1-3997 | bgcargo,glodap |
| BATS | salt | 322757 | 1988-2025 | 0-4728 | bcodmo_bottle,bgcargo,glodap,ncei_btm_mooring,socat |
| BATS | temp | 326482 | 1988-2025 | 0-4728 | bcodmo_bottle,bgcargo,glodap,ncei_btm_mooring,socat |
| HOT | ALK | 6602 | 1988-2024 | 1-4728 | bcodmo_bottle,bcodmo_spots,glodap,hot_co2,scripps_co2 |
| HOT | Chl | 30073 | 1988-2025 | 0-1601 | bcodmo_bottle,bgcargo,ncei_whots_pco2 |
| HOT | Chl_HPLC | 3098 | 1988-2016 | 2-1014 | bcodmo_bottle |
| HOT | Chl_sat | 344 | 1997-2026 | 0-0 | satchl |
| HOT | DIC | 6736 | 1988-2024 | 1-4728 | bcodmo_bottle,bcodmo_spots,glodap,hot_co2,scripps_co2 |
| HOT | DOC | 4757 | 1993-2017 | 3-4715 | bcodmo_bottle,bcodmo_spots |
| HOT | Fe | 193 | 2020-2023 | 15-900 | bcodmo_hot_metals |
| HOT | NO3 | 48191 | 1984-2023 | 1-4730 | bcodmo_bottle,bcodmo_spots,bgcargo,glodap |
| HOT | O2 | 103831 | 1984-2026 | 0-4730 | bcodmo_bottle,bcodmo_spots,bgcargo,glodap,ncei_whots_pco2 |
| HOT | PIC_flux | 178 | 2001-2024 | 150-150 | bcodmo_flux,bcodmo_flux_v4 |
| HOT | PO4 | 17957 | 1984-2024 | 1-4730 | bcodmo_bottle,bcodmo_spots,glodap,hot_co2 |
| HOT | POC | 3649 | 1989-2019 | 1-4426 | bcodmo_bottle,bcodmo_spots |
| HOT | POC_flux | 528 | 1988-2024 | 70-500 | bcodmo_flux,bcodmo_flux_v4 |
| HOT | PON | 3626 | 1989-2019 | 1-4426 | bcodmo_bottle,bcodmo_spots |
| HOT | PP | 2251 | 1988-2024 | 2-178 | bcodmo_pp,bcodmo_pp_v5 |
| HOT | SiO2 | 18002 | 1984-2024 | 1-4730 | bcodmo_bottle,bcodmo_spots,glodap,hot_co2 |
| HOT | bSi_flux | 193 | 1997-2024 | 150-150 | bcodmo_flux,bcodmo_flux_v4 |
| HOT | bbp | 19038 | 2013-2026 | 5-1701 | bgcargo |
| HOT | fCO2 | 92602 | 1993-2025 | 0-5 | ncei_whots_pco2,socat |
| HOT | pCO2 | 40318 | 1988-2025 | 0-15 | hot_co2,ncei_whots_pco2 |
| HOT | pH | 34404 | 1992-2026 | 0-4725 | bcodmo_bottle,bcodmo_spots,bgcargo,hot_co2,ncei_whots_pco2 |
| HOT | salt | 240120 | 1984-2026 | 0-4730 | bcodmo_bottle,bcodmo_spots,bgcargo,glodap,hot_co2,ncei_whots_pco2,scripps_co2,socat |
| HOT | temp | 246121 | 1984-2026 | 0-4730 | bcodmo_bottle,bcodmo_spots,bgcargo,glodap,hot_co2,ncei_whots_pco2,scripps_co2,socat |
| HOT | theta | 77978 | 1988-2016 | 0-4730 | bcodmo_bottle |
| HydroS | ALK | 764 | 1983-2024 | 0-4718 | glodap,scripps_co2 |
| HydroS | Chl | 30614 | 2022-2024 | 1-2002 | bgcargo |
| HydroS | Chl_sat | 346 | 1997-2026 | 0-0 | satchl |
| HydroS | DIC | 806 | 1983-2024 | 0-4719 | glodap,scripps_co2 |
| HydroS | DOC | 232 | 2003-2021 | 3-4719 | glodap |
| HydroS | NO3 | 8522 | 2003-2024 | 3-4719 | bgcargo,glodap |
| HydroS | O2 | 73680 | 1955-2025 | 0-4719 | bcodmo_bottle,bgcargo,glodap |
| HydroS | PO4 | 456 | 2003-2021 | 3-4719 | glodap |
| HydroS | SiO2 | 455 | 2003-2021 | 3-4719 | glodap |
| HydroS | bbp | 30219 | 2022-2024 | 1-2002 | bgcargo |
| HydroS | fCO2 | 218928 | 1993-2025 | 5-5 | socat |
| HydroS | pH | 18624 | 2012-2024 | 1-4719 | bgcargo,glodap |
| HydroS | salt | 303752 | 1955-2025 | 0-4719 | bcodmo_bottle,bgcargo,glodap,socat |
| HydroS | temp | 303836 | 1955-2025 | 0-4719 | bcodmo_bottle,bgcargo,glodap,socat |
| PAP | ALK | 236 | 1981-2017 | 0-4905 | glodap |
| PAP | Chl | 49039 | 1991-2025 | 0-2015 | bgcargo,glodap,oceansites_gdac_ifremer |
| PAP | Chl_sat | 338 | 1997-2026 | 0-0 | satchl |
| PAP | DIC | 228 | 1981-2017 | 2-4905 | glodap |
| PAP | NO3 | 2327 | 1981-2023 | 0-4905 | bgcargo,glodap,pangaea_jc087 |
| PAP | O2 | 64906 | 1981-2025 | 0-4905 | bgcargo,glodap,oceansites_gdac_ifremer,pangaea_jc087 |
| PAP | PIC_flux | 220 | 1995-2005 | 1000-4700 | pangaea_traps |
| PAP | PO4 | 719 | 1981-2017 | 0-4825 | glodap,pangaea_jc087 |
| PAP | POC_flux | 254 | 1995-2005 | 1000-4700 | pangaea_traps |
| PAP | SiO2 | 965 | 1981-2017 | 0-4905 | glodap,pangaea_jc087 |
| PAP | bSi_flux | 198 | 1996-2005 | 1000-4700 | pangaea_traps |
| PAP | bbp | 16718 | 2016-2025 | 0-2015 | bgcargo |
| PAP | fCO2 | 43431 | 1990-2025 | 5-5 | ncei_ocads_jc278,socat |
| PAP | pCO2 | 13430 | 2010-2025 | 1-35 | ncei_ocads_jc278,ncei_ocads_mooring,oceansites_gdac_ifremer |
| PAP | pH | 185 | 1991-2003 | 0-4905 | glodap |
| PAP | salt | 880941 | 1981-2025 | 0-4905 | bgcargo,glodap,ncei_ocads_jc278,oceansites_gdac_ifremer,socat |
| PAP | temp | 978444 | 1981-2025 | 0-4905 | bgcargo,glodap,ncei_ocads_jc278,oceansites_gdac_ifremer,socat |
| Papa | ALK | 1650 | 1992-2021 | 0-4326 | glodap,ncei_linep_carbon |
| Papa | Chl | 54237 | 1958-2025 | 0-2002 | bgcargo,glodap,ios_bot_erddap,ooi_papa_flanking_flma,ooi_papa_flanking_flmb,pmel_co2_mooring |
| Papa | Chl_sat | 327 | 1997-2026 | 0-0 | satchl |
| Papa | DIC | 2724 | 1986-2021 | 0-4326 | glodap,ncei_linep_carbon |
| Papa | DOC | 28 | 2014-2014 | 10-4237 | glodap |
| Papa | NO3 | 36830 | 1959-2024 | 0-4258 | bgcargo,glodap,ios_bot_erddap |
| Papa | O2 | 144527 | 1956-2025 | 0-4300 | bgcargo,glodap,ios_bot_erddap,ooi_papa_flanking_flma,ooi_papa_flanking_flmb,pmel_co2_mooring |
| Papa | PIC_flux | 17 | 1983-1993 | 200-3800 | pangaea_osp_trap |
| Papa | PO4 | 9924 | 1956-2024 | 0-4258 | glodap,ios_bot_erddap |
| Papa | POC_flux | 17 | 1983-1993 | 200-3800 | pangaea_osp_trap |
| Papa | SiO2 | 10262 | 1958-2024 | 0-4300 | glodap,ios_bot_erddap |
| Papa | bSi_flux | 17 | 1983-1993 | 200-3800 | pangaea_osp_trap |
| Papa | bbp | 40320 | 2010-2024 | 2-2002 | bgcargo |
| Papa | fCO2 | 73472 | 1990-2025 | 5-5 | socat |
| Papa | pCO2 | 45004 | 2007-2025 | 0-0 | pmel_co2_mooring |
| Papa | pH | 65864 | 2005-2025 | 0-4237 | bgcargo,glodap,ooi_papa_flanking_flma,ooi_papa_flanking_flmb,pmel_co2_mooring |
| Papa | salt | 343389 | 1956-2026 | 0-4300 | bgcargo,glodap,ios_bot_erddap,oceansites_pmel_daily,pmel_co2_mooring,socat |
| Papa | temp | 365675 | 1956-2026 | 0-4300 | bgcargo,glodap,ios_bot_erddap,oceansites_pmel_daily,pmel_co2_mooring,socat |

## Missing data and sources needing login or manual fetch

- HOT (updated 2026-10-08): PP and 150 m trap fluxes now run to Dec 2024 (BCO-DMO 737163 v5, 737393 v4) and bottle
  nutrients/O2/DIC/ALK/pH/particulate C-N-P to Dec 2019 (SPOTS). **Still missing: HOT bottle chemistry 2020-2024, and all
  HPLC/fluorometric Chl/bSi/N2O/DON/DOP/cell counts after 2016.** BCO-DMO 3773 v4 is "Under revision" (no file), the
  3773_v3 ERDDAP needs a login (401), WHOAS v3 is CAPTCHA-gated and SOEST FTP/HOT-DOGS time out. Manual route: see
  HOT > "Not obtained". The HOT 2-dbar CTD (BCO-DMO 3937) and deep moored traps were not fetched.
- Line P cruise data on waterproperties.ca needs an account (https://www.waterproperties.ca/linep/cruises.php);
  used CIOOS ERDDAP + NCEI OCADS instead. No open Line P PP or Timothy et al. 2013 trap series.
- PAP (re-tried 2026-10-08, nothing new): everything after Jul 2019, nutrient/pH sensors, later pCO2, cruise
  CTD/nutrients and traps after 2005 are in BODC (registration required; www.bodc.ac.uk also unresolvable from here):
  https://www.bodc.ac.uk/data/bodc_database/nodb/search/ (EDMED 5912). Copernicus Marine in-situ needs registration.
  EMODnet Physics ERDDAP timed out. There is still no open PAP nitrate or pH time series (see PAP > "Not obtained").
- GEOTRACES IDP2025 dissolved Fe (re-tried 2026-10-08): BODC host unresolvable from here, webODV interactive only. Open
  BCO-DMO substitutes were added for HOT (dFe 2020-2023) and BATS (BAIT dFe 2019). **No dFe for Papa or PAP.** For GA03,
  GP15, Canadian GEOTRACES P26 and GEOVIDE, see Synthesis > GEOTRACES for the manual steps.
- BATS OFP deep traps: open file (v2) ends 2018; PI (M. Conte) asks users to contact her before use.
- No site provides an MLD product; compute from T/S profiles if needed.

---

# Per-source notes


<!-- from HOT/SOURCES_HOT.md -->
## HOT / Station ALOHA (22.75N, 158.0W) - sources

All accessed 2026-10-07 (sources 7-10 on 2026-10-08). Readers: `obs/scripts/std_hot.py` (+ `std_tm.py` for metals) -> `obs/HOT/HOT_obs.csv`.

### Access problems (read first; re-tried 2026-10-08)
- SOEST servers still unreachable on 2026-10-08: TCP connects to `ftp.soest.hawaii.edu` (ports 21/80/443; the
  water, ctd, primary_production and particle_flux dirs linked from https://hahana.soest.hawaii.edu/hot/dataaccess.html)
  and to `scopeserver.soest.hawaii.edu` (the HOT-DOGS extraction CGI `hd-cgi-bin/mymatweb` that
  https://hahana.soest.hawaii.edu/hot/hot-dogs/bextraction.html posts to) time out. `hahana.soest.hawaii.edu` itself
  works (HOT CO2 product). HOT-DOGS also has a password-only "preliminary data" area; it is not public and was not touched.
- HOT bottle (BCO-DMO 3773): v4 (1988-2024, version date 2026-07-07) is "Under revision" on
  https://www.bco-dmo.org/dataset/3773 with no data file listed. ERDDAP `bcodmo_dataset_3773_v3` (1988-2023) returns
  HTTP 401 ("log in"). The v3 WHOAS archive copy (https://doi.org/10.26008/1912/bco-dmo.3773.3 ->
  darchive.mblwhoilibrary.org) sits behind a CAPTCHA "Security Check". None of these were bypassed.
- Gap fill (2026-10-08): PP v5 and trap flux v4 are open as direct files on datadocs.bco-dmo.org and now run to
  HOT-355 (Dec 2024) (sources 7 and 8). Bottle chemistry for HOT-289 to 317 (Jan 2017 to Dec 2019) comes from SPOTS
  (source 9). **HOT bottle chemistry for 2020-2024 is still missing** (see "Not obtained").

### 1. raw/bcodmo_bottle/ - HOT Niskin bottle + CTD-at-bottle data
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_3773.csv?&STNNBR=2 (Station 2 = ALOHA only)
- Files: `hot_bottle_bcodmo3773_stn2.csv`. That download was cut off on the slow link at 2011-09-27, and the truncated row
  is dropped. `hot_bottle_bcodmo3773_stn2_part2_YYYY.csv` holds yearly chunks for 2011-2016 (overlap de-duplicated).
- Citation: Karl, D. M., Fujieki, L. A. et al. Niskin bottle water samples and CTD measurements from the
  Hawaii Ocean Time-Series cruises (HOT project). BCO-DMO. Version 1 (check author list on the landing page). doi:10.1575/1912/bco-dmo.3773.1 (CC-BY 4.0).
  Plus the HOT-requested ack: "Data obtained via the Hawaii Ocean Time-series HOT-DOGS application; University of
  Hawai'i at Manoa. National Science Foundation Award # 2241005."
- Span 1988-10 to 2016-11, 0 to ~4730 m.
- Variables: temp (CTDTMP, in-situ ITS-90), theta (THETA), salt (CTDSAL), O2 (bottle Winkler), O2_CTD, NO3 (= NO3+NO2),
  NO2 (nmol/kg), NO3_LL / PO4_LL (low-level nmol/kg, upper 200 m), PO4, SiO2, DIC, ALK (ueq/kg), pH (TOT25 = total scale
  at 25 C, NOT in situ; very early data possibly NBS), DOC, DON, DOP, TDN, TDP, POC (HOT "PC"), PON, POP (nmol/kg),
  bSi (PSi, nmol/kg), N2O, Chl (fluorometric, ug/L), Chl_HPLC + 19 HPLC_* pigments (ng/L), BactAbund/ProAbund/SynAbund/
  PicoEukAbund (1e5 cells/mL). Bottle pCO2 and water-column PIC exist in the file, but no samples pass the flags.
- Depth: CTDPRS (dbar) converted to depth with UNESCO 1983 at 22.75N. Time: cast `time` (UTC, ISO) from the
  header (cast begin time).
- Flags: QUALT1..QUAL10 are packed digit strings, one digit per variable in the order of the BCO-DMO dataset
  description. 1 = not QC'd, 2 = good, 3 = suspect, 4 = bad, 5 = missing, 9 = not measured. **Kept 1 and 2**, dropped
  the rest. CTDTMP/THETA have no flag and are kept if finite. Duplicate/replicate bottles are kept as separate rows.
- Units: almost all umol/kg (gravimetric). Model mmol/m3 needs x rho (~1.022-1.050 kg/L), e.g. via in-situ density
  from T/S. Pigments and fluorometric Chl are per litre.

### 2. raw/bcodmo_pp/ - 14C primary production
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_737163.csv ; file `hot_pp_bcodmo737163.csv`
- Citation (verify author list on landing page): Karl, D. M. et al. Primary productivity measurements from the Hawaii Ocean Time-Series
  (HOT) project. BCO-DMO, Version 1, doi:10.1575/1912/bco-dmo.737163.1 (CC-BY 4.0).
- Span 1988-10 to 2016-10 (cruise 287; v1 content), 5-175 m.
- PP = mean of up to 3 light replicates (mg C m-3 per dawn-to-dusk incubation, ~12-15 h). Dark bottles are NOT subtracted
  (PP_dark is given separately, 1995-2000 only). Integrate over depth for mg C m-2 d-1. Incubation types are I (in situ
  array), O (on-deck), R (rosette-sampled) and N, recorded in qc as `inc=`.
- Flags: a 10-digit string (bottle, chl, pheo, L1-3, D1-3, salinity). A light rep is kept when its digit and the bottle
  digit are 1 or 2. qc = `<L1L2L3 flags>|inc=<type>|n=<reps>`.
- **Time zone:** Date/Start_time are local HST (BCO-DMO's `time` column wrongly carries a 'Z'). The incubation start time
  is converted HST to UTC (+10 h). When start time is missing, 06:00 HST is assumed.
- The PP file's Chl/Pheo/flow-cytometry columns are not used (the bottle file already has them).

### 3. raw/bcodmo_flux/ - sediment traps (floating, 150 m standard; also 70-500 m early)
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_737393.csv ; file `hot_flux_bcodmo737393.csv`
- Citation (verify author list on landing page): Karl, D. M. et al. Sediment trap flux measurements from the Hawaii Ocean Time-Series (HOT)
  project at station ALOHA. BCO-DMO, Version 1, doi:10.1575/1912/bco-dmo.737393.1 (CC-BY 4.0).
- Span 1988-12 to 2016-10. Variables: POC_flux (HOT "Carbon", mg C m-2 d-1), PON_flux, POP_flux (mg P), bSi_flux
  (mg Si), mass_flux, PIC_flux (150 m only, from 2001).
- **No dates in the file**, so time = HOT cruise mid-date (from HOT_surface_CO2.txt, by cruise number). The traps are
  free-drifting, deployed ~2.5-3 d during each cruise. qc = `trt=<C combined|I individual|O|W>|n=<replicates>`. No
  quality flags. Rows are replicate summaries; several treatments per depth can share a time.
- Deep moored traps (2800/4000 m, 1992-2004) are not in this dataset (they were on ftp.soest only).

### 4. raw/hot_co2/ - HOT surface CO2 system product (Dore)
- URL: https://hahana.soest.hawaii.edu/hot/hotco2/HOT_surface_CO2.txt (+ HOT_surface_CO2_readme.pdf)
- Citation (requested): "adapted from: Dore, J.E., R. Lukas, D.W. Sadler, M.J. Church, and D.M. Karl. 2009. Physical and
  biogeochemical modulation of ocean acidification in the central North Pacific. Proc Natl Acad Sci USA 106:12235-12240."
- Span 1988-10 to 2024-12 (cruise 1-355, updated 2026-01-01). One value per cruise = mean 0-30 dbar, stored at
  depth_m = 15. Date = mid-cruise date (00:00 UTC used).
- Variables: temp, salt, PO4, SiO2, DIC (umol/kg), ALK (ueq/kg), pH (pHmeas_insitu: measured, total scale, at in-situ T),
  pCO2 (pCO2calc_insitu, **calculated** from DIC/TA with CO2SYS). -999 = missing. Some early DIC/TA values are filled from
  Keeling's series (see `notes` letters, kept in qc). Missing PO4/SiO2 were set to 0.07/1.04 by the PI.

### 5. raw/ncei_whots_pco2/ - WHOTS mooring MAPCO2 (PMEL), NCEI OCADS accession 0100080
- URL: https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0100080/ (all `WHOTS_158W_23N_*.csv`, QF logs, one PI_OME pdf);
  landing https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/Moorings/WHOTS_158W_23N.html
- Citation: Sutton, A. J., Sabine, C. L., et al. High-resolution ocean and atmosphere pCO2 time-series measurements from
  mooring WHOTS_158W_23N. NCEI Accession 0100080 (exact citation/DOI in the PI_OME/xml metadata files).
  Mooring: WHOTS (WHOI/UH; Weller, Plueddemann, Lukas).
- Span 2007-06 to 2025-08, 3-hourly. Gaps: no files for ~Feb 2009-Jul 2009, Jul 2011-Jun 2012, Nov 2012-Jul 2013 and
  Oct 2020-Aug 2021. Mooring sits ~10 km from ALOHA (22.67N, 157.97W).
- Variables: pCO2 and fCO2 (seawater, uatm, sat. at SST), pCO2_air, temp (SST) and salt (SSS) from the MAPCO2/SBE at
  ~0.5 m (depth_m = 0.5; air = 0), pH (total scale, from 2012), Chl (fluorometer, ug/L, from 2013), O2 (optode umol/kg,
  2013-2024; a few high values pass QF, so use with care).
- Flags: WOCE-style 2 = good, 3 = questionable, 4 = bad. **Kept QF = 2 only** for CO2 (SW QF used for pCO2/fCO2, Air QF
  for pCO2_air), pH, CHL and DOXY. SST/SSS have no QF and are kept if not -999. Times are UTC. The 2007-2010 files have
  CR line endings and m/d/yy dates, handled in the reader.

### 6. raw/scripps_co2/ - Scripps (Keeling) surface DIC/ALK at ALOHA
- URL: https://scrippsco2.ucsd.edu/wp-content/uploads/sites/533/2026/01/HAWI.csv
  (page https://scrippsco2.ucsd.edu/data/seawater_carbon/ocean_time_series.html)
- Citation given in the file: Lueker, T.J., C.D. Keeling, P.R. Guenther, M. Whalen, and W.G. Mook. Inorganic Carbon
  Variations in Surface Ocean Water near Bermuda. UC San Diego, SIO. https://escholarship.org/uc/item/8742p2nb, 1998.
  (This is the generic Scripps citation. For HOT also cite Keeling et al. 2004 GBC 18, GB4006.)
- Span 1988-10 to 2023-10, 1-38 m. Variables: DIC, ALK (umol/kg), temp, salt. Longitude in the file is positive deg W
  and is converted to negative. No flags (already screened by the PI).

### 7. raw/bcodmo_pp_v5/ - 14C primary production, BCO-DMO 737163 v5 (HOT-1 to HOT-355), added 2026-10-08
- URL: https://www.bco-dmo.org/dataset/737163 ; file
  https://datadocs.bco-dmo.org/dataset/737163/file/g7zVP0nulGW5vn/737163_v5_prim_prod_hot001_hot355.csv (563 KB,
  size checked). Landing page saved as `landing_page_737163.html`.
- Citation: Karl, D. M., Fujieki, L. A. (2026). Primary productivity measurements for the Hawaii Ocean Time-series (HOT)
  program from October 1988 to December 2024 at Station ALOHA. BCO-DMO. (Version 5) Version Date 2026-04-10.
  doi:10.26008/1912/bco-dmo.737163.5 (CC-BY 4.0). Plus the HOT-DOGS acknowledgement above.
- **De-duplication:** only cruises absent from `bcodmo_pp` (v1, HOT-1 to 287) are used, i.e. HOT-289 to 355
  (2017-01 to 2024-12). Source tag `bcodmo_pp_v5`. PP = mean of light replicates (mg C m-3 per dawn-dusk incubation,
  dark NOT subtracted), the same definition as v1. No post-2016 dark bottles exist, so there are no PP_dark rows.
- Time: `Start_ISO_DateTime_UTC` (true UTC incubation start in v5; this matches the v1 HST+10 h conversion exactly for
  the shared cruises). lat/lon from the file.
- Flags: v5 has one flag per light set (`Flag_Light`), plus `Flag_Bottle` and `Flag_Dark` (HOT scheme). Kept 1/2 for
  both bottle and light. qc = `<Flag_Light>|inc=<type>|n=<reps>`. Caveat: v1 had per-replicate flags. On the shared
  cruise 287 the v5 means differ from v1 by up to ~15% at some depths because v5 can't drop individual flagged
  replicates.

### 8. raw/bcodmo_flux_v4/ - sediment traps, BCO-DMO 737393 v4 (HOT-2 to HOT-355), added 2026-10-08
- URL: https://www.bco-dmo.org/dataset/737393 ; file
  https://datadocs.bco-dmo.org/dataset/737393/file/XY4AqAYiGJgrXw/737393_v4_hot_particle_flux.csv (99 KB).
- Citation: Karl, D. M., Fujieki, L. A. (2026). Sediment trap flux measurements for the Hawaii Ocean Time-series (HOT)
  project from December 1988 to December 2024 at Station ALOHA. BCO-DMO. (Version 4) Version Date 2026-04-23.
  doi:10.26008/1912/bco-dmo.737393.4 (CC-BY 4.0).
- **De-duplication:** only cruises absent from `bcodmo_flux` (v1) are used: HOT-289 to 355, 150 m only, 2017-01 to
  2024-12. The values for the shared cruise 287 are identical to v1. Same variables and units as v1 (POC_flux = HOT
  "Carbon", mg C m-2 d-1; PON_flux, POP_flux (mg P), bSi_flux (mg Si), mass_flux, PIC_flux). v4 has real deployment
  start/end times (UTC), so time = **deployment midpoint**. v1 rows still use the cruise mid-date.
  qc = `trt=<treatment>|n=<replicates>`. No quality flags.

### 9. raw/bcodmo_spots/ - SPOTS synthesis, ALOHA bottle subset (fills HOT bottle chemistry to Dec 2019), added 2026-10-08
- URL: ERDDAP https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_896862_v2.csv?&TimeSeriesSite=%22ALOHA%22&DATE%3E=20160101
  (query in `query_url.txt`; DAS in `bcodmo_dataset_896862_v2.das`); file `spots_v2_ALOHA_2016on.csv` (6.2 MB).
  Landing page https://www.bco-dmo.org/dataset/896862 .
- Citation: Lange, N., Fiedler, B., Álvarez, M., Benoit-Cattin, A., Benway, H., Buttigieg, P., Coppola, L., Currie, K. I.,
  Flecha, S., Gerlach, D. S., Honda, M. C., Huertas, E. I., Kinkade, D., Muller-Karger, F., Lauvset, S. K., Körtzinger, A.,
  O'Brien, K. M., Ólafsdóttir, S., Pacheco, F. C., Rueda-Roa, D., Skjelvan, I., Wakita, M., White, A. E., Tanhua, T. (2024).
  Synthesis Product for Ocean Time Series (SPOTS). BCO-DMO. (Version 2) Version Date 2024-02-22.
  doi:10.26008/1912/bco-dmo.896862.2 ; paper: Lange et al. (2024) ESSD (essd-2023-238). Also cite the HOT data
  (SPOTS `DOI` column = 10.1575/1912/bco-dmo.3773.1) and the HOT-DOGS acknowledgement.
- Content: the HOT bottle file, reformatted to WHP/exchange-style names with WOCE flags. In SPOTS, `STNNBR` = HOT
  cruise number. ALOHA coverage ends at HOT-317 (2019-12-20). **De-duplication:** only cruises absent from
  `bcodmo_bottle` (3773 v1, ends HOT-288) are used, i.e. HOT-289 to 317. Check: on the shared 2016 cruises SPOTS and
  3773 v1 agree exactly (436 NO3 pairs, 0 difference in value or cast time).
- Variables: CTDTMP->temp, CTDSAL->salt, OXYGEN->O2, CTDOXY->O2_CTD, NITRAT->NO3 (HOT NO3+NO2), PHSPHT->PO4,
  SILCAT->SiO2, TCARBN->DIC, ALKALI->ALK (umol/kg as reported), PH_TOT->pH (total scale at 25 C), DOC (2017 only),
  TPC->POC, TPN->PON, TPP->POP (all umol/kg as reported; note v1 POP is nmol/kg). NITRIT, NH4, POC/PON/POP proper and
  PCO2 are empty for ALOHA. HPLC pigments, fluorometric Chl, bSi, N2O, DON/DOP and cell counts are NOT in SPOTS, so
  they still end in 2016.
- Depth: CTDPRS converted with UNESCO 1983 at 22.75N (as for source 1). Time: DATE + TIME (UTC).
- **Date fix:** SPOTS dates are wrong for the January cruises HOT-289, 299 and 309 (the HOT mmddyy date lost its
  leading zero, e.g. "12317" = 23 Jan 2017 became 2017-12-03). The reader rebuilds them (month 1, day = last digit of the
  wrong month + wrong day) and accepts a fix only if it falls within 10 days of the HOT cruise date in
  HOT_surface_CO2.txt. 658 rows were repaired and none dropped.
- Flags: WOCE `*_FLAG_W`, kept 2 (good) and 6 (replicate mean); temp has no flag.

### 10. raw/bcodmo_hot_metals/ - dissolved trace metals at ALOHA, BCO-DMO 962986 v2 (Dec 2020 to Nov 2023), added 2026-10-08
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_962986_v2.csv (public) ; landing
  https://www.bco-dmo.org/dataset/962986 . File `bcodmo_dataset_962986_v2.csv` (53 KB) + .das + landing page.
- Citation: Hawco, N. J., Bates, E. S. (2025). Water column dissolved and total dissolvable metal concentrations from
  Hawaii Ocean Timeseries (HOT) R/V Kilo Moana cruises at station ALOHA, North Pacific Subtropical Gyre, from December
  2020 to November 2023. BCO-DMO. (Version 2) Version Date 2025-08-29. doi:10.26008/1912/bco-dmo.962986.2
- 21 HOT cruises (HOT-325 onward), 15-900 m, trace-metal clean sampling. Reader `std_tm.py`. Variables (dissolved):
  Fe, Mn, Co (labile dissolved Co, not UV-oxidised), Ni, Cu, Zn, Cd, Pb, Ti, all nmol/L as reported. Total dissolvable
  (td*) columns are not used.
- Flags (GEOTRACES/SeaDataNet): 1 good, 2 probably good, 3 probably bad, 4 bad, 6 below detection, 9 missing. Kept
  1/2 (PI recommendation). Time = BCO-DMO `time` (UTC, converted from HST rosette deployment time).
- This is the only open dFe time series at ALOHA found. The GEOTRACES GP15 ALOHA station (2018) is in IDP2025 only
  (not fetched, see Synthesis > GEOTRACES).

### Not obtained
- HOT CTD 2-dbar profiles (BCO-DMO 3937_v2, 1988-2023, public on ERDDAP). Too large for the slow link (~10^7 rows), so
  skipped. The bottle file already has CTD T/S/O2 at bottle depths (to 2016). A depth-strided subset would extend T/S/O2
  to 2023 if needed.
- MLD: HOT doesn't distribute an MLD time series in these files. Compute it from the CTD profiles if needed.
- HOT bottle chemistry 2020-2024 (HOT-318 to 355): still missing. BCO-DMO 3773 v4 is "Under revision" with no file,
  ERDDAP 3773_v3 is 401 "log in", the WHOAS v3 archive is behind a CAPTCHA and the SOEST FTP/HOT-DOGS hosts time out.
  **Manual step:** open https://www.bco-dmo.org/dataset/3773 in a browser and download the v3/v4 bottle CSV once it is
  listed (or open https://doi.org/10.26008/1912/bco-dmo.3773.3, pass the WHOAS "Security Check" CAPTCHA and download).
  Alternatively use HOT-DOGS bottle extraction (https://hahana.soest.hawaii.edu/hot/hot-dogs/bextraction.html, Station 2,
  all public variables) from a network that can reach scopeserver.soest.hawaii.edu. Put the file in a new dir
  `HOT/raw/bcodmo_bottle_v3/` (the column names will need a small reader, de-duplicated by cruise like SPOTS). HPLC, Chl,
  bSi, N2O, DON/DOP and cell counts after 2016 need the same file.
- WHOTS physical mooring (UOP/OceanSITES) subsurface T/S, and HOT ACO, were not fetched (out of scope / time).


<!-- from BATS/SOURCES_BATS.md -->
## BATS (Bermuda Atlantic Time-series Study) - sources

Model column: 31.67N, 64.17W. Accessed 2026-10-07. Reader: obs/scripts/std_bios.py (site 'BATS').
The old BIOS page (bats.bios.edu/bats-data) now redirects to bios.asu.edu and links all BATS data to
BCO-DMO (project 2124, https://www.bco-dmo.org/project/2124). The BCO-DMO versions are the current releases.

Spatial filter: samples > 100 km from 31.67N 64.17W dropped (BATS casts spread around the BATS box; ~10% of
bottles are 50-100 km away). lat/lon of each sample kept in the tidy file for tighter filtering.
Cruise types kept: BATS Core (monthly) + Bloom A/B (biweekly Feb-Apr). BATS Validation (BVAL) spatial survey
cruises (datasets 917255, 926534, 939210) NOT used.

Flag convention (all BCO-DMO BATS files): per-parameter QF 1=unverified, 2=verified/acceptable,
3=questionable, 4=bad, 9=no data; bottle QF -3=suspect bottle. KEPT: QF in {1,2}, bottle QF != -3,
depth QF not 3/4. Missing (-999/NA/blank) dropped. qc column = the kept QF.

### 1. bcodmo_bottle - BATS discrete bottle file (v10)
- URL: https://www.bco-dmo.org/dataset/3782 ; file https://datadocs.bco-dmo.org/dataset/3782/file/QAD9QkVh9ORw3x/3782_v10_bats_bottle.csv
- Files: raw/bcodmo_bottle/3782_v10_bats_bottle.csv (17.7 MB), bats_bottle_release_v010_update.txt, Dataset_description.pdf
- Citation: Johnson, R. J., Bates, N. R., Lomas, M. W., Smith, D., Lethaby, P. J., Bakker, R., Davey, E., Derbyshire, L.,
  Enright, M., Garley, R., Hayden, M. G., Lomas, D., May, R., Medley, C., Stuart, E., Chambers, E. (2026) Discrete bottle
  samples collected at the Bermuda Atlantic Time-series Study (BATS) site in the Sargasso Sea from October 1988 through
  December 2025. (Version 10) Version Date 2026-07-24. BCO-DMO. doi:10.26008/1912/bco-dmo.3782.10
- Time: ISO_DateTime_UTC (cast time, UTC; also yyyymmdd + hhmm UTC + decimal_year in file). Depth in m.
- Variables -> tidy: Temperature->temp (ITS-90 in-situ, CTD at bottle firing); Salinity (bottle) -> salt, with
  CTD_Salinity filling gaps (bottle salinity only on ~23% of bottles); Oxygen_1->O2; CO2->DIC; Alkalinity->ALK;
  NO3_plus_NO2->NO3 (it is NO3+NO2); NO2; PO4; Silicate->SiO2 (all umol/kg); POC, PON (ug/kg); POP (umol/kg);
  TOC (umol/kg; total organic C, ~DOC); TN->TDN (umol/kg; total dissolved N incl. NO3; all TN data are from cruise
  >=122 so no DON rows emitted); TDP, SRP (low-level) -> TDP, SRP_LL (nmol/kg); Bio_Si->bSi, Litho_Si->lSi (umol/kg);
  Bact_Enum->BactAbund (1e8 cells/kg); Prochlorococcus, Synechococcus, Picoeuk, Nanoeuk (cells/mL).
- Units are per kg: multiply by in-situ density (~1.025-1.05 kg/L) x 1000 to compare with model mmol/m^3
  (umol/kg * rho[kg/m^3] / 1000 = mmol/m^3). POC/PON ug/kg -> divide by 12.011/14.007 for umol/kg.
- Span 1988-10 to 2025-12; 0-4730 m. Exact duplicate rows removed by standardizer.

### 2. bcodmo_pp - BATS 14C primary production (v8)
- URL: https://www.bco-dmo.org/dataset/893182 ; file 893182_v8_bats_primary_production.csv (1.0 MB) + Dataset_description.pdf
- Citation: Johnson, R. J., Bates, N. R., Lethaby, P. J., Smith, D., Medley, C., Stuart, E., May, R. (2026) Primary
  productivity estimates from the incubation of seawater collected at the BATS site from December 1988 through
  December 2025. (Version 8) Version Date 2026-07-31. BCO-DMO. doi:10.26008/1912/bco-dmo.893182.8
- pp (mean light - dark, mgC/m^3/day) -> PP; QF10_pp in {1,2} kept. Time = in-situ array deployment time (UTC) when
  given, else CTD cast time, else date. Depth 0-160 m (8 depths, 14C in-situ incubation dawn-dusk).
- Caveat: post-2013 array longitudes are reported positive in the source (sign error); reader forces W.

### 3. bcodmo_pigments - HPLC + fluorometric pigments (v10)
- URL: https://www.bco-dmo.org/dataset/893521 ; file 893521_v10_bats_pigments.csv (1.4 MB) + Dataset_description.pdf
- Citation: Bates, N. R., Johnson, R. J., Lethaby, P. J., Medley, C., Smith, D., Stuart, E., May, R., Derbyshire, L.
  (2026) HPLC and fluorometric derived phytoplankton pigment concentrations from seawater collected at the BATS site
  from October 1988 through December 2025. (Version 10) Version Date 2026-07-31. BCO-DMO. doi:10.26008/1912/bco-dmo.893521.10
- p16_Chl -> Chl (fluorometric chl a, ug/kg); p17_Phae -> Phaeo (ug/kg); p14 -> Chl_HPLC (HPLC chl a, MV+DV
  combined, ng/kg); other HPLC pigments -> HPLC_<name> (ng/kg): chl_c3, chlide_a, chl_c1c2, perid, but_fucox, fucox,
  hex_fucox, prasino, diadino, allo, diato, zea_lut, chl_b, ab_carot, lut, zea, a_carot, b_carot.
- NOTE units: Chl in ug/kg but Chl_HPLC in ng/kg (x1e-3 to compare). 0-500 m.

### 4. bcodmo_flux - BATS PITS sediment traps (v8)
- URL: https://www.bco-dmo.org/dataset/894099 ; file 894099_v8_bats_particle_flux.csv (264 KB) + Dataset_description.pdf
- Citation: Johnson, R. J., Bates, N. R., Lomas, M. W., Steinberg, D. K., Derbyshire, L., Hayden, M. G., Lomas, D.,
  Lethaby, P. J., Lopez, P. Z., May, R., Smith, D., Stuart, E., Enright, M. (2026) Determination of carbon, nitrogen, and
  phosphorus content in sinking particles at the BATS site from December 1988 to December 2025 using a Particle
  Interceptor Trap System (PITS). (Version 8) Version Date 2026-07-30. BCO-DMO. doi:10.26008/1912/bco-dmo.894099.8
- Surface-tethered PITS, ~3-day deployments at 150, 200, 300 (and a few 400) m. Time = deployment midpoint (UTC).
  M_avg->mass_flux (mg/m2/d), C_avg->POC_flux (mgC/m2/d), N_avg->PON_flux (mgN/m2/d), P_avg->POP_flux (mmolP/m2/d).
  Averages of replicate tubes, NOT field-blank corrected (FBC_avg/FBN_avg exist in raw for recent years only).
  No QC flags in file. Model flux comparisons: mgC -> mmolC divide by 12.011.

### 5. ncei_ofp - Oceanic Flux Program deep moored traps (v2, 1978-2018)
- URL: https://www.ncei.noaa.gov/archive/accession/0291437 (archived copy of BCO-DMO dataset 704722 v2);
  files from https://www.ncei.noaa.gov/data/oceans/archive/arc0227/0291437/1.1/data/0-data/
- Files: raw/ncei_ofp/dataset-704722_ofp-primary-particle-flux__v2.tsv (405 KB), *_README.txt, Dataset_description.pdf
- Citation: Conte, M. H. (2019) Primary particle flux data (500, 1500, and 3200m depths) of the OFP sediment trap
  time-series in the northern Sargasso Sea from 1978-2018 (version 2, 2019-12-13). BCO-DMO dataset 704722 / NCEI
  Accession 0291437. Also cite Conte, M. H., Ralph, N., & Ross, E. H. (2001) DSR II 48, 1471-1505.
- **PI request: "Please contact the PI (Maureen Conte, BIOS) prior to any use of these data."** Data are open
  (no login) but contact before publication use.
- Moored traps at 500, 1500, 3200 m near 31.83N 64.17W (~15 km from BATS model point). Time = sample start date +
  duration/2 (durations ~2 weeks-2 months). MassFlux->mass_flux, CorgFlux->POC_flux (mgC/m2/d), Nflux->PON_flux
  (mgN/m2/d), CarbFlux->PIC_flux (mg CaCO3/m2/d; x0.12 for mgC), PtotalFlux->POP_flux (mgP/m2/d, total P),
  OpalFlux->bSi_flux (mg opal/m2/d), LithFlux->lith_flux. 'nd' = not determined, dropped. 500/1500 m fluxes are
  semi-quantitative Oct 1992-Jan 1996 (per PI). No flags.
- Newer version (v4, 1978-2024) listed at https://www.bco-dmo.org/dataset/704722 but "Under revision" with no
  downloadable file on 2026-10-07. Bulk-composition/element flux companion (NCEI 0278626, 2000-2015) not fetched.

### 6. ncei_btm_mooring - Bermuda Testbed Mooring MAPCO2 (2005-2007)
- URL: https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/Moorings/BTM_64W_32N.html ;
  data https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0100065/ (NCEI accession 0100065)
- Files: BTM_64W_32N_Oct05_Jul06.csv, BTM_64W_32N_Jul06_Mar07.csv, BTM_64W_32N_Mar07_Oct07.csv, README (630 KB total)
- Citation (from README): Sabine, C., N. Bates, S. Maenner, R. Bott, and A. Sutton. 2010. High-resolution ocean and
  atmosphere pCO2 time-series measurements from mooring BTM_64W_32N. CDIAC, ORNL. doi:10.3334/CDIAC/otg.TSM_BTM_64W_32N
- 3-hourly, 31.78N 64.2W, UTC. fCO2_SW_sat -> fCO2 (uatm), fCO2_Air_sat -> pCO2_air (uatm), SST -> temp, SSS -> salt.
  Depth set to 0.5 m (surface buoy). Flag: xCO2 QF 2 = good kept (3/4 dropped); SST/SSS unflagged.
- Only 2005-10 to 2007-10 exists for BTM pCO2.

### 7. scripps_co2 - Scripps CO2 program surface DIC/ALK at BATS (BATS.csv)
- URL: https://scrippsco2.ucsd.edu/data/seawater-carbon-data/ocean-time-series-data/ ;
  file https://scrippsco2.ucsd.edu/wp-content/uploads/sites/533/2026/01/BATS.csv (created 12/12/2024)
- Citation (from file): Lueker, T.J., C.D. Keeling, P.R. Guenther, M. Whalen, and W.G. Mook. Inorganic Carbon
  Variations in Surface Ocean Water near Bermuda. UC San Diego: SIO. https://escholarship.org/uc/item/8742p2nb, 1998.
  Acknowledge HOT/BATS staff for sample collection (per page); queries to Ralph Keeling.
- Independent (Keeling/Dickson lab) surface DIC, ALK (umol/kg), d13C-DIC (permil) 1989-2024, 1-21 m. Date only
  (time set 00:00 UTC). Longitudes given as positive degW in file (reader negates). Pre-screened, no flags.

### 8. bcodmo_bait_fe - BAIT 2019 dissolved Fe + d56Fe (Conway lab), BCO-DMO 936824 v1, added 2026-10-08
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_936824_v1.csv ; landing https://www.bco-dmo.org/dataset/936824
- Citation: Conway, T. M., Boiteau, R. M., Toth, E. (2024). Dissolved iron concentrations and stable isotope ratios from
  water column samples collected during four Bermuda Atlantic Iron Time-series (BAIT) cruises EN631, AE1909, AE1921,
  AE1930 in the Western Subtropical North Atlantic Gyre in 2019. BCO-DMO. (Version 1) Version Date 2024-09-17.
  doi:10.26008/1912/bco-dmo.936824.1 . BAIT is GEOTRACES process study GApr13.
- Fe = Fe_D_CONC_BOTTLE (GO-Flo, nmol/kg) and Fe_D_CONC_BOAT_PUMP (surface pump, nmol/kg; see units column), Mar,
  May, Aug and Nov 2019, 0.5-1678 m. Time = sampling date (00:00 UTC). d56Fe not used. Flags SeaDataNet (data the PI
  considers accurate = 2, 9 = missing): kept 1/2. 100 km BATS radius applied (all samples pass). Reader `std_tm.py`.

### 9. bcodmo_bait_tm_bottle - BAIT 2019 trace-metal rosette bottle data (Sedwick lab), BCO-DMO 937302 v1, added 2026-10-08
- URL: https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_937302_v1.csv ; landing https://www.bco-dmo.org/dataset/937302
- Citation: Sedwick, P. N., Sohst, B., Johnson, R. J., Williams, T. E. (2024). Concentrations of trace metals and dissolved
  macronutrients and CTD sensor data from four cruises in the Bermuda Atlantic Time-series Study (BATS) region in March,
  May, August and November 2019. BCO-DMO. (Version 1) Version Date 2024-09-26. doi:10.26008/1912/bco-dmo.937302.1
- Fe (DFe), Mn (DMn), Al (DAl), nmol/L as reported, same BAIT casts as source 8 (an independent lab on the same
  samples, so **two Fe estimates per bottle**: compare or pick one, do not average blindly). Time = hydrocast UTC.
  Flags: 1 good, 2 likely contaminated/questionable, 3 questionable (Sc internal standard), 4 not determined, 5 below
  detection; **kept 1 only**. Soluble Fe/Mn, macronutrients and CTD sensor columns not used.
- GA03 (KN204-01, 2011) BATS station dissolved Fe is in GEOTRACES IDP2025 only (not fetched, see Synthesis > GEOTRACES).
  CCHDO's GA03 bottle file (316N20111106_gt_hy1.csv) was checked and has hydrography only, no Fe.

### Not fetched / notes
- BATS CTD 2-dbar profiles (https://www.bco-dmo.org/dataset/3918, 659 MB CSV; T, S, O2, fluorescence, PAR, beam):
  too large for budget; can be subset via ERDDAP https://erddap.bco-dmo.org/erddap/tabledap/ (see dataset page)
  if continuous T/S/O2/fluorescence profiles or MLD are needed. No MLD product is provided by BATS.
- No surface pH time series at BATS in these files (pH could be computed from DIC/ALK with CO2SYS).
- No login-required sources encountered.


<!-- from HydroS/SOURCES_HydroS.md -->
## Hydrostation S (Bermuda) - sources

Model column: 32.17N, 64.5W. Accessed 2026-10-07. Reader: obs/scripts/std_bios.py (site 'HydroS').
Spatial filter: samples > 50 km from 32.17N 64.5W dropped (~1.5% of bottles); rows with missing position kept.

### 1. bcodmo_bottle - Hydrostation S discrete bottle file (v8)
- URL: https://www.bco-dmo.org/dataset/859990 ;
  file https://datadocs.bco-dmo.org/dataset/859990/file/XYNNL2zCqnNXQw/859990_v8_hydrostation_s_bottle.csv
- Files: raw/bcodmo_bottle/859990_v8_hydrostation_s_bottle.csv (6.0 MB), Dataset_description.pdf
- Citation: Bates, N. R., Johnson, R. J., Lethaby, P. J., Smith, D., Medley, C., May, R., Derbyshire, L., Goncalves Neto,
  A., Bakker, R., Stuart, E., Chambers, E. (2026) Discrete bottle data from Hydrostation S in the Sargasso Sea from
  January 1955 through December 2025. (Version 8) Version Date 2026-08-19. BCO-DMO. doi:10.26008/1912/bco-dmo.859990.8
- Variables: Temperature -> temp (in-situ; reversing thermometers before the CTD era), Salinity_1 (bottle, PSS-78)
  -> salt with CTD_Salinity filling gaps, Oxygen -> O2 (umol/kg). No nutrients/carbon in this file.
- Time: ISO_DateTime_UTC (UTC). ~Biweekly 1955-2025, 0-4200 m.
- Flags: 1=unverified, 2=verified, 3=questionable, 4=bad, 9=no data; bottle QF -3=suspect. Kept QF in {1,2},
  bottle != -3, depth QF not 3/4.
- O2 in umol/kg: x rho/1000 for mmol/m^3.

### 2. scripps_co2 - Scripps CO2 program surface DIC/ALK at Hydrostation S (BERM.csv)
- URL: https://scrippsco2.ucsd.edu/data/seawater-carbon-data/ocean-time-series-data/ ;
  file https://scrippsco2.ucsd.edu/wp-content/uploads/sites/533/2026/01/BERM.csv (created 12/12/2024)
- Citation (from file): Lueker, T.J., C.D. Keeling, P.R. Guenther, M. Whalen, and W.G. Mook. Inorganic Carbon
  Variations in Surface Ocean Water near Bermuda. UC San Diego: SIO. https://escholarship.org/uc/item/8742p2nb, 1998.
  Acknowledge HOT/BATS staff for sample collection (per page); queries to Ralph Keeling.
- Surface (0-25 m) DIC, ALK (umol/kg), d13C_DIC (permil), 1983-09 to 2024-02 (~monthly). Date only, time set
  00:00 UTC; a few dates in mm/dd/yy format handled. Longitude positive = degW in file (reader negates). No flags
  (provider-screened).
- This is the only open carbonate record found for Hydro S; the Bates Hydro S/BATS DIC record is in the BATS bottle file
  (BATS site), not a separate Hydro S BCO-DMO dataset.

### Not fetched
- Hydrostation S CTD 2-dbar profiles: https://www.bco-dmo.org/dataset/860014 (large; ERDDAP subset possible).
- No nutrients, chl, PP, pH or MLD products for Hydro S. No login-required sources encountered.


<!-- from Papa/SOURCES_Papa.md -->
## Ocean Station Papa / Line P P26 (model column 50.1N, 144.9W): source notes

All files accessed 2026-10-07. Reader: `obs/scripts/std_papa.py`. Output: `obs/Papa/Papa_obs.csv`.

### 1. DFO-IOS bottle profiles (Line P / Station P), `raw/ios_bot_erddap/`
- URL: CIOOS Pacific ERDDAP, dataset `IOS_BOT_Profiles`
  https://data.cioospacific.ca/erddap/tabledap/IOS_BOT_Profiles.html
  Subset: 49.9-50.3N, 145.2-144.6W (P26 nominal 50.0N 145.0W).
- Files: `IOS_BOT_Profiles_P26.csv` (row 2 = units), `IOS_BOT_Profiles_metadata.csv`.
- Variables: temp (in-situ; CTD TEMPS901/902/601/ST01, reversing thermometer TEMPRTN1, or the
  harmonized sea_water_temperature), salt (PSS-78; bottle or CTD), O2 (umol/kg DOXMZZ01; mL/L
  DOXYZZ01 kept *only* where umol/kg is absent, units column says which), NO3 (NTRZAAZ1 =
  nitrate+nitrite, umol/L), SiO2 (umol/L), PO4 (umol/L), Chl (CPHLFLP1, extracted fluorometric,
  mg/m^3).
- Span 1956-2024 (Weather-ship era 1956-1981 very dense; Line P after 1981 typically
  2-3 cruises/yr: Feb, May/Jun, Aug/Sep). Depths 0-~4300 m.
- QC: the ERDDAP dataset does not expose the IOS flag channels; values are used as archived by
  IOS (qc column empty). Expect occasional outliers.
- Citation/ack: Institute of Ocean Sciences, Fisheries and Oceans Canada (DFO); served by CIOOS
  Pacific (Hakai Institute acknowledgement). License: free use/redistribution, no warranty.
  Contact DFO.PAC.SCI.IOSData-DonneesISO.SCI.PAC.MPO@dfo-mpo.gc.ca.
- Login needed: the waterproperties.ca Line P cruise index (https://www.waterproperties.ca/linep/cruises.php)
  requires an account; individual cruise pages (e.g. /linep/2011-27/index.php) are open but
  one-file-per-cast. Not used.

### 2. Line P carbonate system (NCEI OCADS), `raw/ncei_linep_carbon/`
- URLs: https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/RepeatSections/clivar_line_p.html
  - 0234342: https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0234342/ -> `LineP_for_Data_Synthesis_1990-2019_v1.csv`, `0234342_metadata.html`
  - 0300980 (cruise 2021-008): `0300980_2021-008_data.xlsx`
  - 0310718 (cruise 2024-002): `0310718_18DD20240124_Data.xlsx`
  - Not fetched: 0110260 (per-cruise WHP files 1985-2017, overlaps 0234342), 0302482 (2008).
- Variables taken: DIC, ALK (umol/kg). Nutrients/O2 in these files duplicate the IOS bottles
  and are NOT used (avoids double counting).
- Same 49.9-50.3N, 145.2-144.6W box. Kept data span 1990-2021 (DIC n=1291, ALK n=885); the
  2024-002 P26 DIC/TA are all flagged 3 and dropped.
- Flags: WOCE bottle flags; kept 2 (good) and 6 (mean of replicates); dropped 3, 4, 9 and -999.
  Depth -999 replaced by CTD pressure (dbar ~ m).
- Citation: Franco, A.C.; Ianson, D.; Ross, T.; Hamme, R.C.; Monahan, A.H.; Christian, J.R.;
  Davelaar, M.; Johnson, W.K.; Miller, L.A.; Robert, M.; Tortell, P.D. (2021). A compilation of
  inorganic carbon system and other hydrographic and chemical discrete profile measurements
  obtained during the fifty five Line P cruises in the Northeast Pacific Ocean over the period
  from 1990 to 2019 (NCEI Accession 0234342). NOAA NCEI. https://doi.org/10.25921/zrw8-kn24
  For 0300980 / 0310718 cite the NCEI accession (DFO IOS, D. Ianson et al.).

### 3. NOAA PMEL Papa CO2/OA surface mooring, `raw/pmel_co2_mooring/`
- URL: https://data.pmel.noaa.gov/pmel/erddap/tabledap/pmel_co2_moorings_cba8_5413_09f9.html
  (same data as NCEI OCADS moored time series). Files: `papa_co2_mooring.csv`, `papa_co2_mooring_metadata.csv`.
- 3-hourly, 2007-06-08 to 2025-05-30, at 50.13N 144.84W, sensors at 0.5 m.
- Variables: temp (SST), salt (SSS), pCO2 (seawater, uatm), pCO2_air (uatm; depth set to 0),
  pH (total scale), O2 (umol/kg, salinity compensated; ERDDAP unit attribute "PSU" is a
  metadata error), Chl (nighttime fluorescence, ug/L, calibration bias of 2 applied, Roesler et
  al. 2017). xCO2_air and NTU not used. O2 only 2021-2025, Chl 2010-2025, pH/pCO2 2007-2025.
- QC: published file contains only final-QC'd data (bad values already removed); NaN dropped.
- Citation: Sutton, A.J., et al. (2019) Autonomous seawater pCO2 and pH time series from 40
  surface buoys and the emergence of anthropogenic trends, ESSD 11, 421-439,
  https://doi.org/10.5194/essd-11-421-2019. PMEL asks to be informed of publication use and
  manuscripts sent to PMEL for review; co-authorship may be appropriate.

### 4. NOAA PMEL OCS Papa daily gridded T/S (OceanSITES), `raw/oceansites_pmel_daily/`
- URL: OceanSITES GDAC DATA_GRIDDED/PAPA, file `OS_PAPA_200706_M_TSVM_50N145W_dy.nc`
  (ftp://ftp.ifremer.fr/ifremer/oceansites/DATA_GRIDDED/PAPA/ ; mirror
  https://dods.ndbc.noaa.gov/thredds/catalog/oceansites/DATA_GRIDDED/PAPA/catalog.html - the NDBC
  download was very slow and truncated, so the Ifremer copy was used).
- Daily mean temp at 30 levels (1-300 m) and salt at 28 levels (1-300 m), 2007-06 onward.
- QC: OceanSITES <VAR>_QC; kept 1 (good), 2 (probably good); dropped others.
- Citation: "These data were collected and made available by the Ocean Climate Station Project
  Office of NOAA/PMEL." (PMEL OCS acknowledgement; https://www.pmel.noaa.gov/ocs/Papa)
- Higher-rate (10-min/hourly) per-deployment files exist in DATA/PAPA (~20-60 MB each) and the
  PMEL OCS ERDDAP (papa_hourly_temp/psal); not downloaded (daily is enough for a daily model).
- No MLD product, no nitrate sensor on the PMEL Papa surface mooring.

### 5. OOI Global Station Papa flanking moorings A/B (riser BGC), `raw/ooi_papa_flanking/`
- URL: OOI ERDDAP https://erddap.dataexplorer.oceanobservatories.org/erddap/ datasets
  ooi-gp03flm{a,b}-ris01-03-dostad000 (O2), -04-phsenf000 (pH), -05-flortd000 (chl).
- Daily means computed server-side (orderByMean "time/1day") of samples with QARTOD
  qc_agg <= 2 (1 PASS, 2 NOT_EVALUATED); 3 SUSPECT / 4 FAIL / 9 MISSING excluded.
- FLMA ~49.98N 144.25W, FLMB ~50.33N 144.40W (~45 km E/NE of P26). Depth: daily mean
  instrument pressure from the co-mounted dostad (dbar ~ m; mostly 25-40 m, nominal 30 m,
  mooring blow-down to ~65 m); 30 m assumed where pressure missing.
- Span 2013-07 to 2025-10 (gaps between deployments).
- Variables: O2 (umol/kg), pH (total scale), Chl (ug/L, factory calibration, no bias correction).
- Citation: Ocean Observatories Initiative, funded by the NSF; "OOI data were obtained from the
  NSF Ocean Observatories Initiative Data Portal, http://ooinet.oceanobservatories.org".
- Not fetched: OOI hybrid profiler mooring GP02HYPM (wire-following profiler 150-2400 m with O2,
  T, S) and glider data - large; could be added via the same ERDDAP with daily orderByMean.

### 6. Ocean Station P sediment traps (Wong et al. 1999), `raw/pangaea_osp_trap/`
- URL: https://doi.org/10.1594/PANGAEA.92552 (file `PANGAEA_92552.tab`), CC-BY-3.0.
- Annual-mean fluxes at 200, 1000, 3800 m, 1982-1993: mass_flux (g/m2/yr), POC_flux,
  bSi_flux (opal), PIC_flux (mol/m2/yr), PON_flux (g/m2/yr). time = start of averaging period;
  qc column carries the averaging duration. Multi-year composite rows (1834-4253 d) dropped.
  Only 1989-1993 overlaps the model period.
- Citation: Wong, C.S. (2003): Carbon, Nitrogen and Silica Particle Flux of OSP_trap. PANGAEA,
  https://doi.org/10.1594/PANGAEA.92552; Wong, C.S. et al. (1999) Deep-Sea Res. II 46,
  2735-2760.
- Not found openly: Timothy et al. (2013) 1982-2006 trap series (published tables only);
  newer Line P/Papa trap work (e.g. EXPORTS 2018, BCO-DMO) not fetched.

### Caveats
- Units: nutrients from IOS bottles are umol/L (no conversion needed vs model mmol/m^3 ~ umol/L);
  O2, DIC, ALK in umol/kg need x rho/1000 (~1.026) to compare with model mmol/m^3. Some old O2
  is mL/L (x 44.66 -> umol/L).
- NO3 at Line P = nitrate + nitrite.
- All times UTC (IOS/NCEI/PMEL/OOI report UTC).
- Line P P26 occupied ~2-3x per year since 1981 (seasonal aliasing); weather-ship era pre-1981.
- Mooring products are point measurements near-surface (0.5 m) or ~30-40 m (OOI), not profiles.


<!-- from PAP/SOURCES_PAP.md -->
## PAP (Porcupine Abyssal Plain Sustained Observatory) – sources

Model column: 49.0N, 16.5W. PAP-SO mooring nominal position ~49.0N, 16.3–16.5W (water depth ~4850 m).
All files accessed 2026-10-07. Reader: `obs/scripts/std_pap.py` (called by `standardize_obs.py`).
Tidy output: `obs/PAP/PAP_obs.csv`. All times UTC.

### 1. OceanSITES PAP-1 / PAP-2 / PAP-3 mooring NetCDF (primary)
- `raw/oceansites_gdac_ifremer/` – OceanSITES GDAC, ftp://ftp.ifremer.fr/ifremer/oceansites/DATA/PAP/ (52 files, ~13 MB; 17 files re-issued Nov 2023).
  (https://data-oceansites.ifremer.fr did not respond on 2026-10-07; the FTP did, but was very slow and gave truncated
  transfers – every file was size-checked against the FTP listing.)
- `raw/oceansites_ndbc/` – same file set from the NDBC mirror,
  https://dods.ndbc.noaa.gov/thredds/catalog/oceansites/DATA/PAP/catalog.html . Used only as a fallback when an
  Ifremer copy can't be read. Note: the NDBC copies of the three PAP-2 files are truncated/older.
- Files: `OS_PAP-1_<YYYYMM>_<D|R|P>_<product>.nc` (product = CTD, CTDO, Chl, Wetlabs_Chl_30m, Cyclops_Chl_30m,
  O, PCO2, PCO2_1m, PCO2_30m, ISUS_N_30m), `OS_PAP-2_*_D_CTD.nc` (2002–2004 T/S chain 10–1000 m),
  `OS_PAP-3_201205_P_deepTS.nc` (4850 m T/S).
- Variables used: TEMP->temp (in-situ, degC), PSAL->salt, DOXY/DOXM->O2 (umol/kg), PCO2XXXX/PCO2->pCO2,
  CPHLPS01/CPHLPM01/Chl->Chl (fluorometer chl-a, ug/L, Wetlabs and Cyclops sensors both kept). Only a file's own
  product variables are taken (PCO2/Chl files also carry ancillary CTD copies, which are skipped to avoid duplicates).
  One file per deployment+product: delayed-mode (D) > real-time (R) > provisional (P).
- Positions: early files store longitude as positive (16.4 meaning 16.4W), so it is negated; the 200307 CTD file has a
  bad longitude (−0.69), so lat/lon are left blank for it (it is at the mooring).
- Time span 2002-10 to 2019-07. **No OceanSITES PAP files after the Jul 2019 deployment exist on either GDAC.**
- Sensor depths: until 2007 a CTD chain (PAP-1 30–1000 m; PAP-2 10–1000 m); from 2010 sensors sit at **~1 m
  (buoy keel)** and on a **~30 m frame** (measured pressure ~30–37 m because of mooring knockdown); 1000 m
  microcat in 2007/2009; 4850 m deep T/S 2012. `depth_m` = measured pressure (dbar ≈ m) where valid and QC'd,
  otherwise the nominal DEPTH. 2002–2004 PAP-1 sensors were knocked down well below nominal (e.g. 40 m nominal,
  ~76 m median pressure).
- QC: kept OceanSITES QC 0 (no QC performed), 1 (good), 2 (probably good); dropped 3, 4, 5, 8, 9 and fills.
  Many 2014–2019 files only carry QC=0, so gross range checks were added: temp −2–35, salt 30–40, O2 50–450 umol/kg,
  Chl 0–30 ug/L, pCO2 150–700 uatm (15–71 Pa for the Pa-labelled 2017 file).
- Unit caveats: pCO2 labelled "microPa" (2014, 2015) and "Pa" (2019) has uatm magnitudes (~300–400), so units are
  written as `uatm`; the 2010–2013 "microAtmospheres" files are written as `uatm` too. The 2017 PCO2 file (true Pa)
  has almost no unmasked values, so in practice no 2017–18 pCO2. The 2014 CTDO file's DOXY is a copy of pressure
  (units "decibar") and is skipped; the 2015 CTDO DOXM is all NaN.
- **Nitrate: OceanSITES `OS_PAP-1_201507_R_ISUS_N_30m.nc` holds only CTD variables, so no NO3 time series.**
  pH: none in the OceanSITES files.
- Citation/acknowledgement (from the file metadata): "These data were collected and made freely available from the
  Porcupine Abyssal Plain (PAP) Observatory and the UK national programs that contribute to it, together with
  European Projects (FixO3, EuroSITES, MERSEA, ANIMATE) that have supported it, and the OceanSITES project."
  PIs: S. Hartman, R. Lampitt (NOC). Also cite the OceanSITES GDAC (OceanSITES, 2026, https://doi.org/10.17882/49516 ).
- Suggested reference: Hartman, S.E. et al. (2012), The Porcupine Abyssal Plain fixed-point sustained observatory
  (PAP-SO): variations and trends from the Northeast Atlantic fixed-point time-series, ICES J. Mar. Sci. 69(5), 776–783
  (check the citation before using it).

### 2. NCEI OCADS PAP-SO mooring pCO2 (accession 0312034)
- https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0312034/ ; landing page
  https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/Coastal/PAP.html
- `raw/ncei_ocads_mooring/747F20130425.csv`, `747F20150702.csv`, `Metadata.csv`.
- ProOceanus CO2-Pro at 1 m: pCO2 (uatm), 2013-04-25..2013-12 and 2015-07-02..2016-03-01, 2-6 values per day.
  QC: 2013 file flag 1 = good kept (5 bad / 9 missing dropped); 2015 file has no flags (range check only).
- Mixed date formats (dd/mm/yyyy and dd/mm/yy) parsed day-first. T/S from these files not used (same microcat as the
  OceanSITES 1 m CTD); the 2015 `sbo_37_ox` oxygen has no units (probably ml/L), so it isn't used.
- Overlap: this is the archived version of the OceanSITES 1 m pCO2 for 2013 and 2015-16; OceanSITES 1 m pCO2
  rows inside those windows are dropped in favour of NCEI.
- Citation: Hartman, S.; Lampitt, R. Surface pCO2 measurements from the PAP-SO mooring (NCEI Accession 0312034).
  NOAA National Centers for Environmental Information. Dataset. Accessed 2026-10-07.

### 3. NCEI OCADS RRS James Cook JC278 underway pCO2 (accession 0312066, EXPOCODE 740H20250530)
- https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0312066/ ; `raw/ncei_ocads_jc278/740H20250530.csv` (+ PI_OME.xml).
- Underway SST, SSS, pCO2/fCO2 at SST (wet, uatm), 2025-05-30..06-22. Only points **within 100 km of 49N 16.5W**
  are kept (2025-06-04..06-21 on station at PAP-SO), intake depth 5 m. Flags (SOCAT/WOCE): kept 2 = good only.
  Atmospheric xCO2 columns are empty (NaN) in this file, so there are no pCO2_air rows.
- Citation: Flohr, A.; Hartman, S. Surface underway pCO2 from RRS James Cook cruise JC278 (EXPOCODE 740H20250530)
  (NCEI Accession 0312066). NOAA NCEI. Accessed 2026-10-07. These data will also appear in SOCAT (another fork
  handles SOCAT; watch for duplicates).

### 4. PANGAEA – PAP sediment traps (Lampitt et al.)
- `raw/pangaea_traps/PANGAEA.<id>.tab` from https://doi.pangaea.de/10.1594/PANGAEA.<id>?format=textfile
  - 108311–108319: Lampitt et al. (2001), traps PAP-I, III, V, XV, XVIII, XIX, XX, XXIIIa, XXV (1989–1999);
    parent https://doi.org/10.1594/PANGAEA.724295
  - 794531–794535: Lampitt et al. (2010/2012), traps PAP-XXVI, XXVII, XXVIII, XXXI, XXXIV (1999–2005), 3000 m.
  - 807946: Torres-Valdés et al. (2013) EURO-BASIN Atlantic trap compilation. Downloaded but **not used** (its PAP
    entries duplicate the Lampitt sets above).
- The 1989–90 deployments (PAP-I/III/V, 108311–108313) were at ~47.8N 19.5W (~260 km from the model column), so the
  reader leaves them out (kept in raw/). Used: 1995–2005 deployments at 48.99–49.07N 16.2–16.4W.
- Depths 1000, 3000 and ~4700 m (100 m above bottom); mostly 3000 m after 1998. Variables: mass_flux (mg/m2/d),
  POC_flux (mmol C/m2/d in the 2001 sets ["TOC flux", organic C after carbonate removal]; mg C/m2/d in the 2010/12 sets),
  PIC_flux (mmol or mg C/m2/d), bSi_flux (mg opal (SiO2)/m2/d), PON_flux (mg N/m2/d; total N in the 2001 sets).
  **Units differ between sets: check the units column.**
- `time` = midpoint of each cup's collection interval (start/end dates are in the raw files). No flags.
- No open trap data after 2005 were found (later PAP trap fluxes sit in BODC, see below). License CC-BY-3.0.
- Citation: Lampitt, R.S. et al. (2001) Particle flux from sediment traps PAP-… [datasets], PANGAEA,
  https://doi.org/10.1594/PANGAEA.724295 ; Lampitt, R.S., Salters, V.J.M., de Cuevas, B. et al. (2010/2012) Particle flux
  from sediment trap PAP-XXVI…XXXIV, PANGAEA, doi:10.1594/PANGAEA.794531–794535. Supplement to Lampitt et al. (2001)
  Prog. Oceanogr. 50, 27–63 and Lampitt et al. (2010) Deep-Sea Res. II 57, 1346–1361.

### 5. PANGAEA – JC087 bottle nutrients/O2 at PAP (Jun 2013)
- https://doi.org/10.1594/PANGAEA.832864 ; `raw/pangaea_jc087/PANGAEA.832864.tab`
- CTD bottles 2–4800 m, 2013-06-03..06-14, ~48.5–48.7N 16.0–16.5W. Variables: NO3 (**NO3+NO2**), NO2, NH4, PO4,
  SiO2 (umol/L), O2 (Winkler, umol/L, mean of replicates). Values below detection ("<0.01") are dropped.
- Citation: Stinchcombe, M.C.; Davey, E. (2014): Nutrients and oxygen concentrations measured in water samples from the
  North Atlantic during the James Cook cruise JC087 in March 2013. PANGAEA, https://doi.org/10.1594/PANGAEA.832864
  (the title says March, but the sample dates are June 2013). CC-BY-3.0.

### 6. BODC NODB – PAP cruise bottle data, 2013-2019 (source tag `bodc_nodb_<cruise>`; added 2026-10-08)
- BODC NODB search (https://www.bodc.ac.uk/data/bodc_database/nodb/search/), Site = Porcupine Abyssal Plain (PAP),
  541 series in the CLASS, Porcupine Abyssal Plain Observatory and RAGNARoCC collections, all Unrestricted. Requested
  the 209 water-column-chemistry, water-sample, CTD, fluorescence and PAR series from 1992 on (request RN-2506,
  2026-10-08, after Dustin logged in and signed the BODC licence). ODV text plus per-series HTML docs in
  `raw/bodc_nodb/unz/` (CF-netCDF copies kept; QXF binaries deleted). Reader: `scripts/std_bodc.py`.
- Ingested (bottle data only, SeaDataNet QV 0/1/2 or blank): NO3 (**NO3+NO2**), NO2, PO4, SiO2 (umol/L; DY032 2015,
  JC165 2018, DY103 2019), Winkler O2 (umol/L; DY032, DY050, DY077, DY103), DIC and ALK (umol/kg; JC165 2018,
  DY103 2019), extracted GF/F chl (JC085, JC087 2013). 1878 rows. JC087 nutrients are skipped because they duplicate
  PANGAEA.832864 (section 5).
- Not ingested: moored/CTD sensor series (pCO2 2007-2018, optode O2 2009-2021, ISUS/SUNA NO3 2009 and 2016-2017,
  fluorescence, T/S, PAR, currents), which overlap the OceanSITES files in section 1. This pull has **no mooring series
  after 2019**.
- Citation: data supplied by the British Oceanographic Data Centre (BODC), National Oceanography Centre, UK; cite the
  originating cruise/PI per series as given in the HTML documentation.

### Not obtained / needs login
- **BODC PAP collection** (partly done 2026-10-08: bottle data 2013-2019 now in section 6) (all post-2019 mooring data, recovered sensor data incl. SUNA/ISUS nitrate, pH (SAMI/SeaFET),
  pCO2 2016+, PAR, CTD casts and nutrients from PAP cruises, sediment trap fluxes after 2005): BODC asks you to register
  or log in before downloading. NOC says to search https://www.bodc.ac.uk/data/bodc_database/nodb/search/ for
  "Porcupine" or collection #5192 (the direct link https://www.bodc.ac.uk/data/bodc_database/nodb/data_collection/5912
  returns 404). EDMED record: https://www.bodc.ac.uk/resources/inventories/edmed/report/5912/
- NOC near-real-time data FTP (ftp://ftp.noc.soton.ac.uk/pub/animate/pap/…) – path no longer exists (550);
  https://apps.noc.ac.uk/pap/ returns 403 (the deployment pages at https://projects.noc.ac.uk/pap/outputs/time-series-data
  only show PNG plots, 2002–2025).
- Copernicus Marine in-situ TAC (NRT PAP buoy) – needs registration; not tried.
- **Re-tried 2026-10-08 (all failed or gated):**
  - OceanSITES GDAC (ftp://ftp.ifremer.fr/ifremer/oceansites/DATA/PAP/): unchanged, last deployment file is
    `OS_PAP-1_201907_R_{CTDO,PCO2}.nc`; no PAP dir in DATA_GRIDDED. `data-oceansites.ifremer.fr` doesn't resolve;
    NCEI has no OceanSITES PAP mirror (404). Ifremer ERDDAP (erddap.ifremer.fr) has no PAP/OceanSITES mooring datasets.
  - EMODnet Physics ERDDAP (https://erddap.emodnet-physics.eu/erddap/): TCP connects but every request timed out (60-90 s).
  - BODC (www.bodc.ac.uk, incl. the open Published Data Library): local DNS returns SERVFAIL for www.bodc.ac.uk, so it
    was not reachable from this machine. The PAP collection also needs a BODC login.
  - NCEI OCADS PAP page: still only accessions 0312034 and 0312066 (both already used). PANGAEA search: nothing at PAP
    after 2005 (traps, 1980s-90s hydrography).
  - Copernicus Marine in-situ (NRT PAP buoy incl. nitrate/pH if any) needs a free Copernicus Marine account; not attempted.
- So there is still **no open PAP-SO nitrate/phosphate/silicate or pH/carbonate time series, and no mooring data after
  2019-07**. Manual steps: (1) log in at BODC (https://www.bodc.ac.uk/data/bodc_database/nodb/search/, search
  "Porcupine Abyssal Plain") and request the PAP-SO mooring series (SUNA/ISUS nitrate, SeaFET/SAMI pH, pCO2, CTD-O 2019+)
  and the PAP cruise bottle nutrients/DIC/TA; (2) or register at https://data.marine.copernicus.eu/register and use the
  In Situ TAC global in-situ product (INSITU_GLO_PHYBGCWAV_DISCRETE_MYNRT_013_030; look for the PAP platform); (3) or contact the PAP-SO PIs at NOC
  (S. Hartman, R. Lampitt) for the processed nutrient/pH time series.
- Near-PAP open alternatives found but **not ingested** (outside the requested scope or radius): EXPORTS North Atlantic,
  RRS Discovery DY131, May 2021, ~48.8-49.2N 14.5-15.0W (110-140 km E of the column): euphotic-zone NO3/SiO2/PO4 +
  14C/32Si production, BCO-DMO 893293 (public ERDDAP). Its dissolved TM data (BCO-DMO 954941) need a BCO-DMO login.
  PAP-SO April 2017 upper-ocean trap fluxes, Th-234 and in-situ pump POC (BCO-DMO 765835, 765859, 765850, public).
- GLODAP / SOCAT / Argo / satellite / GEOTRACES near PAP: handled by the synthesis fork.

### Caveats for model comparison
- O2 (OceanSITES) is umol/kg; JC087 nutrients and O2 are umol/L. Both need density conversion to compare with model mmol/m3.
- Chl is uncalibrated/factory-calibrated fluorescence-derived chl-a (subject to non-photochemical quenching in daytime).
- Native sampling: CTD ~30 min–2 h, pCO2 2–6 per day, Chl 2–24 per day. Gaps between deployments are common.
- Two Chl sensors (Wetlabs, Cyclops) overlap at 30 m in 2014–2016; both are kept.


<!-- from SOURCES_synthesis.md -->
## Synthesis-product subsets (all sites)

Accessed 2026-10-07. Reader: `scripts/std_synth.py` (called by `standardize_obs.py` for every site).
Download helpers: `scripts/fetch_synth.sh` (SOCAT, BGC-Argo; chunked + retried ERDDAP requests,
because this network silently truncates long transfers), `scripts/subset_glodap.py` (GLODAP).
Site centres: HOT 22.75N 158.0W; BATS 31.67N 64.17W; HydroS 32.17N 64.5W; Papa 50.1N 144.9W; PAP 49.0N 16.5W.
BATS and HydroS are ~60 km apart, so their synthesis subsets overlap (same samples can appear in both CSVs).
All synthesis rows carry their own `lat`/`lon` (sample position), so the offset from the model column
can be computed; filter on distance if a tighter match is wanted.

### GLODAPv2.2023 (source tag `glodap`)
- Files: per-ocean CSVs `GLODAPv2.2023_Atlantic_Ocean.csv.zip`, `GLODAPv2.2023_Pacific_Ocean.csv.zip`
  from https://glodap.info/glodap_files/v2.2023/ (also at NCEI OCADS accession 0283442).
  Zips kept only in the session scratchpad; the per-site subsets (all original columns, unmodified rows,
  |dlat|<=1 deg and |dlon|<=1 deg) are in `<SITE>/raw/glodap/GLODAPv2.2023_<Ocean>_within1deg.csv`.
- Radius: 1 deg box. GLODAP does NOT contain the HOT/BATS time-series bottle data themselves (only
  repeat-hydrography/cruise data passing nearby), so counts are small at BATS/HOT.
- Variables: temp (in-situ, no flag), salt, O2, NO3, NO2, PO4, SiO2, DIC, ALK, pH (two flavours:
  total scale @25C/0dbar and total scale in situ; see units column), DOC, DON, TDN, Chl (ug/kg) where present.
- Units: umol/kg (needs density to compare with model mmol/m^3). Depth = G2depth (m). Time UTC; missing
  hour/minute set to 00:00.
- Flags: kept only GLODAP flag f == 2 (acceptable). Flag 0 (interpolated/calculated) and 9 dropped.
  These are the *adjusted* (bias-corrected) GLODAP values.
- Citation: Lauvset, S. K., et al. (2024). The annual update GLODAPv2.2023: the global interior ocean
  biogeochemical data product. Earth Syst. Sci. Data, 16, 2047-2072, https://doi.org/10.5194/essd-16-2047-2024 ;
  data: Lauvset et al. (2023) GLODAPv2.2023, NCEI Accession 0283442, https://doi.org/10.25921/zyrq-ht66 .
  GLODAP asks users to also cite the original cruise data DOIs where used substantially (G2doi column in the raw subset).

### SOCAT v2026 (source tag `socat`)
- ERDDAP: https://data.pmel.noaa.gov/socat/erddap/tabledap/socat_v2026_fulldata (query template in
  `<SITE>/raw/socat/query_url.txt`); files `socat_v2026_box1deg_<y0>-<y1>.csv` (multi-year chunks for early years, yearly afterwards; years with no data have no file; multi-year requests often failed with HTTP 408 on this server).
- Query: lat/lon within +-1 deg of site, time >= 1990, WOCE_CO2_water == 2. Variables fetched:
  expocode, platform, dataset QC flag, time, lat, lon, depth, sal, temp, fCO2_recommended, WOCE flag.
- Tidy variables: fCO2 (uatm at SST), temp (SST, degC), salt. depth_m = SOCAT sample depth, or 5 m if blank.
- Flags: dataset QC flag A-D kept (E dropped), WOCE 2 only. qc column = "<datasetQC>/<WOCE>".
- Caveats: includes both ships and moorings. SOCAT contains the PMEL moorings (WHOTS/MOSEAN at ALOHA,
  BTM/Hog Reef near Bermuda, Station Papa, PAP), so it will partly DUPLICATE mooring pCO2 fetched by the
  site-specific readers; de-duplicate by source when combining. fCO2 ~ pCO2 - ~1 uatm. High-frequency
  (minutes-3 h) records: average to daily before comparing with the model.
- Citation: Bakker, D. C. E., et al. Surface Ocean CO2 Atlas Database Version 2026 (SOCATv2026),
  NOAA NCEI / www.socat.info. Acknowledgement requested: "The Surface Ocean CO2 Atlas (SOCAT) is an
  international effort, endorsed by IOCCP, SOLAS and IMBER, to deliver a uniformly quality-controlled
  surface ocean CO2 database. The many researchers and funding agencies responsible for the collection
  of data and quality control are thanked for their contributions to SOCAT." (www.socat.info)

### BGC-Argo synthetic profiles (source tag `bgcargo`)
- ERDDAP: https://erddap.ifremer.fr/erddap/tabledap/ArgoFloats-synthetic-BGC (Ifremer; query template in
  `<SITE>/raw/bgcargo/query_url.txt`); files `argo_sprof_box100km_<y>-<y+1>.csv` (one per year).
- Query: +-0.9 deg lat by +-0.9/cos(lat) deg lon box; reader then keeps profiles with great-circle
  distance <= 100 km of the site.
- Variables: *_ADJUSTED fields only: temp, salt (psal), O2 (doxy, umol/kg), NO3 (umol/kg), Chl
  (chla_adjusted, mg/m3, fluorescence-derived; known factor-of-~2 regional bias vs HPLC), pH (total scale,
  in situ), bbp (bbp700, m-1). depth_m = PRES_ADJUSTED in dbar (not converted to m; ~1% difference at depth).
  To keep the tidy files manageable, Argo temp/salt are written only on levels that also carry a good
  BGC value (raw files keep the full ~1-2 dbar CTD profiles).
- Flags: adjusted QC in {1,2,5,8} and position QC in {1,2,5,8}. Real-time-only parameters with no
  adjusted value (e.g. raw chla, unadjusted nitrate) are excluded, so recent profiles may be missing
  BGC variables until delayed-mode adjustment. Includes floats whose only BGC sensor is O2.
- Citation: Argo (2000). Argo float data and metadata from Global Data Assembly Centre (Argo GDAC).
  SEANOE. https://doi.org/10.17882/42182 . Acknowledgement requested: "These data were collected and made
  freely available by the International Argo Program and the national programs that contribute to it
  (https://argo.ucsd.edu, https://www.ocean-ops.org). The Argo Program is part of the Global Ocean
  Observing System."

### Satellite chlorophyll: ESA OC-CCI v6.0 monthly (source tag `satchl`)
- ERDDAP: https://oceanwatch.pifsc.noaa.gov/erddap/griddap/esa-cci-chla-monthly-v6-0 (NOAA PIFSC OceanWatch;
  CoastWatch-West pfeg ERDDAP was unreachable on access date). File `<SITE>/raw/satchl/occci_v6_monthly_3x3.csv`
  (actually a 4x4 block of 4-km pixels, +-0.0625 deg, Sep 1997 - latest month). Query in `query_url.txt`.
- Tidy: Chl_sat = median of valid pixels per month (mg m-3), depth 0, lat/lon = site; qc = "npix=<n>".
  Time stamp is the ERDDAP time of the monthly composite (early in the month), not the month centre.
- Caveats: OCx-blended surface (first optical depth) chl; cloud gaps (esp. Papa/PAP winter); compare to
  model surface-layer Chl.
- Citation: Sathyendranath, S., et al. (2019) An ocean-colour time series for use in climate studies:
  the experience of the Ocean-Colour Climate Change Initiative (OC-CCI). Sensors 19, 4285,
  doi:10.3390/s19194285; and "ESA Ocean Colour CCI dataset, Version 6.0, European Space Agency,
  available online at https://esa-oceancolour-cci.org/", accessed via NOAA PIFSC OceanWatch ERDDAP.

### GEOTRACES dissolved Fe: IDP2025 NOT downloaded (re-tried 2026-10-08); open BCO-DMO substitutes used
- IDP2025 (released Nov 2025; CC-BY 4.0 + Fair Data Use Statement; no login) is distributed as whole packages from the
  BODC Published Data Library (https://www.bodc.ac.uk/geotraces/data/idp2025/ ,
  doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c). On 2026-10-08 www.bodc.ac.uk didn't resolve from this machine
  (local DNS SERVFAIL), so the package was not fetched. webODV (https://geotraces.webodv.awi.de/) is reachable but is an
  interactive session app (no clean scripted export), so it was not used.
- Searched instead for open copies of the original datasets on BCO-DMO ERDDAP (2,700 datasets listed): **GA03 and
  GP15 dissolved Fe are not on BCO-DMO ERDDAP**. GP15 Fe-ligand data are "log in". There is no Canadian GEOTRACES
  P26 or GEOVIDE dFe there either (only GA01 dissolved Pb). CCHDO GA03/GP15/GEOVIDE files are hydrography only.
- Used instead (see per-site notes, reader `std_tm.py`): HOT `bcodmo_hot_metals` (dFe + 8 metals at ALOHA, 2020-2023)
  and BATS `bcodmo_bait_fe` / `bcodmo_bait_tm_bottle` (BAIT/GApr13 2019 dFe, Mn, Al). **Papa and PAP have no dFe.**
- **Manual step for IDP2025:** in a browser open https://www.bodc.ac.uk/geotraces/data/idp2025/ , choose the discrete
  sample data (seawater) in ODV or CSV/NetCDF form, accept the GEOTRACES Fair Data Use Statement and download (no
  account needed). Or log in/register at https://geotraces.webodv.awi.de/ and use the Data Extractor. Subset to
  within 1 deg of: BATS 31.67N 64.17W (GA03 2011 leg KN204-01 BATS station),
  ALOHA 22.75N 158W (GP15 ALOHA station, Oct 2018), Papa 50.1N 144.9W (Canadian GEOTRACES Line P occupations of P26),
  PAP 49N 16.5W (GA01 GEOVIDE 2014, check whether any station is within 1 deg). Variable Fe_D_CONC_BOTTLE (nmol/kg) and,
  if wanted, Mn/Zn/Cu/Ni/Cd/Co_D_CONC. Put the file in `<SITE>/raw/geotraces_idp2025/` and add a reader to `std_tm.py`.
- Citation: GEOTRACES Intermediate Data Product Group (2025). The GEOTRACES Intermediate Data Product 2025
  (IDP2025). NERC EDS British Oceanographic Data Centre NOC. doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c.

### Tidy-row counts from these products (after QC), 2026-10-07
- HOT: Argo O2 65k, NO3 30k, pH 19k, Chl 16k, bbp 19k (2002-2026); SOCAT fCO2 50k (1993-2025);
  GLODAP ~40-100 per var (1984-2002, few cruises); Chl_sat 344 months.
- BATS: Argo O2 31k, NO3 5.7k, Chl 22k, pH 12.5k (2007-2025, BGC mostly 2022+); SOCAT fCO2 210k;
  GLODAP 150-350 per var (2003-2021); Chl_sat 346.
- HydroS: as BATS (overlapping boxes): Argo O2 41k, NO3 8k; SOCAT fCO2 219k; GLODAP 300-530; Chl_sat 346.
- Papa: Argo O2 100k, NO3 27k, Chl 41k, pH 34k (2006-2025); SOCAT fCO2 73k; GLODAP 2.3-3.2k for
  nutrients/O2/T/S, DIC 1.4k, ALK 765 (1985-2019, Line P cruises); Chl_sat 327.
- PAP: Argo O2 27k, Chl 17k, NO3 1.4k (2011-2025), no adjusted pH; SOCAT fCO2 33k; GLODAP 200-1000
  (1981-2017); Chl_sat 338.


## GEOTRACES IDP2025 dissolved Fe (added 2026-10-08)
- Extracted with the GEOTRACES webODV data extractor (https://geotraces.webodv.awi.de/, IDP2025 > seawater;
  variables DEPTH, CTDPRS_UP_T_VALUE, Fe_D_CONC; 1926 stations), ASCII spreadsheet (1.1 MB zip, 34 MB txt),
  kept in `_geotraces_idp2025_raw/`. `scripts/std_geotraces.py split` writes per-site subsets of stations within
  120 km to `<SITE>/raw/geotraces_idp2025/`; the reader keeps SeaDataNet flags 1-2, units nmol/kg as reported.
- Papa: 144 Fe (GP02 2017 to 4238 m; GPpr07 Line P process study 2012-2020, upper 600 m).
  BATS and Hydrostation S: 56 each (GA02 2010, GA03 2011). GApr13 (BAIT 2019) is skipped here because it is
  already ingested from BCO-DMO by std_tm. HOT and PAP: no IDP2025 Fe within 120 km (nearest 370 km, 338 km).
- Citation: GEOTRACES Intermediate Data Product Group (2025). The GEOTRACES Intermediate Data Product 2025
  (IDP2025). NERC EDS British Oceanographic Data Centre NOC. doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c
  (CC-BY 4.0, GEOTRACES Fair Data Use Statement).

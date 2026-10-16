# Synthesis-product subsets (all sites)

Accessed 2026-10-07. Reader: `scripts/std_synth.py` (called by `standardize_obs.py` for every site).
Download helpers: `scripts/fetch_synth.sh` (SOCAT, BGC-Argo; chunked + retried ERDDAP requests,
because this network silently truncates long transfers), `scripts/subset_glodap.py` (GLODAP).
Site centres: HOT 22.75N 158.0W; BATS 31.67N 64.17W; HydroS 32.17N 64.5W; Papa 50.1N 144.9W; PAP 49.0N 16.5W.
Added 2026-10-09: SOTS 46.94S 142.28E; K2 47.0N 160.0E; CARIACO 10.50N 64.67W; DYFAMED 43.4167N 7.8667E
(see "Update 2026-10-09" below).
Added later on 2026-10-09: MUNIDA 45.83S 171.53E; OWSM 66.0N 2.0E (see "Update 2026-10-09 (b)" below).
BATS and HydroS are ~60 km apart, so their synthesis subsets overlap (same samples can appear in both CSVs).
All synthesis rows carry their own `lat`/`lon` (sample position), so the offset from the model column
can be computed; filter on distance if a tighter match is wanted.

## GLODAPv2.2023 (source tag `glodap`)
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

## SOCAT v2026 (source tag `socat`)
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

## BGC-Argo synthetic profiles (source tag `bgcargo`)
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

## Satellite chlorophyll: ESA OC-CCI v6.0 monthly (source tag `satchl`)
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

## GEOTRACES dissolved Fe: NOT downloaded
- IDP2025 (released Nov 2025; supersedes IDP2021) is open access (CC-BY 4.0, Fair Data Use Statement; no
  login), but only as whole-product packages from the BODC Published Data Library
  (https://www.bodc.ac.uk/geotraces/data/idp2025/ , doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c), with
  the download served through a JavaScript page, or via the interactive webODV subsetting service
  (https://geotraces.webodv.awi.de/). Neither is scriptable cleanly here; skipped.
- Relevant stations if fetched by hand: GA03 (BATS, 2010-11), GP15 (Station ALOHA, 2018), GP02/Canadian
  GEOTRACES at P26/Papa, GA01 GEOVIDE (2014) passes near PAP. Fe_D_CONC (nmol/kg).
- Citation: GEOTRACES Intermediate Data Product Group (2025). The GEOTRACES Intermediate Data Product 2025
  (IDP2025). NERC EDS British Oceanographic Data Centre NOC. doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c.

## Tidy-row counts from these products (after QC), 2026-10-07
- HOT: Argo O2 65k, NO3 30k, pH 19k, Chl 16k, bbp 19k (2002-2026); SOCAT fCO2 50k (1993-2025);
  GLODAP ~40-100 per var (1984-2002, few cruises); Chl_sat 344 months.
- BATS: Argo O2 31k, NO3 5.7k, Chl 22k, pH 12.5k (2007-2025, BGC mostly 2022+); SOCAT fCO2 210k;
  GLODAP 150-350 per var (2003-2021); Chl_sat 346.
- HydroS: as BATS (overlapping boxes): Argo O2 41k, NO3 8k; SOCAT fCO2 219k; GLODAP 300-530; Chl_sat 346.
- Papa: Argo O2 100k, NO3 27k, Chl 41k, pH 34k (2006-2025); SOCAT fCO2 73k; GLODAP 2.3-3.2k for
  nutrients/O2/T/S, DIC 1.4k, ALK 765 (1985-2019, Line P cruises); Chl_sat 327.
- PAP: Argo O2 27k, Chl 17k, NO3 1.4k (2011-2025), no adjusted pH; SOCAT fCO2 33k; GLODAP 200-1000
  (1981-2017); Chl_sat 338.

## Update 2026-10-08

## GEOTRACES dissolved Fe: IDP2025 NOT downloaded (re-tried 2026-10-08); open BCO-DMO substitutes used
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

## Update 2026-10-09: four new sites (SOTS, K2, CARIACO, DYFAMED)

- Same products, queries, boxes and QC as above; sites added to the tables in `fetch_synth.sh`, `subset_glodap.py`,
  `std_synth.py` and `std_geotraces.py`. SOCAT and BGC-Argo via `fetch_synth.sh <root> socat|bgcargo SOTS K2 CARIACO
  DYFAMED`; OC-CCI with the same +-0.0625 deg query form (`<SITE>/raw/satchl/query_url.txt`). Accessed 2026-10-09.
- `subset_glodap.py` and `std_geotraces.py split` now take an optional `OBS_SITES="SOTS K2 ..."` environment variable to
  subset only those sites (default: all sites), so the existing sites' raw subsets were not rewritten.
- GLODAPv2.2023: SOTS rows are in the **Indian** Ocean file (548 rows; nothing in the Pacific file), K2 in the Pacific
  file (1903 rows), CARIACO in the Atlantic file (72 rows). The Mediterranean is in the Atlantic file (plus one 1977
  cruise in the Indian file), but **no GLODAP station lies within 1 deg of DYFAMED** (nearest ~440 km, 2011).
- GEOTRACES IDP2025 Fe within 120 km: SOTS 159 raw rows (GIpr10 2020, GPpr11 2016, GSc01 2017), K2 446 (GP02, GP18
  2011). **None at CARIACO (nearest 920 km) or DYFAMED (nearest 210 km, GApr09 2017).**
- No BGC-Argo profiles within 100 km of CARIACO (enclosed basin); SOCAT at CARIACO is a single 2011 year.

## Tidy-row counts from these products (after QC), 2026-10-09
- SOTS: Argo O2 126k, Chl 63k, bbp 23k, pH 15k, NO3 6.1k (2005-2026, BGC 2010+); SOCAT fCO2 148k (1993-2025);
  GLODAP 45-513 per var (1998-2016; DIC 88, ALK 87); IDP2025 Fe 54 (2016-2020); Chl_sat 321 months.
- K2: Argo O2 101k, NO3 7.1k (2017-2026; no adjusted Chl, pH or bbp); SOCAT fCO2 5.1k (1995-2025); GLODAP 1.1-2.1k
  per var for nutrients/O2/T/S/DIC/ALK/pH (1985-2021), Chl 99, DOC 66; IDP2025 Fe 85 (2011-2017); Chl_sat 326.
- CARIACO: no Argo; SOCAT fCO2 928 (2011 only); GLODAP 5-71 per var (1997-2003); no IDP2025 Fe; Chl_sat 346.
- DYFAMED: Argo Chl 274k, O2 261k, bbp 243k, NO3 13k (2012-2026), no adjusted pH; SOCAT fCO2 50k (1998-2023);
  no GLODAP; no IDP2025 Fe; Chl_sat 346.


## Update 2026-10-09 (b): two more sites (MUNIDA, OWSM)

- Same products, queries, boxes and QC as above; MUNIDA (45.83S 171.53E) and OWSM (66.0N 2.0E) added to the site tables
  in `fetch_synth.sh`, `subset_glodap.py`, `std_synth.py` and `std_geotraces.py`, and to `standardize_obs.py`
  (readers `std_munida` / `std_owsm` + `std_synth` + `std_geotraces`). SOCAT and BGC-Argo via
  `fetch_synth.sh <root> socat|bgcargo MUNIDA OWSM`; OC-CCI with the same +-0.0625 deg query form
  (`<SITE>/raw/satchl/query_url.txt`); GLODAP and GEOTRACES subsets with `OBS_SITES="MUNIDA OWSM"`, so no other site's
  raw subsets or tidy CSVs were rewritten (md5 of the other nine `<SITE>_obs.csv` checked unchanged). Accessed 2026-10-09.
- GLODAPv2.2023: OWSM rows are in the **Arctic** Ocean file (GLODAPv2.2023_Arctic_Ocean.csv.zip, which holds the Nordic
  Seas; nothing in the Atlantic file), 2419 rows from 9 cruises (1994-2021). 2098 of them are one compiled "cruise"
  58P320011031 (weather ship Polarfront, Station M bottle series 2001-2021), very likely the same Skjelvan DIC/ALK/nutrient
  bottles as the site source `ncei_ocads_owsm` (NCEI 0112884; not checked row by row): **deduplicate glodap vs ncei_ocads_owsm at OWSM**.
  **No GLODAP station lies within 1 deg of MUNIDA** (nearest ~540 km, a 2009 cruise at 41.9S 175.5E).
- GEOTRACES IDP2025 Fe within 120 km: **none at MUNIDA (nearest ~900 km, GIPY01 2006 / GPpr02 2010) or OWSM (nearest
  ~950 km, GN04 2015)**; the per-site TSVs are header-only.
- SOCAT at OWSM: M/V Polarfront (weather ship on station, 2006-2009, 347k rows), the Station M mooring (2011-2021, 12.9k)
  and G.O. Sars (2006-2025, 6.5k); 348k of the 366k fCO2 rows are within 30 km of the station. **This is the only surface
  pCO2/fCO2 at OWSM**: the site reader's NorEMSO HydroC fCO2 (`nmdc_noremso_mooring`) is at ~470-580 m.
- SOCAT at MUNIDA: the 1 deg box covers the whole Otago shelf/Subtropical Front, so only 4352 of the 45.9k fCO2 rows are
  within 15 km of the station (the `std_munida` cut). SOCAT 2004-2006 holds 4 Munida (61MU) and 12 Polaris II (61P2)
  cruises whose dates match the NCEI underway files read by `std_munida` (`ncei_ocads_0100218` M040127-M060730,
  `ncei_ocads_0112883` P060830, P061024), so those **2004-2006 SOCAT rows (10.5k in the box, 2300 within 15 km) very
  likely duplicate the site's NCEI underway pCO2** (not checked row by row). Polaris II (to 2025) and Tangaroa (2013-2024) SOCAT
  cruises continue, extending the underway record past the open NCEI files.
- BGC-Argo at MUNIDA is sparse (O2 1112 levels 2016-2026; Chl/NO3/pH/bbp 313 each, 2017-2026). At OWSM, adjusted NO3 and pH
  come from a single day (2022-09-11).

## Tidy-row counts from these products (after QC), 2026-10-09 (b)
- MUNIDA: Argo O2 1.1k, Chl/NO3/pH/bbp 313 each (2016-2026); SOCAT fCO2 45.9k (1994-2025); no GLODAP; no IDP2025 Fe;
  Chl_sat 339 months.
- OWSM: Argo T/S 70k, Chl 61k (2013-2022), O2 55k (2006-2025), bbp 54k, NO3 201, pH 447 (2022 only); SOCAT fCO2 366k
  (1995-2025); GLODAP 2.4k T/S, DIC 2.1k, ALK 1.5k, NO3/SiO2/PO4/O2 1.1-1.3k (1994-2021; mostly the Station M bottles),
  pH 42, Chl 12; no IDP2025 Fe; Chl_sat 226 months (winter gaps at 66N).

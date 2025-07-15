import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

def open_CGEM(filepath):
    # open C-GEM outputs
    try:
        # Assuming data is in a structured format like CSV or similar
        df = pd.read_csv(filepath, delimiter='\t', header=None)
        # Process the DataFrame
        # put timestep (s) as index
        df.index = df.values[:, 0]
        # To delete first and last columns
        df.drop(columns=df.columns[0], axis=1, inplace=True)
        df.drop(columns=df.columns[-1], axis=1, inplace=True)
        if(filepath[len(filepath)-9:len(filepath)]=='depth.dat'):
            df.drop(columns=df.columns[0], axis=1, inplace=True)
        elif(filepath[len(filepath)-9:len(filepath)]=='width.dat'):
            df.drop(columns=df.columns[0], axis=1, inplace=True)
        elif (filepath[len(filepath) - 5:len(filepath)] == 'U.dat'):
            df.drop(columns=df.columns[0], axis=1, inplace=True)
        df.columns = np.arange(0, np.size(df,1), 1)
        print(df)
    except FileNotFoundError:
        print(f"File '{filepath}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return df
def compute_flux(U,depth,width,concentration,start,end):
    # compute upstream and downstream fluxes in Tg yr-1
    # negative: towards the ocean
    # positive: towards the river source
    # U : horizontal velocity in m/s
    # depth: depth in m
    # width = channel width in m
    # concentration: solutes concentration in mmol m-3
    # start: idx of beginning of period of interest
    # end: idx of end of period of interest
    sec_in_yr = 60 * 60 * 24 * 365
    mol2g_C = 12.01070
    mg2Tg = 1E-15
    up_idx = np.size(concentration, 1) - 1
    dw_idx = 0
    saving_ts = concentration.index[1]-concentration.index[0]
    up_flx = (np.sum(U[start:end][[up_idx]] *
                     concentration.loc[start:end][[up_idx]] *
                     depth[start:end][[up_idx]] *
                     width[start:end][[up_idx]]) *
              sec_in_yr * mol2g_C * mg2Tg) / saving_ts
    dw_flx = (np.sum(U[start:end][[dw_idx]] *
                     concentration.loc[start:end][[dw_idx]] *
                     depth[start:end][[dw_idx]] *
                     width[start:end][[dw_idx]]) *
              sec_in_yr * mol2g_C * mg2Tg) / saving_ts
    return float(up_flx), float(dw_flx)
def integrate_CO2(depth,width,FCO2,start,end,DELXI):
    # integrate CO2 fluxes over the entire estuary [Tg C yr-1]
    # negative output means CO2 outgassing to the atmosphere
    # depth: depth in m
    # width = channel width in m
    # FCO2: CO2 flux in [mmol C m^−3 s^−1]
    # start: idx of beginning of period of interest
    # end: idx of end of period of interest
    # DELXI: length of estuarine section in m
    sec_in_yr = 60 * 60 * 24 * 365
    mol2g_C = 12.01070
    mg2Tg = 1E-15
    saving_ts = FCO2.index[1]-FCO2.index[0]
    int_CO2 = np.sum(np.sum(FCO2.loc[start:end]*
                     depth[start:end] *
                     width[start:end] *
                     DELXI) * sec_in_yr * mol2g_C * mg2Tg / saving_ts)
    return int_CO2
def make_arrows_flx(up_flx,dw_flx,max_scale,order,variable_name,title):
    # up_flx, dw_flx: upstream and downstream fluxes
    # max_scale: scale factor for arrow width
    # order: set vertical position 1=top, 2 =middle, 3=bottom
    # variable_name: variable name as string
    # estuary
    rectangle = mpatches.Rectangle((0.2, 0), 0.6, 0.6, edgecolor = 'k', facecolor = 'w')
    axs.add_patch(rectangle)
    axs.annotate('Estuary', (0.5, 0.3), color='k', weight='bold',
                fontsize=12, ha='center', va='center')
    axs.annotate('Upstream', (0.1, 0.8), color="C1", weight='bold',
                 fontsize=12, ha='center', va='center')
    axs.annotate('Downstream', (0.95, 0.8), color="C2", weight='bold',
                 fontsize=12, ha='center', va='center')
    up = up_flx/max_scale
    dw = dw_flx/max_scale
    # upstream
    if (up < 0):
        x_tail = 0.
        x_head = 0.2
        x = x_head
    elif (up > 0):
        x_tail = 0.2
        x_head = 0.
        x = x_head - 0.1
    if (order == 1):
        y_tail = 0.5
        y_head = 0.5
    elif (order == 2):
        y_tail = 0.5 - .2
        y_head = 0.5 - .2
    elif (order == 3):
        y_tail = 0.5 - .4
        y_head = 0.5 - .4
    dx = x_head - x_tail
    dy = y_head - y_tail
    width = np.max([abs(up), 0.002])
    arrow = mpatches.FancyArrow(x_tail, y_tail, dx, dy,
                                width=width, length_includes_head=True, color="C1"
                                )
    axs.add_patch(arrow)
    axs.annotate(variable_name, (0.05, y_tail+dy+width*3), color='k', weight='bold',
                 fontsize=12, va='center')
    axs.annotate(str(round(abs(up_flx),2)), (x, y_head), color='k', weight='bold',
                 fontsize=12, va='center',annotation_clip=False)
    # downstream
    if (dw < 0):
        x_tail = 0.8
        x_head = 1
        x = x_head
    elif (dw > 0):
        x_tail = 1.0
        x_head = 0.8
        x = x_head - 0.1
    if (order == 1):
        y_tail = 0.5
        y_head = 0.5
    elif (order == 2):
        y_tail = 0.5 - .2
        y_head = 0.5 - .2
    elif (order == 3):
        y_tail = 0.5 - .4
        y_head = 0.5 - .4
    dx = x_head - x_tail
    dy = y_head - y_tail
    width = np.max([abs(dw), 0.002])
    arrow = mpatches.FancyArrow(x_tail, y_tail, dx, dy,
                                width=width, length_includes_head=True, color="C2"
                                )
    axs.add_patch(arrow)
    axs.annotate(variable_name, (0.85, y_tail+dy+width*3), color='k', weight='bold',
                 fontsize=12, va='center')
    axs.annotate(str(round(abs(dw_flx),2)), (x, y_head), color='k', weight='bold',
                 fontsize=12, va='center',annotation_clip=False)
    if (title != None):
        plt.title(title, weight='bold',
                  fontsize=15, ha='center', va='center')
    # hide x-axis
    axs.get_xaxis().set_visible(False)
    # hide y-axis
    axs.get_yaxis().set_visible(False)
    # hide axes and borders
    plt.axis('off')
def make_arrows_airwaterflx(airwater_flx,max_scale,variable_name):
    # flx: air-water flux
    # max_scale: scale factor for arrow width
    # variable_name: variable name as string
    flx = airwater_flx/max_scale
    # air-water
    x_head = 0.5
    x_tail = 0.5
    if (flx < 0):
        y_tail = 0.6
        y_head = 0.8
        y = y_head
    elif (flx > 0):
        x_tail = 0.8
        x_head = 0.6
        y = 0.4
    dx = x_head - x_tail
    dy = y_head - y_tail
    width = np.max([abs(flx), 0.002])
    arrow = mpatches.FancyArrow(x_tail, y_tail, dx, dy,
                                width=width, length_includes_head=True, color="C3"
                                )
    axs.add_patch(arrow)
    axs.annotate(variable_name, (0.35, 0.7), color='k', weight='bold',
                 fontsize=12)
    axs.annotate(str(round(abs(flx),3)), (0.45, y), color='k', weight='bold',
                 fontsize=12,annotation_clip=False)

#load depth
depth = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/depth.dat")
#load width
width = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/width.dat")
#load horizontal velocity
U = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/U.dat")

#load CO2 flx [mmol C m^−3 s^−1]
CO2flx = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/FCO2.dat")
idx_lastyear = CO2flx.index[0]
idx_final = CO2flx.index[len(CO2flx)-1]
DELXI = 1000
FCO2_int = integrate_CO2(depth,width,CO2flx,idx_lastyear,idx_final,DELXI)

#load DIC concentration [mmol m^-3]
DIC = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/DIC.dat")
#compute flux [Tg yr-1]
up_DIC, dw_DIC = compute_flux(U,depth,width,DIC,idx_lastyear,idx_final)
#load TOC concentration [mmol m^-3]
TOC = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/TOC.dat")
#compute flux [Tg yr-1]
up_TOC, dw_TOC = compute_flux(U,depth,width,TOC,idx_lastyear,idx_final)
#load ALK concentration [mmol m^-3]
ALK = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/ALK.dat")
#compute flux [Tg yr-1]
up_ALK, dw_ALK = compute_flux(U,depth,width,ALK,idx_lastyear,idx_final)

fig, axs = plt.subplots(nrows=1)
make_arrows_flx(up_DIC,dw_DIC,300,1,'DIC','Carbon budget (flux in Tg C yr$^{-1}$)')
make_arrows_flx(up_ALK,dw_ALK,300,2,'ALK',None)
make_arrows_flx(up_TOC,dw_TOC,300,3,'TOC',None)
make_arrows_airwaterflx(FCO2_int,300,'Air-water CO$_2$ flux')
plt.savefig('Carbon_budget.png', dpi=180)


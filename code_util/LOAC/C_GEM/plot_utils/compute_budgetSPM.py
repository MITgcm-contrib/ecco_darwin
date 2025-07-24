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
    # concentration: solutes concentration in [g/l]
    # = solutes concentration in [kg/m3]
    # start: idx of beginning of period of interest
    # end: idx of end of period of interest
    sec_in_yr = 60 * 60 * 24 * 365
    kg2Tg = 1E-9
    up_idx = np.size(concentration, 1) - 1
    dw_idx = 0
    saving_ts = concentration.index[1]-concentration.index[0]
    up_flx = (np.sum(U[start:end][[up_idx]] *
                     concentration.loc[start:end][[up_idx]] *
                     depth[start:end][[up_idx]] *
                     width[start:end][[up_idx]]) *
              sec_in_yr * kg2Tg) / len(concentration.loc[start:end])
    dw_flx = (np.sum(U[start:end][[dw_idx]] *
                     concentration.loc[start:end][[dw_idx]] *
                     depth[start:end][[dw_idx]] *
                     width[start:end][[dw_idx]]) *
              sec_in_yr * kg2Tg) / len(concentration.loc[start:end])
    return float(up_flx), float(dw_flx)
def integrate_rates(width,rates,start,end,DELXI):
    # integrate fluxes/rates over the entire estuary [Tg C yr-1]
    # depth: depth in m
    # width = channel width in m
    # rates: flux/rates in [kg m^-2 s^-1]
    # start: idx of beginning of period of interest
    # end: idx of end of period of interest
    # DELXI: length of estuarine section in m
    sec_in_yr = 60 * 60 * 24 * 365
    kg2Tg = 1E-9
    saving_ts = rates.index[1]-rates.index[0]
    int_rates = np.sum(np.sum(rates.loc[start:end]*
                     width[start:end] *
                     DELXI) * sec_in_yr * kg2Tg / len(rates.loc[start:end]))
    return int_rates
def make_arrows_flx(up_flx,dw_flx,max_scale,variable_name,title):
    # up_flx, dw_flx: upstream and downstream fluxes
    # max_scale: scale factor for arrow width
    # order: set vertical position 1=right, 2 =left
    # variable_name: variable name as string
    # estuary
    rectangle = mpatches.Rectangle((0.2, 0), 0.6, 0.6, edgecolor = 'k', facecolor = 'w')
    axs.add_patch(rectangle)
    axs.annotate('Estuary', (0.5, 0.3), color='k', weight='bold',
                fontsize=12, ha='center', va='center')
    axs.annotate('Upstream', (0.1, 0.8), color="C1", weight='bold',
                 fontsize=12, ha='center', va='center')
    axs.annotate('Downstream', (0.95, 0.8), color="b", weight='bold',
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
    y_tail = 0.5
    y_head = 0.5
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
    y_tail = 0.5
    y_head = 0.5
    dx = x_head - x_tail
    dy = y_head - y_tail
    width = np.max([abs(dw), 0.002])
    arrow = mpatches.FancyArrow(x_tail, y_tail, dx, dy,
                                width=width, length_includes_head=True, color="b"
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
def make_rates(rates,max_scale,variable_name,order):
    # rates: deposition/erosion of SPM
    # max_scale: scale factor for arrow width
    # variable_name: variable name as string
    flx = rates/max_scale
    # water-sediment
    if(order==1):
        x_head = 0.3
        x_tail = 0.3
        x = 0.25
    elif(order==2):
        x_head = 0.7
        x_tail = 0.7
        x = 0.65
    if (rates < 0):
        y_tail = 0.
        y_head = 0.2
        y = y_head
    elif (rates > 0):
        y_tail = 0.2
        y_head = 0.
        y = 0.55
    dx = x_head - x_tail
    dy = y_head - y_tail
    width = np.max([abs(flx), 0.002])
    arrow = mpatches.FancyArrow(x_tail, y_tail, dx, dy,
                                width=width, length_includes_head=True, color="gray"
                                )
    axs.add_patch(arrow)
    axs.annotate(variable_name, (x-0.05, y_tail), color='k', weight='bold',
                 fontsize=12, va='center')
    axs.annotate(str(round(abs(rates),2)), (x, y_head), color='k', weight='bold',
                 fontsize=12, va='center',annotation_clip=False)

#load depth
depth = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/depth.dat")
#load width
width = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/width.dat")
#load horizontal velocity
U = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/U.dat")

#load erosion [mg m^-2 s^-1]
ero = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/erosion.dat")
idx_lastyear = ero.index[0]
idx_final = ero.index[len(ero) - 1]
DELXI = 1000
ero_int = -integrate_rates(width,ero,idx_lastyear,idx_final,DELXI)
#load deposition [mg m^-2 s^-1]
depo = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/deposition.dat")
depo_int = integrate_rates(width,depo,idx_lastyear,idx_final,DELXI)

#load SPM concentration [g/l]
SPM = open_CGEM("/Users/rsavelli/Documents/CMS_LOAC/Guayas/outputs/SPM.dat")
#compute flux [Gg yr-1]
up_SPM, dw_SPM = compute_flux(U,depth,width,SPM,idx_lastyear,idx_final)

fig, axs = plt.subplots(nrows=1)
make_arrows_flx(up_SPM,dw_SPM,1E3,'SPM','Sediment budget (flux in Tg yr$^{-1}$)')
make_rates(depo_int,1E3,'Deposition',1)
make_rates(ero_int,1E3,'Erosion',2)
plt.savefig('Sediment_budget.png', dpi=180)


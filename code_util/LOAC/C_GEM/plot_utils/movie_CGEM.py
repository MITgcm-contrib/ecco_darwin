from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.preprocessing import normalize
import subprocess
import shlex

# also requires ffmpeg installed on your machine

# open C-GEM outputs
def open_CGEM(filepath):
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
# create colorline plot
def plot_colourline(x,y,c,line_width,minval,maxval):
    col = cm.viridis((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=col[i], linewidth = line_width[i],solid_capstyle='round')
    im = ax.scatter(x, y, c=c, s=0, cmap=cm.viridis, vmin = minval, vmax= maxval)
    return im
# create animation
def make_CGEManim(filepath, # path to C-GEM outputs
                  shapepath, # path to river estuary shapefile.
                             # It must be equidistant N points with N matching model grid.
                             # Point must be labelled with id field increasing from upstream to downstream.
                             # Example, estuary 74000m-long with DELXI = 1000m is a 74-point grid.
                             # Shapefile must contain width of the river.
                             # (source: SWORD database https://www.swordexplorer.com)
                  lat0, # center map on lat0 (degree)
                  lon0, # center map on lon0 (degree)
                  width, # width of map window (meters)
                  height, # height of map window (meters)
                  rivername, # Name of the river for title
                  scale_riverwidth, # scale factor for width of river segments
                  variable_name, # Variable from C-GEM (ex Salinity...) for colorbar label
                  units, # Variable units (ex PSU) for colorbar label
                  minval, # Variable minimal bound for colorbar range (ex 0 for Salinity)
                          # Leave to None for adjustment to min value of output
                          # !! Output will be scaled to this range
                  maxval, # Variable maximal bound for colorbar range (ex 34 for Salinity)
                          # Leave to None for adjustment to max value of output
                          # !! Output will be scaled to this range
                  parallels_grid, # Define grid for parallels (degree)
                  meridians_grid, # Define grid for meridians (degree)
                  frame_count, # Number of frame to print
                  output_path,  # Output directory
                  movie_name # Name of final movie file
                  ):
    df = open_CGEM(filepath)
    shapefile = gpd.read_file(shapepath)
    gif_indx = 0
    if (frame_count == None):
        duration = len(df)
    else:
        if isinstance(frame_count, int):
            duration = frame_count
        else:
            raise Exception("ERROR: frame_count must be an integer")
    for pp in range(0,duration):
        fig, ax0 = plt.subplots(figsize=(6.5, 4.5))
        m1 = Basemap(width=width,height=height,projection='lcc',
                  lat_0=lat0, lon_0=lon0,resolution='f')
        m1.fillcontinents()
        # plot river lines
        order = np.argsort(shapefile.id)
        xs = np.array(shapefile.x)[order]
        ys = np.array(shapefile.y)[order]
        riverwidth = normalize(np.array(shapefile.width)[order].reshape(1, -1), norm="l2")
        if (minval == None):
            minval = np.min(df)
        if (maxval == None):
            maxval = np.max(df)
        im = plot_colourline(*m1(xs, ys), np.flip(df.values[pp,:]), riverwidth[0,:]*scale_riverwidth, minval, maxval)
        cbar = m1.colorbar(im)
        cbar.set_label(variable_name + ' ' + units,font='Arial',
                       size=14,weight="bold")
        if (round(df.index[pp]/(60*60*24),2) < 1) :
            plt.title(rivername + ' Estuary (t = {:0>4} h)'.format(str(round(df.index[pp] / (60 * 60), 2))))
        elif (round(df.index[pp]/(60*60*24),2) > 1) :
            plt.title(rivername + ' Estuary (t = {:0>4} h)'.format(str(round(df.index[pp]/(60*60*24),2))))
        elif (round(df.index[pp]/(60*60*24),2) > 365) :
            plt.title(rivername + ' Estuary (t = {:0>4} h)'.format(str(round(df.index[pp] / (60 * 60 * 24), 2)-365)))
            break
        # draw parallels and meridians.
        # label parallels on right and top
        # meridians on bottom and left
        parallels = parallels_grid
        # labels = [left,right,top,bottom]
        m1.drawparallels(parallels,labels=[True,False,True,False])
        meridians = meridians_grid
        m1.drawmeridians(meridians,labels=[True,False,False,True])
        #plt.pause(0.1)
        # gif_maker('blue_marble_rotating_globe.gif','.',gif_indx,len(dates)-1,dpi=90)
        gif_indx += 1
        plt.savefig(output_path + '/frame_' + '{:0>4}'.format(str(gif_indx)) + '.png', dpi=180)
        plt.close()

    # compute animation movie from individual frames with ffmpeg
    subprocess.run(shlex.split('ffmpeg -r 5 -s 1135x1072 -pattern_type glob -i '+ output_path+'"/frame_*.png" -vcodec '
                               'libx264  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -crf 25  -pix_fmt yuv420p '
                               + output_path+'/'+movie_name))

make_CGEManim(filepath = '/Users/rsavelli/Documents/C_GEM/S.dat',
              shapepath = "/Users/rsavelli/Documents/CMS_LOAC/Guayas/Guayas2.shp",
              lat0 = -2.4,
              lon0 = -80,
              width = 120000,
              height = 95000,
              rivername = 'Guayas',
              scale_riverwidth = 50,
              variable_name = 'Salinity',
              units = '(PSU)',
              minval = None,
              maxval = None,
              parallels_grid = np.arange(-3,-2,.3),
              meridians_grid = np.arange(-82.,-79.,.3),
              frame_count = 500,
              output_path = './Guayas/movie_frame',
              movie_name = 'Salinity.mp4')

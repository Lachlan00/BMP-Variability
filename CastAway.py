"""
Custom class methods for CastAway CTD surveys
---------------------------------------------
Author: Lachlan Phillips
Email: lachlan.phillips.00@gmail.com

This class object is specifically desinged for wokring with data speciifc to our 
Linkage Grant study and will not easily translate to othr CastAway datasets without
some modification. Feel free to contact me should you like to adpat this code for other
projects and I can help you out.

Note: The duplicate meta data can be confusing as ctd objects contain meta information
      but I also pass in more complete meta data data frames. This is mainly because
      python-ctd  handles metadata as pandas data frame attributes where I have my own
      system of meta data handling I use. Consequently, a Casts object will have a metadata
      attached to each cast and a complete set attached to the Casts object itself. 
"""
from itertools import compress
from datetime import datetime
import scipy.interpolate
import numpy as np
import pandas as pd
import numpy as np
import pytz
import math
import cmocean
import re
import pyreadr
from scipy.spatial import KDTree
from progressbar import ProgressBar

import ctd
from oceans.datasets import etopo_subset

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import pdb

class Casts:
    """
    data:
        A list of python-ctd casts (pandas data frames).
    meta:
        A pandas data frame containing meta data.
    tz:
        Local timezone of casts.
    """
    def __init__(self, data, meta, tz='Australia/Sydney'):
        self.data = data
        self.meta = meta

        # format data
        self.meta.loc[:,'cast_time_UTC'] = [pytz.utc.localize(datetime.strptime(t, '%Y-%m-%d %H:%M:%S')) for t in self.meta['cast_time_UTC']]
        self.meta.loc[:,'cast_time_local'] = [pytz.timezone(tz).localize(datetime.strptime(t, '%Y-%m-%d %H:%M:%S')) for t in self.meta['cast_time_local']]

        # format metadata (and update from meta information)
        for i in range(0, len(self.data)):
            self.data[i]._metadata['Start longitude'] = float(self.meta.iloc[i]['start_longitude'])
            self.data[i]._metadata['Start latitude'] = float(self.meta.iloc[i]['start_latitude'])


        # Remove meta with no casts and vice versa
        # ____add later____

    # Get max depth for all casts
    def max_depth(self):
        return max([max(df.index.values.astype(float)) for df in self.data])

    # Get max and min variable
    def variable_range(self, variable):
        minVal = min([min(df[variable]) for df in self.data])
        maxVal = max([max(df[variable]) for df in self.data])
        return (minVal, maxVal)


    # def extract station and transects
    def expand_transect_stations(self, meta):
        """
        Extract out transect codes into trasects and stations
        e.g. T3S1  tp T3, S1
        """
        meta['transect'] = [t[0:2] for t in meta['transect_id']]
        meta['station'] = [t[2:4] for t in meta['transect_id']]
        return meta

    def return_section_names(self, order=True):
        meta = self.meta
        meta = self.expand_transect_stations(self.meta)
        # join surveys and transects
        sections = set([x[0]+'-'+x[1] for x in zip(meta.survey_id, meta.transect)])
        surveys = [s.split('-')[0] for s in sections]
        transects = [s.split('-')[1] for s in sections]

        df = pd.DataFrame.from_dict({
            'survey': surveys,
            'transect': transects,
            'year': [int(s.split('-')[0].split('_')[0]) for s in sections]
            })

        if order:
            df['repeat'] = [int(re.sub("[^0-9]", "", s.split('-')[0].split('_')[1])) for s in sections]
            df['tval'] = [int(s.split('-')[1][1]) for s in sections]
            df = df.sort_values(['year', 'repeat', 'tval'], ascending=[True, True, False])
            df = df.drop(['repeat', 'tval'], axis=1)

        return(df.reset_index(drop=True))
        

    def transect_runs(self, survey, transect):
        """
        For a given Survey and transect return all runs on that transect
        """
        meta = self.meta
        meta = self.expand_transect_stations(meta)
        # Filter to survey and transect
        meta = meta[meta['survey_id'] == survey]
        meta = meta[meta['transect'] == transect]

        # work out how many repeats
        runs = list(meta.groupby([t.date() for t in meta['cast_time_local']]))

        # HACK 2016_S2 T4 with multiple runs on same day
        # Later separate into better transects
        if survey == '2016_S2' and transect == 'T4':
            runs[2] = (runs[2][0], runs[2][1][0:4])

        # remove runs with only one cast
        runs = [run for run in runs if len(run[1]) > 1]

        return runs

    def swarm_environment_data(self, interp_limits=(1.875, 5), set_crossshore_dist=15):
        """
        Returns xsection interporaltion data
        """
        sections = self.return_section_names()
        data = self.data
        meta = self.meta

        # dict to hold outputs
        output = [None]*10000

        # calc max depth
        maxDepth = self.max_depth()
        outi = 0
        print('Interporlating xsections..')
        for idx, row in sections.iterrows():
            # get runs
            runs = self.transect_runs(row.survey, row.transect)
            print(row.survey+' '+row.transect)
            for run in runs:
                # get casts
                transect_casts = [list(compress(data, meta['file_name'] == fn))[0] for fn in run[1]['file_name']]
                meta_section = [df._metadata['File name'] for df in transect_casts]
                meta_section = meta[meta['file_name'].isin(meta_section)]
                # get longitudes
                lons = meta_section['start_longitude']
                # just get one lat as all lats the same more or less (will need when we calculate distance)
                lat = round(meta_section['start_latitude'].mean(), 4)
                # calculate distance
                x = [round(harversine(min(lons), lat, lon, lat), 2) for lon in lons]
                # add distance to each cast
                for i in range(0, len(transect_casts)):
                    transect_casts[i]['x'] = x[i]  
                # Create a master data set
                df = pd.concat(transect_casts)    

                # get max depths
                cast_pos = x
                max_depths = [math.ceil(max(cast.index.values.astype(float))) for cast in transect_casts]  

                # extract data
                x = np.asanyarray(df.x)
                y = df.index.values.astype(float)
                val_temp = np.asanyarray(df['temperature'])
                val_salt = np.asanyarray(df['salinity'])

                # set plot limits (so all cross sections equal distance)
                maxDist = set_crossshore_dist
                if set_crossshore_dist is None:
                    maxDist = x.max()

                # create interporlations
                xi, yi = np.mgrid[x.min():maxDist:80j, y.min():maxDepth:80j]

                # creat option to choose interporaltion method later
                #val_interp = normal_interp(x, y, val, xi, yi)
                val_temp_interp = rescaled_interp(x, y, val_temp, xi, yi)
                val_temp_interp = np.ma.masked_invalid(val_temp_interp)
                val_salt_interp = rescaled_interp(x, y, val_salt, xi, yi)
                val_salt_interp = np.ma.masked_invalid(val_salt_interp)

                # make a mask for values far away from real values
                cast_times = meta_section.cast_time_UTC.reset_index(drop=True)
                dtUTC = [None]*xi.size
                k = 0
                for i in range(0, xi.shape[0]):
                    for j in range(0, xi.shape[1]):
                        nearest_cast_depth = max_depths[np.argmin(abs(cast_pos - xi[i, j]))]
                        dtUTC[k] = cast_times[np.argmin(abs(cast_pos - xi[i, j]))]
                        if (xi[i, j] < (x.max() + interp_limits[0])) & (yi[i, j] < (nearest_cast_depth + interp_limits[1])):
                            continue
                        else:
                            val_temp_interp[i, j] = np.nan
                            val_salt_interp[i, j] = np.nan
                        k += 1

                # make data frame
                output[outi] = pd.DataFrame.from_dict({
                    'transect_distance':xi.ravel(),
                    'depth':yi.ravel(),
                    'temp':val_temp_interp.ravel(),
                    'salt':val_salt_interp.ravel(),
                    'lat':lat,
                    'lon':[reverse_harversine(min(lons), lat, d) for d in xi.ravel()],
                    'dtUTC':dtUTC,
                    'survey':row.survey,
                    'transect':row.transect
                })
                outi += 1

        # Make into single dataframe
        output = [df for df in output if df is not None]
        df_interp = pd.concat(output)

        # load and get values for acoustics
        # load acoustic data
        result = pyreadr.read_r('./data/surveys/acoustics/agg.rds')
        agg = result[None]
        agg['dtUTC'] = [pytz.utc.localize(datetime.strptime(str(d)+'_'+str(t).replace(' ',''), "%Y%m%d_%H:%M:%S.%f")) for d, t in zip(agg['Date_M'], agg['Time_M'])]
        agg['dt_local'] = [dt.replace(tzinfo=pytz.utc).astimezone(pytz.timezone("Australia/Sydney")) for dt in agg.dtUTC]
        agg['CTD_temp'] = None
        agg['CTD_salt'] = None
        
        # For each swarm locate nearest values
        print('\nCalculating nearest environmental data for swarms..')
        pbar = ProgressBar(max_value=len(agg))
        for idx, row in agg.iterrows():
            # filter survey (isin is much faster than list comprehension])
            dat = df_interp[df_interp['survey'].isin([row.survey])].reset_index(drop=True)
            #Find nearest time
            times = set(dat.dtUTC)
            times = [t for t in times if not pd.isnull(t)]
            tidx = np.argmin([t - row['dtUTC'] for t in times])
            mintime = times[tidx]
            # filter to time
            dat = dat[[t == mintime for t in dat['dtUTC']]].reset_index(drop=True)

            # build KDtree and find nearest point (ignore lat as all the same)
            # Probably not the best as different units but... works really well.. 
            kdt = KDTree(list(zip(dat.lon, dat.depth)))
            vals = dat.iloc[kdt.query([[row.Lon_M, row.Depth_mean]])[1]]
            agg['CTD_temp'].iloc[idx] = float(vals['temp'])
            agg['CTD_salt'].iloc[idx] = float(vals['salt'])
            if (idx%50 == 0) or (idx == len(agg)-1):
                pbar.update(idx)

        return agg


    def plot_xsection(self, survey, transect,
        variable='temperature', vmin=None, vmax=None,
        figsize=(12, 6), set_crossshore_dist=None,
        set_depth=None, interp_limits=(1.875, 5), edge_offset=.1, bathy=True,
        acoustics=False):
        """
        CTD cross-section plot for a transect
        interp limits in km (length) and m (depth)
        """
        # work out how many repeats
        runs = self.transect_runs(survey, transect)
        data = self.data
        meta = self.meta

        print(survey+' '+transect+' '+variable)

        # cmaps
        if variable == 'temperature':
            cmap = cmocean.cm.thermal
            cbar_lab = "Temperature ($^\circ$C)"
        elif variable == 'salinity':
            cmap = cmocean.cm.haline
            cbar_lab = "Salinity (PSU)"

        # plot each run
        fig_out = [None]*len(runs)
        idx = 0
        for run in runs:
            print(idx)
            # get casts
            transect_casts = [list(compress(data, meta['file_name'] == fn))[0] for fn in run[1]['file_name']]
            # Do plotting here...
            meta_section = [df._metadata['File name'] for df in transect_casts]
            meta_section = meta[meta['file_name'].isin(meta_section)]
            # get longitudes
            lons = meta_section['start_longitude']
            # just get one lat as all lats the same more or less (will need when we calculate distance)
            lat = round(meta_section['start_latitude'].mean(), 4)
            # calculate distance
            x = [round(harversine(min(lons), lat, lon, lat), 2) for lon in lons] # Getting lons is not defined error (scoping ?)
            # add distance to each cast
            for i in range(0, len(transect_casts)):
                transect_casts[i]['x'] = x[i]  
            # Create a master data set
            df = pd.concat(transect_casts)    

            # get max depths
            cast_pos = x
            cast_pos[0] = cast_pos[0] + edge_offset
            if cast_pos[-1] > 14:
                cast_pos[-1] = cast_pos[-1] - edge_offset
            max_depths = [math.ceil(max(cast.index.values.astype(float))) for cast in transect_casts]  

            # extract data
            x = np.asanyarray(df.x)
            y = df.index.values.astype(float)
            val = np.asanyarray(df[variable])

            # set plot limits (so all cross sections equal distance)
            maxDist = set_crossshore_dist
            if set_crossshore_dist is None:
                maxDist = x.max()
            maxDepth = set_depth
            if set_depth is None:
                maxDepth = y.max()
            # create interporlations
            xi, yi = np.mgrid[x.min():maxDist:100j, y.min():maxDepth:100j]

            # creat option to choose interporaltion method later
            #val_interp = normal_interp(x, y, val, xi, yi)
            val_interp = rescaled_interp(x, y, val, xi, yi)
            val_interp = np.ma.masked_invalid(val_interp)

            # make a mask for values far away from real values
            for i in range(0, xi.shape[0]):
                for j in range(0, xi.shape[1]):
                    nearest_cast_depth = max_depths[np.argmin(abs(cast_pos - xi[i, j]))]
                    if (xi[i, j] < (x.max() + interp_limits[0])) & (yi[i, j] < (nearest_cast_depth + interp_limits[1])):
                        continue
                    else:
                        val_interp[i, j] = np.nan

            # __Figure__
            fig, ax = plt.subplots(figsize=figsize)
            plt.title(survey+", "+transect+" ("+str(meta_section.iloc[1]['cast_time_local'].date())+')')
            # TOPO (add later from kriging)
            #xm, hm = gen_topomask(h, lon, lat, dx=dx, kind=kind)
            #ax.plot(xm, hm, color="black", linewidth=linewidth, zorder=3)
            #ax.fill_between(xm, hm, y2=hm.max(), color="0.9", zorder=3)
            
            # Station markers
            offset = .5
            ax.plot(cast_pos, [offset]*len(cast_pos), 'wo',  ms=10)
            # show depth of cast
            for j in range(0, len(max_depths)):
                ax.plot([cast_pos[j], cast_pos[j]], [offset, max_depths[j]], c='w', ls='--')
            ax.set_xlabel("Cross-shore distance (km)", fontsize=12)
            ax.set_ylabel("Depth (m)", fontsize=12)
            ax.set_ylim(0, maxDepth)
            ax.set_xlim(0, maxDist)
            ax.invert_yaxis()
            #ax.xaxis.set_ticks_position("top")
            #ax.xaxis.set_label_position("top")
            ax.yaxis.set_ticks_position("left")
            ax.yaxis.set_label_position("left")
            ax.xaxis.set_tick_params(tickdir="out", labelsize=12, pad=1)
            ax.yaxis.set_tick_params(tickdir="out", labelsize=12, pad=1)

            # set vmin and vmax
            if vmin is None:
                vmin = np.nanmin(val_interp)
            if vmax is None:
                vmax = np.nanmax(val_interp)

            #cs = ax.contourf(xi, yi, val_interp, extend="both", zorder=2) # levels=levels1, alpha=1.0,
            cs = ax.pcolor(xi, yi, val_interp, zorder=2, cmap=cmap, vmin=vmin, vmax=vmax)
            cbar = fig.colorbar(cs, ax=ax)#, pad=0.05, shrink=0.95)
            cbar.ax.set_ylabel(cbar_lab, rotation=90, labelpad=16, size=12)
            con = ax.contour(xi, yi, val_interp, 5, colors='k', alpha=.8, linestyles='dashed')
            ax.clabel(con, fontsize=9, inline=True)

            if bathy:
                # load bathymetry data
                # result is a dictionary where keys are the name of objects and the values python
                # objects. In the case of Rds there is only one object with None as key
                result = pyreadr.read_r('./data/surveys/kriging/krig_bathymetry-df.rds')
                bathy_df = result[None]
                # extract nearest lat line
                bathy_df = bathy_df[~np.isnan(bathy_df['depth'])]
                bathy_df = bathy_df[round(bathy_df['lat'],2) == round(lat, 2)]
                bathy_df['lat'] = round(bathy_df['lat'], 2)
                depths = [round(np.mean(df[1].depth)) for df in bathy_df.groupby(round(bathy_df['lon'], 2))]
                bathy_lons = [np.mean(df[1].lon) for df in bathy_df.groupby(round(bathy_df['lon'], 2))]
                minlon = min(meta_section['start_longitude'])
                distance = [round(harversine(minlon, lat, lon, lat), 2) for lon in bathy_lons]
                # make distances negtaive west of first cast
                dis_idx = [i for i in range(len(bathy_lons)) if bathy_lons[i] < minlon]
                if len(dis_idx) > 0:
                    for i in dis_idx:
                        distance[i] = -distance[i]

                ax.fill_between(distance, depths, [1e8]*len(depths), interpolate=True, 
                    color='#808080', zorder=100) 
                ax.plot(distance, depths, color='#4a4a4a', zorder=110)

            if acoustics:
                # load acoustic data
                result = pyreadr.read_r('./data/surveys/acoustics/agg.rds')
                agg = result[None]
                agg = agg[agg['survey'] == survey]
                agg.dt = [pytz.utc.localize(datetime.strptime(str(d)+'_'+str(t).replace(' ',''), "%Y%m%d_%H:%M:%S.%f")) for d, t in zip(agg['Date_M'], agg['Time_M'])]
                agg.dt_local = [dt.replace(tzinfo=pytz.utc).astimezone(pytz.timezone("Australia/Sydney")) for dt in agg.dt]
                # filter to date
                agg = agg[[a == b for a, b in zip([dtl.date() for dtl in  agg.dt_local], [meta_section.iloc[1]['cast_time_local'].date()]*len(agg))]]

                agg = agg[round(agg['Lat_M'],2) == round(lat, 2)]
                agg['Lat_M'] = round(agg['Lat_M'], 2)
                agg = agg.reset_index(drop=True)
                if len(agg) > 0:
                    distance = [round(harversine(minlon, lat, lon, lat), 2) for lon in agg.Lon_M]
                    # make distances negtaive west of first cast
                    dis_idx = [i for i in range(0, len(agg)) if agg['Lon_M'][i] < minlon]
                    if len(dis_idx) > 0:
                        for i in dis_idx:
                            distance[i] = -distance[i]
                    agg_size = (((agg.Sv_mean - min(agg.Sv_mean))/(max(agg.Sv_mean) - min(agg.Sv_mean)))*50)+15
                    ax.scatter(distance, agg.Depth_mean, s=agg_size, 
                        facecolors='none', edgecolor='#4287f5', zorder=150)
                else:
                    print('No aggrgations for transect...')

            fig_out[idx] = fig
            idx += 1

        return fig_out


    def plot_transect_map(self, survey, transect,
                          inset=True, inset_loc=4,
                          margin=0.2, figsize=(6, 6)):
        # work out how many repeats
        runs = self.transect_runs(survey, transect)

        # get just first run (all maps the same)n
        transect_casts = [list(compress(self.data, self.meta['file_name'] == fn))[0] for fn in runs[0][1]['file_name']]

        # get lon/lat for casts
        lonlat = [(df._metadata['Start longitude'], df._metadata['Start latitude']) for df in transect_casts]
        # extent = [minlon, maxlon, minlat, maxlat]
        extent = [min([x[0] for x in lonlat]) - margin, max([x[0] for x in lonlat]) + margin,
                  min([x[1] for x in lonlat]) - margin, max([x[1] for x in lonlat]) + margin]

        fig, ax = plt.subplots(figsize=figsize)
        dateset = ', '.join(set([str(run[1].iloc[1]['cast_time_local'].date()) for run in runs]))
        plt.title(survey+", "+transect+" ("+dateset+')')
        m = Basemap(llcrnrlon=extent[0], urcrnrlon=extent[1],
            llcrnrlat=extent[2], urcrnrlat=extent[3],
            projection='merc', resolution='h')
        m.drawcoastlines()
        m.fillcontinents(color='0.85')
        meridians = np.arange(extent[0], extent[1] + 1, round((extent[1] - extent[0]) / 4, 2))
        parallels = np.arange(extent[2], extent[3] + 1, round((extent[3] - extent[2]) / 4, 2))
        m.drawparallels(parallels, linewidth=0, labels=[1, 0, 0, 0])
        m.drawmeridians(meridians, linewidth=0, labels=[0, 0, 0, 1])
        m.ax = ax

        if inset:
            axin = inset_axes(m.ax, width="30%", height="30%", loc=inset_loc)
            # Global inset map.
            inmap = Basemap(lon_0=np.mean(m.boundarylons),
                            lat_0=np.mean(m.boundarylats),
                            projection='ortho', ax=axin, anchor='NE')
            inmap.drawcountries(color='white')
            inmap.fillcontinents(color='gray')
            bx, by = inmap(m.boundarylons, m.boundarylats)
            xy = list(zip(bx, by))
            mapboundary = Polygon(xy, edgecolor='k', linewidth=1, fill=False)
            inmap.ax.add_patch(mapboundary)

        # limits
        limits = [min(m.boundarylons), max(m.boundarylons), min(m.boundarylats), max(m.boundarylats)]

        # Add stations and topo
        x, y, topo = etopo_subset(limits, smoo=True, tfile=None)
        topo = np.where(topo > -1., 1.e10, topo)
        topo = np.ma.masked_values(topo, 1.e10)
        cs = m.contour(x, y, -topo, (100, 200, 500, 1000), colors='k',
                        latlon=True, alpha=0.5)
        m.ax.clabel(cs, fmt='%1.0f m', fontsize=8, inline=1)
        m.plot([x[0] for x in lonlat], [x[1] for x in lonlat], 'k.', latlon=True)

        return fig

# End Cast Class definitions #
#----------------------------#
def kdtree_process(kdt, lon0, lat0, array_shape):
    """
    Adapted from:
    https://github.com/Unidata/python-workshop/blob/fall-2016/notebooks/netcdf-by-coordinates.ipynb
    """
    lat0_rad = lat0 * math.pi/180.0
    lon0_rad = lon0 * math.pi/180.0
    clat0,clon0 = np.cos(lat0_rad), np.cos(lon0_rad)
    slat0,slon0 = np.sin(lat0_rad), np.sin(lon0_rad)
    dist_sq_min, minindex_1d = kdt.query([clat0*clon0, clat0*slon0, slat0])
    iy_min, ix_min = np.unravel_index(minindex_1d, array_shape)
    
    return (iy_min, ix_min)

def read_CastAway_csv(dir, metafile=None):
    """
    Read a set of CastAway casts
    Pass directory of casts
    Also can pass metafile of casts to read
    """
    # Read in cast meta data
    if metafile is not None:
        meta = pd.read_csv(metafile)
        data = list(meta.file_name)
    else:
        raise ValueError("Haven't coded this in yet!")
        """
        read all files names in directory
        """

    for i in range(0, len(data)):
        data[i] = ctd.read.from_castaway_csv(dir+'/'+data[i]+'.csv')

    casts = Casts(data, meta)

    return casts

def harversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    # harversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2.)**2. + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2.)**2.
    c = 2. * math.asin(math.sqrt(a))
    km = 6371. * c # radius of earth
    return km


def reverse_harversine(lon_orig, lat_orig, km, direction='east'):
    """
    Calculate points directly north, south, east and 
    west a certain distance from given coordinates
    """
    if direction != 'east':
        raise ValueError('Only east is supported right now. Soz.')
    # convert decimal degrees to radians
    lon_orig, lat_orig = map(math.radians, [lon_orig, lat_orig])
    # reverse harversine formula
    c = km / 6371.
    a = math.sin(c/2.)**2.
    dlat = 2. * math.asin(math.sqrt(a))
    dlon = 2. * math.asin(math.sqrt(a/(math.cos(lat_orig)**2.)))
    # convert back to decimal degrees 
    lon_orig, lat_orig, dlat, dlon = map(math.degrees, [lon_orig, lat_orig, dlat, dlon])
    # find coordinates
    # north = lat_orig + dlat
    # south = lat_orig - dlat
    east = lon_orig + dlon
    # west = lon_orig - dlon
    # correct over the 0-360 degree line
    # if west > 360:
    #     west = west - 360
    if east > 360:
        east = east - 360

    # export region
    return east

def normal_interp(x, y, a, xi, yi):
    rbf = scipy.interpolate.Rbf(x, y, a)
    ai = rbf(xi, yi)
    return ai

def rescaled_interp(x, y, a, xi, yi):
    a_rescaled = (a - a.min()) / a.ptp()
    ai = normal_interp(x, y, a_rescaled, xi, yi)
    ai = a.ptp() * ai + a.min()
    return ai

def plot(x, y, a, ai, title):
    fig, ax = plt.subplots()

    im = ax.imshow(ai.T, origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
    ax.scatter(x, y, c=a)

    ax.set(xlabel='X', ylabel='Y', title=title)
    fig.colorbar(im)









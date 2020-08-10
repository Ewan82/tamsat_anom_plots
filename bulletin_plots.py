import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import netCDF4 as nc
import datetime as dt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from shapely.geometry.polygon import LinearRing
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.image as mpimg


def plt_pcolor_map(var_arr, lon, lat,
                   label=r'Precipitation Anomaly (mm day$^{-1}$)', title='', vmin=None, vmax=None, col_map='RdBu',
                   levels=np.linspace(-100, 120, 13), extend='max', plt_type='contour'):
    """
    Function that will plot data on a map, given data array and lat/lon values.
    :param var_arr: 2D array of data (arr)
    :param lon: longitude values (arr)
    :param lat: latitude values (arr)
    :param label: colorbar label (str)
    :param title: plot title (str)
    :param vmin: min value for colorbar (flt)
    :param vmax: max value for colorbar (flt)
    :param col_map: Python colormap to use (str)
    :param levels: levels to discretise colorbar (arr)
    :param extend: extend colorbar in this direction (str)
    :param plt_type: typr of plot, either 'contour' or 'mesh' (str)
    :return: figure and ax objects
    """
    sns.set_context('paper')
    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(8,12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    # create colorbar
    cmap = plt.cm.get_cmap(col_map)
    norm = colors.BoundaryNorm(levels, cmap.N)
    # plot either contour of pcolormesh map
    if plt_type=='contour':
        im = ax.contourf(lon, lat,
                         var_arr,
                         levels=levels,
                         norm=norm,
                         transform=ccrs.PlateCarree(),
                         spacing='uniform',
                         extend=extend,
                         vmin=vmin, vmax=vmax,
                         cmap=col_map)
    elif plt_type=='mesh':
        im = ax.pcolormesh(lon, lat,
                         var_arr,
                         transform=ccrs.PlateCarree(),
                         norm=norm,
                         vmin=vmin, vmax=vmax,
                         cmap=col_map)
    # add features to map, borders, counties, etc
    #ax.add_feature(cfeature.OCEAN)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.STATES, edgecolor = 'gray', alpha=0.3)
    ax.add_feature(cfeature.BORDERS)
    #ax.add_feature(cfeature.RIVERS, edgecolor='black')
    #ax.add_feature(cfeature.LAKES, edgecolor='red', facecolor = 'none')
    #rivers_lake_centerlines = cfeature.NaturalEarthFeature('physical', 'rivers', '50m', edgecolor = 'red')
    #rivers_lake_centerlines = cfeature.NaturalEarthFeature(category = 'physical', name = 'rivers_lake_centerlines',
    #scale = '10m', facecolor = 'none', edgecolor = 'black')
    #ax.add_feature(rivers_lake_centerlines)
    lake = cfeature.NaturalEarthFeature(category = 'physical', name = 'lakes',
    scale = '10m', facecolor = 'none', edgecolor = 'black')
    ax.add_feature(lake)
    # add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=plt.Axes)
    cbar = fig.colorbar(im, cax=cax, extend=extend)
    cbar.set_label(label, fontsize=12)
    # add title
    ax.set_title(title)
    # add gridlines
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    ax.xaxis.set_major_formatter(lon_formatter)
    lat_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=1.0, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent([min(lon), max(lon), min(lat), max(lat)],
                   ccrs.PlateCarree())
    return fig, ax


def plot_tamsat_anom_ghana(year=2020, month=6, lat_bnds=(3.75, 11.9), lon_bnds=(-4.93, 2.95)):
    """
    Plot map of monthly TAMSAT rainfall anomaly using data on JASMIN GWS
    :param year: year of data (flt)
    :param month: month of data(flt)
    :param lat_bnds: upper/lower bounds for latitude (tup)
    :param lon_bnds: upper/lower bounds for longitude (tup)
    :return: figure and ax objects
    """
    file_str = '/gws/nopw/j04/tamsat/public/rfe/data/v3.1/monthly-anomalies/{}/{}' \
               '/rfe{}_{}_anom.v3.1.nc'.format(year, str(month).zfill(2), year, str(month).zfill(2))
    dat = nc.Dataset(file_str, 'r')
    lats = dat.variables['lat'][:]
    lons = dat.variables['lon'][:]
    rfe = dat.variables['rfe_filled']
    lat_idx = np.where((lats > lat_bnds[0]) & (lats < lat_bnds[1]))
    lon_idx = np.where((lons > lon_bnds[0]) & (lons < lon_bnds[1]))
    fig, ax = plt_pcolor_map(rfe[0, lat_idx[0], lon_idx[0]], lons[lon_idx[0]], lats[lat_idx[0]],
                             label=r'Precipitation Anomaly (mm)',
                             vmin=-175, vmax=175, title='TAMSAT Monthly Rainfall Anomlay {}/{}'.format(str(month).zfill(2), year),
                             col_map='RdBu',levels=np.linspace(-200,200,17), extend='both', plt_type='contour')
    add_logo(fig)
    dat.close()
    return fig, ax


def add_logo(fig, x_o=400, y_o=25):
    """
    Function that adds GSSTI and NCEO logo to a figure
    :param fig: Figure object to add logo to
    :param x_o: xo position of logo on figure (float)
    :param y_o: yo position of logo on figure (float)
    :return: 'logo added' (str)
    """
    logo_arr = mpimg.imread('gssti_nceo_logo2.png')
    fig.figimage(logo_arr, xo=x_o, yo=y_o)
    return 'logo added'
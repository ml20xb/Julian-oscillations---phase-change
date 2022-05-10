import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.ticker as cticker
import cartopy.io.shapereader as shapereader

# 绘制地图的函数
def contour_map(fig,img_extent):
    lon_W,lon_E,lat_S,lat_N = img_extent
    fig.set_extent(img_extent, crs=ccrs.PlateCarree())
    fig.set_xticks(np.arange(lon_W,lon_E+2,2), crs=ccrs.PlateCarree())
    fig.set_yticks(np.arange(lat_S,lat_N+2,2), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    fig.xaxis.set_major_formatter(lon_formatter)
    fig.yaxis.set_major_formatter(lat_formatter)
    China_map = shapereader.Reader("D:/Specwork/python_plot_data/shp/cnmap/cnhimap.shp").geometries()
    fig.add_geometries(China_map, ccrs.PlateCarree(),facecolor='none', edgecolor='black',zorder = 1,lw=0.7)

main_dir = "D:/Specwork/Num01-50/Num35_MJO/20220322/"

names = ['AGR','DST','INB','INP','NRD','ORD','POW','RES','SHP','SOL']

ff = xr.open_dataset(main_dir+'GRIDCRO2D.nc')
lat = ff['LAT'][0,0,:,:]
lon = ff['LON'][0,0,:,:]

for n in range(len(names)):
	name = names[n]
	f = xr.open_dataset(main_dir + 'yrd_ei_2017/gridded_emission/GEI-YRD-2017-V2.0-'+name+'-TOT.ncf')

	NOX = f['NOX'][0,0,:,:]

	if n==0:
		NOX_sum = NOX 
	else:
		NOX_sum = NOX_sum + NOX

latS = 26.5
latN = 35.5
lonW = 114.5
lonE = 123.
plt.figure(figsize=(8,6))
ax = plt.subplot(1,1,1,projection=ccrs.PlateCarree())
contour_map(ax,[lonW,lonE,latS,latN])
ax.set_ylim(lat[:,0].values.min(),lat[:,0].values.max())
ax.set_xlim(lon[0,:].values.min(),lon[0,:].values.max())
c=ax.pcolormesh(lon[0,:],lat[:,0],NOX_sum,vmin=0,vmax=1500,cmap='Oranges',transform=ccrs.PlateCarree())
cb=plt.colorbar(c)
plt.savefig(main_dir+"out.pdf")
#!/usr/bin/env python3

import numpy as np
import julian
import datetime
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import aacgm
import argparse
from matplotlib.colors import Normalize

from global_vel_class import global_vel
from global_pot_class import global_pot
from vel_col_cmap import vel_col
from noonlon import geo_noonlon
from mag_continents import magContinents
from map_utils import get_local_azm
from date_strings import make_date_str

Re=6357e3
dtor=np.pi/180

def color_contours(fig,ax,lats,lons,pot,cmap,norm,mxval=30000):
    dlat=0
    alpha=.7
    j=0
    num=len(pot)

    while dlat ==0:
        dlat=lats[j+1]-lats[j]
        j+=1

        
    patches=[]
    colors=np.zeros((num))
    colors=[]


    dlon=lons[1]-lons[0]
    dlon_last=dlon
    
    for j in range(1,num):
        
        if np.abs(lats[j]-lats[j-1])<dlat:
            dlon=lons[j]-lons[j-1]

        x=lons[j]-dlon/2
        y=lats[j]-dlat/2
        width=dlon
        # if width > 20:
        #     print('global: ',y,x,width)

        c=cmap(norm(pot[j]))
        patches.append(Rectangle((x,y),width,dlat,color=c,alpha=alpha,ec=None))

        
        
    p=PatchCollection(patches,transform=ccrs.PlateCarree(),cmap=cmap,match_original=True)
    ax.add_collection(p)

    cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
    posn=ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0+.1*posn.height,
                      0.03, posn.height*.8])

    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbar_ax)





def read_pot(fp):

    dateline=fp.readline()
    dtstrs=dateline.split()
    yr=int(dtstrs[0].strip())
    mo=int(dtstrs[1].strip())
    dy=int(dtstrs[2].strip())
    hr=int(dtstrs[3].strip())
    mt=int(dtstrs[4].strip())
    sc=int(dtstrs[5].strip())

    zline=fp.readline()
    npts=int(fp.readline().strip())

    lats=np.zeros(npts)    
    lons=np.zeros(npts)
    pot=np.zeros(npts)
    ee=np.zeros(npts)
    en=np.zeros(npts)
    ve=np.zeros(npts)
    vn=np.zeros(npts)

    for j in range(npts):
        pline=fp.readline()
        pstrs=pline.split()
        lats[j]=float(pstrs[0].strip())
        lons[j]=float(pstrs[1].strip())
        pot[j]=float(pstrs[2].strip())
        ee[j]=float(pstrs[3].strip())
        en[j]=float(pstrs[4].strip())
        ve[j]=float(pstrs[5].strip())
        vn[j]=float(pstrs[6].strip())

    return(lats,lons,pot,ee,en,ve,vn,yr,mo,dy,hr,mt,sc)




colors=[]
v_col=vel_col(ncol=28)
mxval=30000

parser = argparse.ArgumentParser(description='local_df_vel argument parser')
parser.add_argument("pot_file", nargs=1, type=str,
                    help="file must be specified (eg: '201001200000-201001202359.map.pot')")

parser.add_argument('--mxval', required=False, type=int, default=mxval)
parser.add_argument('--local', required=False, dest='local', action='store_true', default=False)
parser.add_argument('--geo', required=False, dest='geo', action='store_true', default=False)
parser.add_argument('--png', required=False, dest='png', action='store_true', default=False)
parser.add_argument('--south', required=False, dest='south', action='store_true', default=False)



args = parser.parse_args()
pot_file = args.pot_file[0]
mxval=args.mxval
local=args.local
geo=args.geo
pngplt=args.png
south=args.south

fp=open(pot_file,"r")

if pngplt:
    xdm=12
    ydm=12
else:
    xdm=6
    ydm=6
    
fig=plt.figure(figsize=(xdm,ydm))

istrs=["00","01","02","03","04","05","06","07","08","09"]

while(True):
    fig.clf()
    ax = fig.add_subplot(1,1,1)

    lats,lons,pot,ee,en,ve,vn,yr,mo,day,hr,mt,sc=read_pot(fp)
    npts=len(lats)

    norm=Normalize(vmax=mxval,vmin=-mxval)

    min_lon=0.
    max_lon=359.
    if south:
        max_lat=-55.
        min_lat=-90
        pole_lat=-90.    
    else:    
        min_lat=55.
        max_lat=90.
        pole_lat=90

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)


    aacgm.set_datetime(yr,mo,day,hr,mt,sc)

    if geo:
        tme=datetime.datetime(yr,mo,day,hr,mt,sc)
        noonlon=geo_noonlon(tme)-180
    else:
        if south:
            noonlon=aacgm.inv_mlt_convert(yr,mo,day,hr,mt,sc,12)
        else:
            noonlon=aacgm.inv_mlt_convert(yr,mo,day,hr,mt,sc,12) - 180
    
    
    
    crs=ccrs.AzimuthalEquidistant(central_longitude=noonlon,
                                  central_latitude=pole_lat,
                                  false_easting=0.0,
                                  false_northing=0.0,
                                  globe=None)

    ax = fig.add_subplot(1,1,1,projection=crs)
    ax.set_extent((min_lon,max_lon,min_lat,max_lat),ccrs.PlateCarree())
    # if geo:
    #     ax.add_feature(cfeature.OCEAN)
    #     ax.add_feature(cfeature.LAND, edgecolor='black')
    # else:
    #     magContinents(ax,boundary=True,alpha=.3)
    
    ax.gridlines()    
    ax.set_boundary(circle, transform=ax.transAxes)
    
    color_contours(fig,ax,lats,lons,pot,v_col,norm,mxval=mxval)

    dlat=0
    j=0
    while dlat ==0:
        dlat=lats[j+1]-lats[j]
        dy=dlat*dtor*Re
        j+=1
    

    dlon=lons[1]-lons[0]
    vlat=np.zeros(npts)
    vlon=np.zeros(npts)
    ev=np.zeros(npts)
    vv=np.zeros(npts)
    angs=np.zeros(npts)
    vangs=np.zeros(npts)
    evv=np.zeros(npts)


    escl=1
    vscl=1/15000

    # for j in range(npts):
    
    #     mag=np.sqrt(ee[j]*ee[j]+en[j]*en[j])*escl
    #     azm=np.arctan2(ee[j],en[j])/dtor
    
    #     ang=get_local_azm(ax,lats[j],lons[j],azm)        
    
    #     angs[j]=ang
    #     ev[j]=mag
    
    # ax.quiver(lons,lats,ev,evv,width=.002,zorder=2,scale=1,angles=angs,scale_units='width',
    #               transform=ccrs.PlateCarree(),headwidth=2,clip_on=False)
    
    for j in range(npts):
        
        mag=np.sqrt(ve[j]*ve[j]+vn[j]*vn[j])*vscl
        azm=np.arctan2(ve[j],vn[j])/dtor
        
        ang=get_local_azm(ax,lats[j],lons[j],azm)        
        
        vangs[j]=ang
        vv[j]=mag

    ax.quiver(lons,lats,vv,evv,width=.002,zorder=2,scale=1,angles=vangs,scale_units='width',
              transform=ccrs.PlateCarree(),headwidth=2,clip_on=False)



    print(yr,mo,day,hr,mt)

    

    if mo < 10: mostr=istrs[mo]
    else: mostr=str(mo)

    if dy < 10: dystr=istrs[day]
    else: dystr=str(day)

    if hr < 10: hrstr=istrs[hr]
    else: hrstr=str(hr)

    if mt < 10: mtstr=istrs[mt]
    else: mtstr=str(mt)
    

    figname=str(yr)+mostr+dystr+hrstr+mtstr+"_pot.png"
    
    plt.savefig(figname,orientation='landscape',bbox_inches='tight',dpi=150)



plt.show()

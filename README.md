#%%cargo librerias
import os 
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
#%%
#SIMULACION 1

#%%cargo variables y seteo ruta
#dir = "/home/nadia/Documentos/circulaciongeneral/practica/p1atm/eb1/" #ruta nanu
dir = 'C:/Users/Alicia/Desktop/circulacion/Practica_atmo/TPRossby/SM1/' # ruta celu 
os.chdir(dir)
dS = xr.open_dataset(dir+"shallow.nc",decode_times=False)
#%% cargo variables
lon=dS["lon"].values
lat=dS["lat"].values
forz=dS["fr"].values
time=dS["time"].values
u=dS["ucomp"].values
v=dS["vcomp"].values
vort=dS["vor"].values
div=dS["div"].values
lons, lats = np.meshgrid(lon, lat) #paso lon y lat a matrices. 
#%%calculo ks
from dy import derivy #cargar la funcion del campus
Uyy=derivy(vort[-1,:,:],128/180*111320) #diferencial y= numero grillas a grados y despues a metros
betha=2*7.29e-5*np.cos(lat*np.pi/180)/6371000 #2*omega*coslat/rt
betha=np.array([betha,]*256).transpose()  #celes fijate que esta matriz queda simetrica respecto al ecuador y toda positiva, eso es lo que queremos
grad_mer_vortabs=betha-Uyy
#por si hay valores negativos de grad_mer_vortabs/u
ks=np.empty((128,256))
ks[:]=np.nan
for i in range(0,128):
    for j in range(0,256):
        if (grad_mer_vortabs[i,j]/u[-1,i,j]>0):
            ks[i,j]=((np.sqrt((grad_mer_vortabs[i,j])/u[-1,i,j]))) 
#%%Obtenesmo Ks transformado en cantidad de ondas por circulo de latitud
Ks=6371000*ks*np.cos(lats*np.pi/180)          
#%%Calculamos los minimos y maximos de cada uno
forz_levs=np.linspace(np.min(forz),np.max(forz),50)
vort_levs=np.linspace(np.min(vort[-1]),np.max(vort[-1]),11)
u_levs=np.linspace(np.min(u[-1]),np.max(u[-1]),11)
grad_mer_vortabs_levs=np.linspace(np.min(grad_mer_vortabs),np.max(grad_mer_vortabs),20)
Ks_levs=np.linspace(np.nanmin(Ks),np.nanmax(Ks),20) 
Uyy_levs =np.linspace(np.nanmin(Uyy),np.nanmax(Uyy),20)        
DATOS=[forz[-1]*1e-4,u[-1],grad_mer_vortabs*1e11,Ks,Uyy,vort[-1]]
LEVEL=[forz_levs*1e-4,u_levs,grad_mer_vortabs_levs*1e11,Ks_levs,Uyy_levs,vort_levs]
Titulo=["Forzante EB1","Viento zonal EB1","Gradiente meridional de vorticidad absoluta EB1","Ks$'$ EB1","Uyy","vort"]
name=['eb1_forzante.jpg','eb1_u.jpg','eb1_grad_vortabs.jpg','eb1_Ks.jpg','UyyEb1.jpg','vortEb1.jpg']
unid=["*1e4 m^2*s^-2","m/s","*1e-11 m-1*s-1","adim.","m^-1*s^-1","ver"]
ch_cmap=["Reds","Reds","Reds","Reds","jet","jet"]
#%% Creamos las figuras
for i in range(0,6):
    fig=plt.figure(figsize=(8,5),dpi=200)

    #Definimos proyección
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0, 359, -90, 90], crs=crs_latlon)
    #Graficamos
    im=ax.contourf(lons,lats,DATOS[i],LEVEL[i],cmap=plt.get_cmap(ch_cmap[i]),extend='both',transform=crs_latlon)

    #Agregamos barra de colores
    cbar=plt.colorbar(im,fraction=0.052, pad=0.04,shrink=0.8,aspect=8,format='%0.2f')
    cbar.ax.set_title(unid[i],fontsize=10) #chequear unidades
    #Características del mapa
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360,40), crs=crs_latlon)
    ax.set_yticks(np.arange(-90, 95,15), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #Características del mapa
    #Titulo
    plt.title(Titulo[i],fontsize=16, y=0.98,loc="center")
    #Guardar figura
    fig.tight_layout()
    plt.savefig(name[i])
#%%zoom de los hemisferios
Ks_level_zoomHN=np.linspace(np.min(Ks[(15,90),:]), np.max(Ks[(15,90),:]),11) #HN   
Ks_level_zoomHS=np.linspace(np.min(Ks[(-90,-15),:]), np.max(Ks[(-90,-15),:]),11) #HS
HEM=[Ks_level_zoomHN,Ks_level_zoomHS]
TITH=["Ks$'$ EB1 (HN)","Ks$'$ EB1 (HS)"]
nameH=["eb1_ks_HN.jpg","eb1_ks_HS.jpg"]
r1=[15,-90]
r2=[90,-15]
r3=[105,0]
for j in range(0,2):
    fig=plt.figure(figsize=(8,5),dpi=200)
    #Definimos proyección
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0, 359, r1[j], r2[j]], crs=crs_latlon)
    #Graficamos
    im=ax.contourf(lons,lats,Ks,HEM[j],cmap=plt.get_cmap("Reds"),extend='both',transform=crs_latlon)
    #Agregamos barra de colores
    cbar=plt.colorbar(im,fraction=0.052, pad=0.04,shrink=0.8,aspect=8,format='%0.2f')
    cbar.ax.set_title("adim.",fontsize=10) #chequear unidades
    #Características del mapa
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360,40), crs=crs_latlon)
    ax.set_yticks(np.arange(r1[j], r3[j],15), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #Características del mapa
    #Titulo
    plt.title(TITH[j],fontsize=16, y=0.98,loc="center")
    #Guardar figura
    fig.tight_layout()
    plt.savefig(nameH[j])

#%%fin simulacion 1
    
#%%    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"                                                                        "
"                                SIMULACION 2                            "
"                                                                        "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
#%%cargo librerias
import os 
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
#%%cargo variables y seteo ruta
#dir = "/home/nadia/Documentos/circulaciongeneral/practica/p1atm/eb2/" #ruta nanu
dir = 'C:/Users/Alicia/Desktop/circulacion/Practica_atmo/TPRossby/SM2/' #ruta celeste
os.chdir(dir)
dS = xr.open_dataset(dir+"shallow.nc",decode_times=False)
#%% cargo variables
lon=dS["lon"].values
lat=dS["lat"].values
forz=dS["fr"].values
time=dS["time"].values
u=dS["ucomp"].values
v=dS["vcomp"].values
vort=dS["vor"].values
div=dS["div"].values
lons, lats = np.meshgrid(lon, lat) #paso lon y lat a matrices. 
#%%calculo ks
from dy import derivy #cargar la funcion del campus
Uyy=derivy(vort[-1],128/180*111320) #diferencial y= numero grillas a grados y despues a metros
betha=2*7.29e-5*np.cos(lat*np.pi/180)/6371000 #2*omega*coslat/rt
betha=np.array([betha,]*256).transpose()  #celes fijate que esta matriz queda simetrica respecto al ecuador y toda positiva, eso es lo que queremos
grad_mer_vortabs=betha-Uyy
#por si hay valores negativos de grad_mer_vortabs/u
ks=np.empty((128,256))
ks[:]=np.nan
for i in range(0,128):
    for j in range(0,256):
        if (grad_mer_vortabs[i,j]/u[-1,i,j]>0):
            ks[i,j]=((np.sqrt((grad_mer_vortabs[i,j])/u[-1,i,j]))) 
#%%Obtenesmo Ks transformado en cantidad de ondas por circulo de latitud
Ks=6371000*ks*np.cos(lats*np.pi/180)          
#%%Calculamos los minimos y maximos de cada uno
forz_levs=np.linspace(np.min(forz),np.max(forz),50)
u_levs=np.linspace(-np.max(u[-1]),np.max(u[-1]),11)
grad_mer_vortabs_levs=np.linspace(np.min(grad_mer_vortabs),np.max(grad_mer_vortabs),20)
Ks_levs=np.linspace(np.nanmin(Ks),np.nanmax(Ks),20)          
DATOS=[forz[-1]*1e-4,u[-1],grad_mer_vortabs*1e11,Ks]
LEVEL=[forz_levs*1e-4,u_levs,grad_mer_vortabs_levs*1e11,Ks_levs]
Titulo=["Forzante EB2","Viento zonal EB2","Gradiente meridional de vorticidad absoluta EB2","Ks$'$ EB2"]
name=['eb2_forzante.jpg','eb2_u.jpg','eb2_grad_vortabs.jpg','eb2_Ks.jpg']
unid=["*1e4 m^2*s^-2","m/s","*1e-11 m-1*s-1","adim."]
ch_cmap=["Reds","seismic","Reds",'PiYG']
#%% Creamos las figuras
for i in range(0,4):
    fig=plt.figure(figsize=(8,5),dpi=200)

    #Definimos proyección
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0, 359, -90, 90], crs=crs_latlon)
    #Graficamos
    im=ax.contourf(lons,lats,DATOS[i],LEVEL[i],cmap=plt.get_cmap(ch_cmap[i]),extend='both',transform=crs_latlon)
    if i==1:
        ax.contour(lons,lats,DATOS[i],levels=0,colors="white",transform=crs_latlon)

    #Agregamos barra de colores
    cbar=plt.colorbar(im,fraction=0.052, pad=0.04,shrink=0.8,aspect=8,format='%0.2f')
    cbar.ax.set_title(unid[i],fontsize=10) #chequear unidades
    #Características del mapa
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360,40), crs=crs_latlon)
    ax.set_yticks(np.arange(-90, 95,15), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #Características del mapa
    #Titulo
    plt.title(Titulo[i],fontsize=16, y=0.98,loc="center")
    #Guardar figura
    fig.tight_layout()
    plt.savefig(name[i])

#%%Hago zoom sobre figura de Ks para ambos hemisderios
#creamos las listas
Ks_level_zoomHN=np.linspace(np.min(Ks[(15,90),:]), np.max(Ks[(15,90),:]),11) #HN   
Ks_level_zoomHS=np.linspace(np.min(Ks[(-90,-15),:]), np.max(Ks[(-90,-15),:]),11) #HS
r1=[15,-90]
r2=[90,-15]
r3=[105,0]
HEM=[Ks_level_zoomHN,Ks_level_zoomHS]
TITH=["Ks$'$ EB2 (HN)","Ks$'$ EB2 (HS)"]
nameH=["eb2_ks_HN.jpg","eb2_ks_HS.jpg"]
#ploteamos
for j in range(0,2):
    fig=plt.figure(figsize=(8,5),dpi=200)
    #Definimos proyección
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0, 359, r1[j],r2[j]], crs=crs_latlon)
    # Graficamos
    im=ax.contourf(lons,lats,Ks,HEM[j],cmap=plt.get_cmap("Reds"),extend='both',transform=crs_latlon)
    #Agregamos barra de colores
    cbar=plt.colorbar(im,fraction=0.052, pad=0.04,shrink=0.8,aspect=8,format='%0.2f')
    cbar.ax.set_title("adim.",fontsize=10) #chequear unidades
    #Características del mapa
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360,40), crs=crs_latlon)
    ax.set_yticks(np.arange(r1[j], r3[j],15), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #Características del mapa
    #Titulo
    plt.title(TITH[j],fontsize=16, y=0.98,loc="center")
    #Guardar figura
    fig.tight_layout()
    plt.savefig(nameH[j])
    
#%% Vemos en una longitud donde se puede detectar las latitudes criticas de la onda de Rossby

Kes=Ks[:,93]
u1=u[-1]
u1_r=u1[:,93]
lat1=lat[63:]
K1=Kes[63:]
u2=u1_r[63:]
plt.plot(lat,Kes)
fig=plt.figure()
plt.subplot(121)
plt.plot(lat1,K1,'m')
plt.xlabel('Latitud en Grados',fontsize=12)
plt.ylabel('Ks$´$(adm)',fontsize=12)
plt.xticks(visible=True,fontsize=12)
plt.yticks(visible=True,fontsize=12)
plt.grid(True)
plt.title('propagación de la onda 40°E(HN)')
plt.subplot(122)
plt.plot(lat1,u1_r[63:],'r')
plt.xlabel('Latitud en Grados',fontsize=12)
plt.ylabel('Viento zonal(m/s)',fontsize=12)
plt.title('Viento zonal a 40°E(HN)',fontsize=12)
plt.xticks(visible=True,fontsize=12)
plt.yticks(visible=True,fontsize=12)
plt.grid(True)
fig.tight_layout()
plt.savefig('propagacion de la onda de ROssby a 40°E.jpg',bbox_inches='tight',pad_inches=0.05,dpi=200)
    
    
    
    

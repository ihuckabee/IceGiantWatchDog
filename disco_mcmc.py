import numpy as np
import os
#os.environ["picaso_refdata"] = "/home/izzyh/Documents/picaso-master/reference"
#os.environ["PYSYN_CDBS"] = "/home/izzyh/Documents/picaso-master/grp/redcat/trds"
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
refdata = os.getenv("picaso_refdata")
import sys
sys.path.append("..")
import pdb
from virga import justdoit as vdi
import pandas as pd
import astropy.units as u
from picaso import justplotit as ppi
from picaso import justdoit as pdi
import pandas as pd
import xarray as xr
import time
from bokeh.plotting import show, figure
import pickle
import math
import copy
import matplotlib.pyplot as plt

def lc_forMCMC(spots_dict): #where spots are a tuple of lat_c, lon_c, dlat, dlon, spot_type
    #unpacking dict 
    #picaso planet object 
    neptune = spots_dict['picaso_object']['planet']
    opacity = spots_dict['picaso_object']['opacity']

    #rotation parameters
    p_rot = spots_dict['rot_params']['p_rot']
    timestep = spots_dict['rot_params']['timestep']
    timedur = p_rot * spots_dict['rot_params']['num_rot']
    rnep = spots_dict['rot_params']['rnep']

    #spatial grid
    lon = spots_dict['spatial_grid']['lon_array']
    lat = spots_dict['spatial_grid']['lat_array']
    lon_len = len(lon)
    lat_len = len(lat)

    #vertical profile
    pres = spots_dict['vertprofile']['pressure']
    temp = spots_dict['vertprofile']['temp_tile']
    chem_ch4 = spots_dict['vertprofile']['CH4']
    chem_h2s = spots_dict['vertprofile']['H2S']
    chem_h2 = spots_dict['vertprofile']['H2']
    chem_he = spots_dict['vertprofile']['He']
    nwno = spots_dict['vertprofile']['nwno']

    #wind profile
    xspeed = spots_dict['windprofile']['xspeed']
    ylat = spots_dict['windprofile']['ylat']

    #cloud_info 
    spots = spots_dict['cloud_info']['spots'] 
    df_base = spots_dict['cloud_info']['cloud_base']
    df_4 = spots_dict['cloud_info']['cloud_4']
    cloud_wnogrid = spots_dict['cloud_info']['wnogrid']

    opd_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); g0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); w0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno))
    test = np.zeros((len(lon), len(lat), len(pres)-1, nwno))

    opd_base = np.array(df_base['opd']).reshape(len(pres)-1,nwno)
    g0_base= np.array(df_base['g0']).reshape(len(pres)-1,nwno)
    w0_base = np.array(df_base['w0']).reshape(len(pres)-1,nwno)

    opd_4 = np.array(df_4['opd']).reshape(len(pres)-1,nwno)
    g0_4 = np.array(df_4['g0']).reshape(len(pres)-1,nwno)
    w0_4 = np.array(df_4['w0']).reshape(len(pres)-1,nwno)

    #defining lat/lon of cloud base
    lat_base = [int(lat_len/8), int(lat_len - lat_len/8)]
    lon_base = [0, lon_len-1]

    num_spot = 0
    num_dspot = 0
    spot_winds = []; dspot_winds = []; lat_spot = []; spot_len = []; lat_dspot = []; dspot_len = []; spot_init_lon = []; dspot_init_lon = []
    for spot in spots:
        lat_c, lon_c, dlat, dlon, spot_type = spot
        avglatspot = lat_c
        ind_spot = np.abs(ylat - avglatspot).argmin()
        vel_band = xspeed[ind_spot] *3600*timestep 
        tick_spot = int(2 * np.pi * rnep * np.cos(avglatspot * np.pi / 180) / lon_len)
        wind_spot = int(np.round(vel_band/ tick_spot))
        if spot_type == 1:
            num_spot+=1
            spot_winds.append(wind_spot)
            lat_spot.append([int(lat_c - dlat/2), int(lat_c + dlat/2)])
            #spot_len.append([int(lon_c - dlon/2), int(lon_c + dlon/2)])
            spot_len.append(dlon)
            spot_init_lon.append(lon_c)
        elif spot_type == 0:
            num_dspot+=1
            dspot_winds.append(wind_spot)
            lat_dspot.append([int(lat_c - dlat/2), int(lat_c + dlat/2)])
            #dspot_len.append([int(lon_c - dlon/2), int(lon_c + dlon/2)])
            dspot_len.append(dlon)
            dspot_init_lon.append(lon_c)
            

    lon_spacing = np.arange(
        0,
        int(len(lon) * (timedur / p_rot)) + 1,
        int(np.round(len(lon) / (timestep * p_rot))),
        dtype=int,
    )

    rot_count = 0
    cols = math.ceil(math.sqrt(len(lon_spacing)))
    rows = math.ceil(len(lon_spacing) / cols)

    lc_timeseries = []

    for i in range(0, len(lon_spacing)):
        index = (lon_spacing[i] - rot_count * lon_len) % lon_len

        move = [(spot_init_lon[j] + i * spot_winds[j]) % lon_len for j in range(len(spot_winds))]
        dmove = [(dspot_init_lon[j] + i * dspot_winds[j]) % lon_len for j in range(len(dspot_winds))]

        #move = [i * w for w in spot_winds]
        #dmove = [i * w for w in dspot_winds]

        start = [(index + m) % lon_len for m in move]
        dstart = [(index + dm) % lon_len for dm in dmove]

        end = [(s + spot_len[j]) % lon_len for j, s in enumerate(start)]
        dend = [(ds + dspot_len[j]) % lon_len for j, ds in enumerate(dstart)]

        test = np.zeros((lon_len, len(lat), len(pres) - 1, nwno))
        opd_array = np.zeros((len(lon), len(lat), len(pres) - 1, nwno))
        g0_array = np.zeros((len(lon), len(lat), len(pres) - 1, nwno))
        w0_array = np.zeros((len(lon), len(lat), len(pres) - 1, nwno))

        # Prepare arrays
        lon_spot, wrap_spot, lon_spot2 = [], [], []
        lon_dspot, wrap_dspot, lon_dspot2 = [], [], []

        for j in range(len(spot_winds)):
            import pdb; pdb.set_trace()
            if end[j] > start[j]:
                lon_spot.append([start[j], end[j]])
                wrap_spot.append(False)
                lon_spot2.append(None)
            else:
                lon_spot.append([start[j], lon_len])
                lon_spot2.append([0, (start[j] + spot_len[j]) % lon_len])
                wrap_spot.append(True)

        for j in range(len(dspot_winds)):
            if dend[j] > dstart[j]:
                lon_dspot.append([dstart[j], dend[j]])
                wrap_dspot.append(False)
                lon_dspot2.append(None)
            else:
                lon_dspot.append([start[j], lon_len])
                lon_dspot2.append([0, (dstart[j] + dspot_len[j]) % lon_len])
                wrap_dspot.append(True)

        if (lon_spacing[i] % lon_len == 0) and (i != 0):
            rot_count += 1

        for press in range(len(pres) - 1):
            for wave in range(nwno):
                test[lon_base[0] : lon_base[1], lat_base[0] : lat_base[1], press, wave] = 1
                opd_array[
                    lon_base[0] : lon_base[1], lat_base[0] : lat_base[1], press, wave
                ] = opd_base[press, wave]
                g0_array[
                    lon_base[0] : lon_base[1], lat_base[0] : lat_base[1], press, wave
                ] = g0_base[press, wave]
                w0_array[
                    lon_base[0] : lon_base[1], lat_base[0] : lat_base[1], press, wave
                ] = w0_base[press, wave] 
                for j in range(len(lon_spot)):
                    test[lon_spot[j][0] : lon_spot[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave] += 2
                    opd_array[
                            lon_spot[j][0] : lon_spot[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += opd_4[press, wave]
                    g0_array[
                            lon_spot[j][0] : lon_spot[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += g0_4[press, wave]
                    w0_array[
                            lon_spot[j][0] : lon_spot[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += w0_4[press, wave]

                    if wrap_spot[j]:
                        test[lon_spot2[j][0] : lon_spot2[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave] += 2
                        opd_array[
                            lon_spot2[j][0] : lon_spot2[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += opd_4[press, wave]
                        g0_array[
                            lon_spot2[j][0] : lon_spot2[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += g0_4[press, wave]
                        w0_array[
                            lon_spot2[j][0] : lon_spot2[j][1], lat_spot[j][0] : lat_spot[j][1], press, wave
                        ] += w0_4[press, wave]

                for j in range(len(lon_dspot)):
                    test[lon_dspot[j][0] : lon_dspot[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave] = 0.0
                    opd_array[
                            lon_dspot[j][0] : lon_dspot[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ] = 0 
                    g0_array[
                            lon_dspot[j][0] : lon_dspot[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ] = 0 
                    w0_array[
                            lon_dspot[j][0] : lon_dspot[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ]  = 0
                    if wrap_dspot[j]:
                        test[lon_dspot2[j][0] : lon_dspot2[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave] = 0.0
                        opd_array[
                            lon_dspot2[j][0] : lon_dspot2[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ] = 0 
                        g0_array[
                            lon_dspot2[j][0] : lon_dspot2[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ] = 0 
                        w0_array[
                            lon_dspot2[j][0] : lon_dspot2[j][1], lat_dspot[j][0] : lat_dspot[j][1], press, wave
                        ] = 0 


        ds = xr.Dataset(
            data_vars=dict(
                temperature=(
                    ["lon", "lat", "pressure"],
                    temp,
                    {"units": "Kelvin"},
                )  # , required
                # kzz = (["x", "y","z"], gcm_out['kzz'])#could add other data components if wanted
            ),
            coords=dict(
                lon=(["lon"], lon, {"units": "degrees"}),  # required
                lat=(["lat"], lat, {"units": "degrees"}),  # required
                pressure=(["pressure"], pres, {"units": "bar"}),  # required*
            ),
            attrs=dict(description="coords with vectors"),
        )

        ds_chem = xr.Dataset(
            data_vars=dict(
                CH4=(["lon", "lat", "pressure"], chem_ch4, {"units": "v/v"}),
                H2S=(["lon", "lat", "pressure"], chem_h2s, {"units": "v/v"}),
                H2=(["lon", "lat", "pressure"], chem_h2, {"units": "v/v"}),
                He=(["lon", "lat", "pressure"], chem_he, {"units": "v/v"}),
            ),
            coords=dict(
                lon=(["lon"], lon, {"units": "degrees"}),  # required
                lat=(["lat"], lat, {"units": "degrees"}),  # required
                pressure=(["pressure"], pres, {"units": "bar"}),  # required*
            ),
            attrs=dict(description="coords with vectors"),
        )

        ds_withchem = ds.update(ds_chem)
        ds_cld = xr.Dataset(
            data_vars=dict(
                opd=(
                    ["lon", "lat", "pressure", "wno"],
                    opd_array,
                    {"units": "depth per layer"},
                ),
                g0=(["lon", "lat", "pressure", "wno"], g0_array, {"units": "none"}),
                w0=(["lon", "lat", "pressure", "wno"], w0_array, {"units": "none"}),
            ),
            coords=dict(
                lon=(["lon"], lon, {"units": "degrees"}),  # required
                lat=(["lat"], lat, {"units": "degrees"}),  # required
                pressure=(["pressure"], pres[:-1], {"units": "bar"}),  # required
                wno=(
                    ["wno"],
                    cloud_wnogrid,
                    {"units": "cm^(-1)"},
                ),  # required for clouds NEEDS CHEWCKING FOR UNITS AND WHATENOT. what do the 24 wv points map to
            ),
            attrs=dict(description="coords with vectors"),
        )
        neptune.atmosphere_3d(ds_withchem, regrid=True, plot=False, verbose=True)
        neptune.clouds_3d(ds_cld)
        out3d = neptune.spectrum(
            opacity, calculation="reflected", dimension="3d", full_output=True
        )
        wno, alb = (
            out3d["full_output"]["wavenumber"],
            out3d["albedo"],
        )  # wno in cm^-1
        lc_timeseries.append(alb)
        wavelength = wno
        print("*********")
        print("***ran one model yippe horray***")
        print("*********")
    return lc_timeseries, wavelength
    

spot1 = 30, 32, 8, 64, 1 #where spots are a tuple of lat_c, lon_c, dlat, dlon, spot_type
spot2 = 16, 80, 12, 16, 1
spots = [spot1, spot2]
lc_timeseries, wavelength = lc_forMCMC(spots)

lc_timeseries = np.array(lc_timeseries)
np.save("lc_timeseries.npy", lc_timeseries)

h_band = [1.3,2.]
j_band = [1.1, 1.4]  # micrometers
k_band = [1.81,2.59]
v_band = [0.4,0.85] 
kep_band = [0.43, 0.89]
bands = [h_band, j_band, k_band, v_band, kep_band]
phot = []


for band in bands:
    band_index = np.where(((wavelength < band[1]) & (wavelength > band[0])))[0]
    band_sum = np.sum(lc_timeseries[0][band_index[0] : band_index[-1]])
    band_width = wavelength[band_index[0]] - wavelength[band_index[-1]]
    del_lam =abs(np.mean(np.diff(wavelength[band_index[0]:band_index[-1]+1])))
    phot_temp = [] 
    for i in range(len(lc_timeseries)):
        band_sum = np.sum(lc_timeseries[i][band_index[0] : band_index[-1]])
        phot_temp.append(band_sum/len(band_index))
        #phot_temp.append(band_sum * del_lam / band_width)
    phot.append(phot_temp)

initial_V = phot[3][0]; percent_changeV = (phot[3] - initial_V) / initial_V * 100

initial_H = phot[0][0]; percent_changeH = (phot[0] - initial_H) / initial_H * 100

initial_J = phot[1][0]; percent_changeJ = (phot[1] - initial_J) / initial_J * 100

initial_K = phot[2][0]; percent_changeK = (phot[2] - initial_K) / initial_K * 100



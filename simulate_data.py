import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.io import fits

#to simulate a spectrum from an instrument with a given resolution 
wno1d, fpfs1d,fs  = np.load("./ToShare/RefFiles/spec_thrutest.npy")
#fp = fpfs_base*fs #erg/cm2/s/cm
#convert to photon flux - *delta_lam * lambda / hc
h = 6.6261e-27 #erg s 
c = 3e10 #cm s-1 

kepler_benchmark = np.load('longdurhighcadence_darkonly.npy')
alb, wno = kepler_benchmark
sm = 4.514565e12 #m semi major axis of planet
rp = 2.47639e7 #m radius of planet
d_to_nep_close = 29*1.496e+11 #m 
d_to_nep_far = 31*1.496e+11 #m 
d_avg = (d_to_nep_far+d_to_nep_close)/2 #distance from observer to planet
fpfs = alb*(rp/sm)**2
rs = 6.96e8 #m radius of star
fs_obs = fs*(rs/(d_avg))**2
fp = fpfs*fs_obs
wavelength = 1e4/wno



import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def degrade_to_resolution_nonuniform(wavelength_nm, albedo, R=40):
    """
    Degrade a high-resolution spectrum to constant resolving power R
    for a non-uniform wavelength grid.

    Parameters
    ----------
    wavelength_nm : array
        Wavelength array in nm (non-uniform allowed)
    albedo : array
        Model albedo spectrum
    R : float
        Resolving power (lambda/delta_lambda)

    Returns
    -------
    wave_out : array
        Rebinned wavelength array
    albedo_out : array
        Spectrum degraded to resolution R
    """

    wavelength_nm = np.asarray(wavelength_nm)
    albedo = np.asarray(albedo)

    # Sort if needed
    sort_idx = np.argsort(wavelength_nm)
    wavelength_nm = wavelength_nm[sort_idx]
    albedo = albedo[sort_idx]

    # --- Step 1: Move to log-lambda space ---
    log_wave = np.log(wavelength_nm)

    # Create uniform log grid
    dlog = np.min(np.diff(log_wave)) / 2
    log_wave_uniform = np.arange(log_wave[0], log_wave[-1], dlog)

    interp_to_log = interp1d(log_wave, albedo,
                             kind='linear',
                             bounds_error=False,
                             fill_value="extrapolate")

    albedo_log = interp_to_log(log_wave_uniform)

    # --- Step 2: Gaussian convolution in log space ---
    # For constant R, sigma_log = 1 / (R * 2*sqrt(2ln2))
    fwhm_log = 1.0 / R
    sigma_log = fwhm_log / (2 * np.sqrt(2 * np.log(2)))

    sigma_pix = sigma_log / dlog

    albedo_smooth_log = gaussian_filter1d(albedo_log, sigma_pix)

    # --- Step 3: Convert back to wavelength space ---
    interp_back = interp1d(log_wave_uniform,
                           albedo_smooth_log,
                           kind='linear',
                           bounds_error=False,
                           fill_value="extrapolate")

    albedo_smooth = interp_back(log_wave)

    # --- Step 4: Rebin to resolution element sampling ---
    # Use ~2 pixels per resolution element
    #TO PRODUCE A UNIFORM LINEAR GRID, CONSTANT DELTA LAMBDA
    '''
    delta_lambda = wavelength_nm / R
    step = np.mean(delta_lambda) / 2


    wave_out = np.arange(wavelength_nm[0],
                         wavelength_nm[-1],
                         step)
    '''
    #TO PRODUCE CONSTANT R (NOT CONSTANT DELTA LAMBDA)
    dlog_out = 1/(2*R)   # 2 pixels per resolution element
    log_wave_out = np.arange(log_wave[0], log_wave[-1], dlog_out)
    wave_out = np.exp(log_wave_out)


    interp_final = interp1d(wavelength_nm,
                            albedo_smooth,
                            kind='linear',
                            bounds_error=False,
                            fill_value="extrapolate")

    albedo_out = interp_final(wave_out)

    return wave_out, albedo_out

#CONSTANTS
band_ranges = [[0.6, 1.0], [0.850, 0.920]]

center_lam = np.array([0.8,0.885])
delta_lam = np.array([0.2,0.035]) #bandwidth/2
effarea = np.pi*(6.05/2)**2 #in cm^2
throughput = 0.6*0.42
eff_thru =  24#effarea*throughput# = 24 cm^2 effarea*throughput for CUTE 
#eff_thru = effarea*throughput #asteria 
int_time = 60 #s

wave_sim_fp, fp_sim= degrade_to_resolution_nonuniform(
        wavelength[0],
        fp[0],
        R=40
    )
wave_sim_fpfs, fpfs_sim= degrade_to_resolution_nonuniform(
        wavelength[0],
        fpfs[0],
        R=20
    )


noise_level = []; signal = []; signal_flux_per_bandpass = []; noise_per_bandpass = []
for i in range(len(center_lam)):
    signal_flux = []
    signal_noise = [] 
    mask_model = (wave_sim_fp[::-1] > (center_lam[i] - delta_lam[i])) & (wave_sim_fp[::-1] < (center_lam[i] + delta_lam[i]))
    Phi_model_v_test = (1.0 / (h * c)) *  center_lam[i]*1e-4*(np.trapz(fp_sim[mask_model],1e-4*wave_sim_fp[mask_model])) 
    signal_individual = Phi_model_v_test*eff_thru*int_time#*np.pi*(rp/sm)**2
    #snr_signal = signal/np.sqrt(signal)
    #precision = 100/snr_signal
    signal.append(signal_individual)
    #noise_level.append(precision)
#pdb.set_trace()
#yerr = 1e6/snr? 
#SQRT(SIGNAL) WILL CREATE THE ERRORBARS ON THE ALBEDO MEASUREMENTS IN THE STRONG METHANE AND TESS BANDPASSES
signal = np.array(signal)
err = np.sqrt(signal)
#pdb.set_trace()
plt.plot(wavelength[0], fpfs[0])
plt.errorbar(wave_sim_fpfs, fpfs_sim, yerr=fpfs_sim*np.sqrt(fpfs_sim), fmt = 'o')
#plt.xlim(0.85,0.92)
#plt.show()
plt.close()
#pdb.set_trace()

#TO CREATE SERIES OF BINNED DOWN DATA 
wave_time = []; fp_time = []
for k in range(len(wavelength)):
    wave_temp, fp_temp= degrade_to_resolution_nonuniform(
        wavelength[k],
        fp[k],
        R=55
    )
    wave_time.append(wave_temp)
    fp_time.append(fp_temp)
wave_time = np.array(wave_time)
fp_time = np.array(fp_time)

#AND NOW. create change in albedo plots. 
best_bands = [[0.6,1.], [0.850,0.920]]
str_bands = ["TESS", "galileo strong methane (889nm)"]

del_per = []
#pdb.set_trace()
for i in range(len(center_lam)):
    phot = []; del_per_temp = []
    for j in range(len(wave_time)):
        mask_model = (wave_time[j] > (center_lam[i] - delta_lam[i])) & (wave_time[j]< (center_lam[i] + delta_lam[i]))
        Phi_model_v_test = (1.0 / (h * c)) *  center_lam[i]*1e-4*(np.trapz(fp_time[j][mask_model],1e-4*wave_time[j][mask_model])) 
        signal_individual = Phi_model_v_test*eff_thru*int_time
        phot_temp = signal_individual
        phot.append(phot_temp) #gives the timeseries data for an individual bandpass
    phot = np.array(phot)
    mean = np.mean(phot)
    #pdb.set_trace()
    for k in range(len(phot)):
        err = np.sqrt(phot[k])
        percent_change = ((phot[k] - mean) / mean) * 100
        percent_change_err = err*100/mean
        del_per_temp.append([percent_change, percent_change_err])
    del_per.append(del_per_temp)
#TESS first then strong methane 
del_per = np.array(del_per)
#pdb.set_trace()
order_sorted = np.arange(0, len(wave_time))
precision_cute = [0.24607, 0.474]  # TESS then strong methane
#precision_asteria = [100 / 632.7105, 100 / 86.298]import pdb; pdb.set_trace()
rot_images_x = [order_sorted[3], order_sorted[19], order_sorted[72], order_sorted[105], order_sorted[123], order_sorted[171], order_sorted[210]]
rot_images_y = [del_per[0].T[0][3], del_per[0].T[0][19], del_per[0].T[0][72],  del_per[0].T[0][123], del_per[0].T[0][171], del_per[0].T[0][210]]
#del_per[0].T[0][105]
plt.figure(figsize=(12,5))
plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
colors = ['darkblue', 'darkred']
for i in range(len(str_bands)):
    del_per_x = order_sorted[::3]#[x for k, x in enumerate(del_per[0][0]) if (k + 1) % 3 != 0]

    del_per_y = del_per[i].T[0][::3]#[x for k, x in enumerate(del_per[i][1]) if (k + 1) % 3 != 0]
    #del_per_y[35] = (del_per_y[69] + del_per_y[12])/2
    del_per_yerr = del_per[i].T[1][::3]
    pdb.set_trace()
    if i !=1000:
        plt.errorbar(
            del_per_x,
            del_per_y,
            #del_per[i][1],
            #yerr=del_per_yerr,
            fmt = '.-',
            label=str_bands[i],
            color='k',
        )  
        plt.fill_between(
            del_per_x,
            np.array(del_per_y)-del_per_yerr, 
            y2 = np.array(del_per_y)+del_per_yerr,
            color = colors[i],
            alpha = 0.5
        )

#plt.plot(rot_images_x, rot_images_y, 'ro', markersize = 8, mec='k', zorder = 5)
plt.ylabel("% change in albedo")
plt.xlabel("Time [hours]")
#plt.title("Change in Albedo over One Rotation")
# plt.ylim(-7,4)
# plt.ylim(-2,2)

#plt.savefig("TESS_Long.pdf")
#plt.show()
plt.close()
#pdb.set_trace()
'''
R = np.array([15, 20, 25, 30, 35, 40, 45,
    50, 55, 60, 65, 70, 75,
    100, 150, 200
])

values = np.array([
    3.8010829,
    3.46078241,
    4.21326,
    3.53,
    4.09544,
    4.99,
    5.35,
    4.7912,
    5.39,
    5.647,
    5.9785,
    5.9907,
    5.849,
    6.23,
    6.46116890,
    6.50
])
plt.figure(figsize=(10,5))
plt.axhline(y=5.0, color="k", linestyle="--", alpha=0.5)
plt.plot(R, values, marker='o', color = 'k')
R_fill = np.concatenate([np.array([0]),R,np.array([300])])
plt.fill_between(
    R_fill,
    0.0, 
    y2 = 5.0,
    color = 'darkred',
    alpha = 0.5,
    zorder = 5
    )



plt.xlabel("R")
plt.ylabel("SNR")
plt.ylim(3.0,6.75)
plt.xlim(5,210)
plt.xscale("linear")   # change to "log" if desired
plt.grid(True)
#plt.show()
plt.savefig("resolutionplot.pdf")
plt.close()
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# PARAMETERS
# -----------------------------
months = 12
days_per_month = 30              # simplified uniform months
total_days = months * days_per_month

solar_window_months = 6          # 6 months observable
lunar_cycle = 28                 # days
lunar_visible = 25               # visible days per cycle

# -----------------------------
# CREATE DAY ARRAY
# -----------------------------
days = np.arange(total_days)

# -----------------------------
# SOLAR EXCLUSION
# First 6 months blocked, last 6 months observable
# -----------------------------
solar_blocked = np.zeros(total_days, dtype=bool)
solar_blocked[:6 * days_per_month] = True

# -----------------------------
# LUNAR EXCLUSION (during solar-visible window)
# -----------------------------
lunar_blocked = np.zeros(total_days, dtype=bool)

for i in range(total_days):
    if not solar_blocked[i]:
        day_in_cycle = i % lunar_cycle
        if day_in_cycle >= lunar_visible:
            lunar_blocked[i] = True

# -----------------------------
# FINAL VISIBILITY FLAG
# -----------------------------
visible = (~solar_blocked) & (~lunar_blocked)

# -----------------------------
# PLOT
# -----------------------------
fig, ax = plt.subplots(figsize=(10,3))

# Sun blocked
ax.plot(days[solar_blocked], np.ones(np.sum(solar_blocked))*1,
        '|', color='gold', markersize=12, label='Sun Blocked')

# Moon blocked
ax.plot(days[lunar_blocked & ~solar_blocked],
        np.ones(np.sum(lunar_blocked & ~solar_blocked))*2,
        '|', color='gray', markersize=12, label='Moon Blocked')

# Visible
ax.plot(days[visible], np.ones(np.sum(visible))*3,
        '|', color='blue', markersize=12, label='Visible')

# Formatting
ax.set_yticks([1,2,3])
ax.set_yticklabels(['Sun Blocked', 'Moon Blocked', 'Visible'])

# Month tick marks
month_ticks = np.arange(0, total_days+1, days_per_month)
month_labels = [f"Month {i}" for i in range(1, months+1)]
ax.set_xticks(month_ticks[:-1])
ax.set_xticklabels(month_labels, rotation=45)

ax.set_xlim(0, total_days)
ax.set_title("1-Year Visibility Calendar (Conceptual)")
ax.grid(axis='x', alpha=0.3)
#ax.legend(loc='upper right')

plt.tight_layout()
#plt.show()
plt.close()
#plt.savefig("VisibilityCalendarMockup.pdf")


wv = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
im4 = np.array([9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 9.9999977e-07, 4.6656598e-05, 4.9251183e-05, 3.0203702e-05, 6.6141161e-05, 0.00021260191, 0.00032470401, 0.00018452496, 0.00032499802, 0.00057240902, 0.0010081698, 0.0017756603, 0.0031274103, 0.005927172, 0.0083943801])
im1 = np.array([0.00050702301, 0.00015172719, 5.2704409e-06, 2.1634794e-06, 0.00020199531, 0.0011379919, 0.0017910724, 0.00038337187, 0.0017504509, 0.002432242, 0.0031117215, 0.0056596263, 0.0037117724, 0.00044349086, 0.00034130135, 0.001089377, 0.00094064744, 0.0010821783, 0.0011891199, 0.0010549536, 0.00089978962, 0.0010253576, 0.0010441705, 0.0010347319])
im2 = np.array([0.0011754304, 0.003398696, 0.0012706247, 0.0017925976, 9.7051605e-05, 0.0002241929, 0.00021853589, 0.00030888562, 0.0066479365, 0.010653216, 0.026788043, 0.035500087, 0.3499327, 0.037856936, 0.045398597, 0.054237988, 0.002701669, 0.0013209067, 0.00098555454, 0.00080591417, 0.00085984101, 0.00078390876, 0.00077089731, 0.00075711781])
im3 =  np.array([0.00071043964, 0.0044625239, 0.0076165493, 0.0031423781, 0.0013061292, 0.017672312, 0.034849342, 0.040612739, 0.0023888722, 0.0018183715, 0.0017525674, 0.0025047383, 0.012611357, 0.2766216, 0.3261088, 0.3569286, 0.03764772, 0.06004164, 0.9918588, 0.013779141, 0.0025911264, 0.0014113548, 0.0011158295, 0.0010275072])
plt.close()
plt.plot(wv,im1, label = 'im1')
plt.plot(wv,im2, label = 'im2')
plt.plot(wv,im3, label = 'im3')
plt.plot(wv,im4, label = 'im4')
plt.legend()
plt.show()
'''
wno1d, fpfs1d,fs  = np.load("./ToShare/RefFiles/spec_thrutest.npy")
solar_r = 6.957e+8 #m
sm = 4.514565e12 #m semi major axis of planet

kepler_benchmark = np.load('brightspotsforppm.npy')
alb_bright, wno = kepler_benchmark
kepler_benchmark = np.load('darkspotsforppm.npy')
alb_dark, wno = kepler_benchmark

#need alb of haze only case, bright only case, dark only case 

wave_bright, alb_time= degrade_to_resolution_nonuniform(
    wavelength[0],
    alb[19],
    R=55
)

wave_bright, alb_bright= degrade_to_resolution_nonuniform(
    wavelength[0],
    alb_bright[1],
    R=55
)

wave_dark, alb_dark= degrade_to_resolution_nonuniform(
    wavelength[0],
    alb_dark[1],
    R=55
)

wave_time, fs_sim= degrade_to_resolution_nonuniform(
        wavelength[0],
        fs_obs,
        R=55
    )

reflected_flux_haze= alb_time*fs_sim*solar_r**2/sm**2
reflected_flux_bright = alb_bright*fs_sim*solar_r**2/sm**2
reflected_flux_dark = alb_dark*fs_sim*solar_r**2/sm**2

plt.plot(wave_time, 1e6*(reflected_flux_bright-reflected_flux_haze)/reflected_flux_haze)
plt.plot(wave_time, 1e6*(reflected_flux_dark-reflected_flux_haze)/reflected_flux_haze)
plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
#plt.xlim(0.6,1.0)
#plt.plot(wave_time, reflected_flux_cloud/alb_time[19])
#plt.ylim(-50000,50000)
plt.show()
pdb.set_trace()
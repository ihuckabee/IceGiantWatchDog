import numpy as np
import matplotlib.pyplot as plt
import math 
import pdb

nep_ptprofile = "/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/nepforpicaso_salr_abridged.txt"
f = open(nep_ptprofile, "r")
lines = f.readlines()[1:]
pressure = []
for i in range(len(lines)):
    pressure.append(float(lines[i].split()[0]))

lon = np.linspace(-180, 180, 128)
lat = np.linspace(-90, 90, 128)
pres = np.array(pressure)
nwno = 24
p_rot = 16.11  # hours
timestep = 2  # hours
timedur = 3*p_rot  # hours
rnep = 0.3463884070945 * 7.1492e07  # meters
# defining windspeed per latitude
v_spot = 300 * 3600  # m/hour
v_dspot = 150 * 3600 
lat_base = [16, 112]
lon_base = [0, 128]

lat_spot = [60, 70]
lat_dspot = [30, 50]

avglatspot = (lat[lat_spot[1]] + lat[lat_spot[0]]) / 2
avglatdspot = (lat[lat_dspot[1]] + lat[lat_dspot[0]]) / 2

tick_spot = int(2 * np.pi * rnep * np.cos(avglatspot * np.pi / 180) / 128)
tick_dspot = int(2 * np.pi * rnep * np.cos(avglatdspot * np.pi / 180) / 128)

spot_width = 64  # [32,96]
dspot_width = 20  # [50,70]

spot_wind = int(np.round(v_spot/tick_spot))
dspot_wind = int(np.round(v_dspot/tick_spot))
lon_spacing = np.arange(0, int(len(lon)*(timedur/p_rot)) + 1, int(np.round(len(lon) / (timestep * p_rot))), dtype=int)

cols = math.ceil(math.sqrt(len(lon_spacing)))
rows = math.ceil(len(lon_spacing) / cols)
pdb.set_trace()
# Create figure with subplots
fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))
axes = axes.flatten()  # Flatten to 1D array for easy looping
#import pdb; pdb.set_trace()
rot_count = 0
lon_len = len(lon)
'''
for i in range(len(lon_spacing)):
    lon_spacing[i] = lon_spacing[i] - rot_count*len(lon)
    wind_move = iter*spot_wind
    test = np.zeros((len(lon), len(lat), len(pres) - 1, nwno))    
    if (lon_spacing[i] + spot_width + wind_move) < len(lon):
        lon_spot = [lon_spacing[i]+wind_move, lon_spacing[i] + spot_width + wind_move]
        wrap_spot = False
    if (lon_spacing[i] + wind_move) >= len(lon):
        lon_spot = [np.abs(len(lon) - 1 - lon_spacing[i] - wind_move), lon_spacing[i] + spot_width + wind_move  - len(lon)]
        wrap_spot = False
        if lon_spot[0] > len(lon):
           lon_spot = [lon_spot[0] - len(lon), lon_spot[1] - len(lon)]
        if (lon_spot[1] > len(lon)) and (lon_spot[0] < len(lon)):
           lon_spot = [lon_spot[0], len(lon) - 1]
           lon_spot2 = [0, lon_spot[1]-len(lon)]
           wrap_spot = True
    elif ((lon_spacing[i] + wind_move + spot_width) >= len(lon)) and ((lon_spacing[i] + spot_width) < len(lon)):
        lon_spot = [lon_spacing[i] + wind_move, len(lon) - 1]
        lon_spot2 = [0, lon_spacing[i] + spot_width + wind_move  - len(lon)]
        wrap_spot = True
    if (lon_spacing[i] + dspot_width + dspot_wind) < len(lon):
        lon_dspot = [lon_spacing[i] + dspot_wind, lon_spacing[i] + dspot_width + dspot_wind]
        wrap_dspot = False
    elif ((lon_spacing[i] + dspot_wind + dspot_width) >= len(lon)) and ((lon_spacing[i] + dspot_width) < len(lon)):
        lon_dspot = [lon_spacing[i] + dspot_wind, len(lon) - 1]
        lon_dspot2 = [0, lon_spacing[i] + dspot_width + dspot_wind - len(lon)]
        wrap_dspot = True
    if (lon_spacing[i] % 128 == 0) and (i != 0):
       rot_count +=1
       #pdb.set_trace()
    
    print(lon_spot)
    print(iter)
    print(wrap_spot)
    pdb.set_trace()
    #if iter > 6:
        #iter = 0
    #else:
    iter +=1
'''
for i in range(len(lon_spacing)):
    index = (lon_spacing[i] - rot_count * lon_len) % lon_len
    move = i * spot_wind
    start = (index + move) % lon_len
    end = (start + spot_width) % lon_len

    dmove = i*dspot_wind
    dstart = (index + dmove) % lon_len
    dend = (dstart + dspot_width) % lon_len
    # Allocate test array
    test = np.zeros((lon_len, len(lat), len(pres) - 1, nwno))

    # Determine if wrapping is needed
    if end > start:
        lon_spot = [start, end]
        wrap_spot = False
    else:
        # Wrapped around: split into two parts
        first_len = lon_len - start
        second_len = spot_width - first_len
        lon_spot = [start, lon_len]
        lon_spot2 = [0, second_len]
        wrap_spot = True

    if dend > dstart:
        lon_dspot = [dstart, dend]
        wrap_dspot = False
    else:
        # Wrapped around: split into two parts
        first_len = lon_len - dstart
        second_len = dspot_width - first_len
        lon_dspot = [dstart, lon_len]
        lon_dspot2 = [0, second_len]
        wrap_dspot = True
    # Update rot_count when lon_spacing loops
    if (lon_spacing[i] % lon_len == 0) and (i != 0):
        rot_count += 1

    # Optional debug print
    print(f"Iter {i} | lon_spot: {lon_spot}" + (f", lon_spot2: {lon_spot2}" if wrap_spot else ""))
    

    for press in range(len(pres) - 1):
        for wave in range(nwno):
            test[lon_base[0] : lon_base[1], lat_base[0] : lat_base[1], press, wave] = 1
 # test[press,wav]
            if wrap_spot == False:
                test[lon_spot[0] : lon_spot[1], lat_spot[0] : lat_spot[1], press, wave] += 2
            elif wrap_spot == True:
                test[
                    lon_spot[0] : lon_spot[1], lat_spot[0] : lat_spot[1], press, wave
                ] += 2  # test[press,wav]
                test[
                    lon_spot2[0] : lon_spot2[1], lat_spot[0] : lat_spot[1], press, wave
                ] += 2
            if wrap_dspot == False:
                test[
                    lon_dspot[0] : lon_dspot[1],
                    lat_dspot[0] : lat_dspot[1],
                    press,
                    wave,
                ] = 0.0
            elif wrap_dspot == True:
                test[
                    lon_dspot[0] : lon_dspot[1],
                    lat_dspot[0] : lat_dspot[1],
                    press,
                    wave,
                ] = 0.0  # test[press,wav]
                test[
                    lon_dspot2[0] : lon_dspot2[1],
                    lat_dspot[0] : lat_dspot[1],
                    press,
                    wave,
                ] = 0.0  # test[press,wav]

    im = axes[i].imshow(
        test[:, :, 0, 0].T, extent=[-180, 180, -90, 90], origin="lower", cmap="bone"
    )
    #axes[i].colorbar()
    axes[i].set_title(f"Timestep {i}")
    axes[i].axis("off")


for j in range(len(lon_spacing)):
    axes[j].axis("off")
plt.tight_layout()
plt.show()


x = [84.74830200809845, 5.702539541166175, -302.9630437924129, -0.6110706417716187, 21.88028088015409, 57.00209686885279, 100.23475493766904, 135.2164847512512, 173.70833173585942, 213.37306865402252, 226.26247135606366, 226.88356921879108, 200.1515879873893, 173.16391262161693, 130.78622259840688, 100.37248176667839, 64.31135110904012, 30.584498506239015, -24.794663159105426, -63.35817528172066, -115.88665701300883, -168.7454471993087, -209.97000681261187, -249.703459686148, -268.7985324867213, -290.42577349687537, -324.6896722573339, -350.36348675964456, -376.7138437582024, -378.480993992515, -385.6395450001032, -375.2470677751661, -387.6672555201326, -365.9923852104389, -356.5974690114635, -330.6384685484418, -349.6350386490384, -296.16134388354266, -320.98078017677403, -265.86851953980954, -219.49055819179713, -187.07491218808707, -134.07102221619277, -86.29367370833174, -55.43165711623024, -28.96185822099403, 9.811340786900928, 46.117548521427295, 100.47570315886952, 139.88769512529603, 175.63709717851697, 221.9543053642683, 202.63450484698166, 230.59010195324367, 226.18756211715925, 215.9680544537082, 193.2466666863279, 170.8287498119896, 146.13848771812894, 117.34856684469901, 45.68519834727806]
y = [85.73089045952634, 90.09580754182505, -90.05805191802098, -90.18797360666262, -88.93039376331984, -87.37117824284935, -84.93721784728332, -82.48475366749034, -79.8885640607079, -75.84068130434821, -73.32667898631756, -68.99693461972716, -64.4541638493705, -61.720543771941294, -58.03256032450861, -55.74603740417893, -53.553360144661674, -51.46370280968492, -47.39443696171282, -44.912389309959714, -41.734509040060715, -37.740585068674065, -34.69395322709431, -31.303441752481973, -29.255929573646874, -26.365796557016438, -21.996391142616915, -17.032818066430607, -9.145830567658749, -4.585233497165817, 3.7909680863434687, 10.982673422793127, -1.0032198529185195, -13.452865102959151, 14.276657528734077, 20.12563879731833, 17.217792394566743, 25.074257823627846, 22.63755220606933, 29.383861152772113, 33.4666907968099, 36.83219887420023, 40.30223710108322, 43.84222556032992, 45.79815035865491, 47.393703020614566, 50.224732607458265, 52.46521483812484, 55.92755374060218, 58.841906813474765, 62.12347037775305, 67.92221196673711, 64.91511414469946, 70.69777960078153, 72.80295202714765, 74.84768369836056, 77.39196814381984, 79.95918119147339, 82.03796349536938, 84.00004916464408, 88.4174889067861]
plt.plot(x,y,'.')
plt.show()
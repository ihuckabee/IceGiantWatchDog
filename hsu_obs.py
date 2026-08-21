import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import numpy as np
from datetime import datetime, timezone, timedelta
import matplotlib.dates as mdates
def utc_to_jd(date, time):
    t = Time(f"{date}T{time}", scale='utc')
    return t.jd

def convert_dataset_to_jd(raw):
    jd_list = []
    for date, time, filt in raw:
        jd = utc_to_jd(date, time)
        jd_list.append((date, time, filt, jd))
    return jd_list


# --- helper: datetime → JD ---
def datetime_to_jd(dt):
    # Unix epoch to JD conversion
    return dt.timestamp() / 86400.0 + 2440587.5
# JD → datetime
def jd_to_datetime(jd):
    return datetime(1970,1,1) + timedelta(days=(jd - 2440587.5))

def build_data_dict(hst_dark):
    data = {
        "HST": {},
        "VLT": {}
        }
    
    def process(raw, obs_name):
        for date_str, time_str, filt in raw:
            dt = datetime.fromisoformat(f"{date_str}T{time_str}").replace(tzinfo=timezone.utc)
            jd = datetime_to_jd(dt)

            if filt not in data[obs_name]:
                data[obs_name][filt] = []

            data[obs_name][filt].append(jd)

    process(hst_dark, "HST")

    # convert lists → numpy arrays
    for obs in data:
        for filt in data[obs]:
            data[obs][filt] = np.array(data[obs][filt])

    return data

hst_dark = [
    ("1994-06-28","15:01:00","F467M"),
    ("1994-06-28","23:00:00","F467M"),
    ("1994-10-10","14:15:00","F410M"),
    ("1994-10-10","14:19:00","F467M"),
    ("1994-10-10","20:31:00","F410M"),
    ("1994-10-10","20:35:00","F467M"),
    ("1994-10-11","04:34:00","F410M"),
    ("1994-10-11","04:38:00","F467M"),
    ("1994-10-18","15:14:00","F410M"),
    ("1994-10-18","15:18:00","F467M"),
    ("1994-10-18","23:25:00","F410M"),
    ("1994-10-18","23:29:00","F467M"),
    ("1994-10-19","04:07:00","F410M"),
    ("1994-10-19","04:11:00","F467M"),
    ("1994-11-01","23:30:00","F410M"),
    ("1994-11-01","23:34:00","F467M"),
    ("1994-11-02","06:01:00","F410M"),
    ("1994-11-02","06:05:00","F467M"),
    ("1994-11-02","12:28:00","F410M"),
    ("1994-11-02","12:32:00","F467M"),

    ("1995-09-01","01:30:00","F467M"),
    ("1995-09-01","07:56:00","F467M"),
    ("1995-09-01","14:24:00","F467M"),
    ("1995-09-01","22:25:00","F467M"),
    ("1995-09-02","04:51:00","F467M"),
    ("1995-09-02","11:17:00","F467M"),

    ("1995-11-22","04:01:00","F410M"),
    ("1995-11-22","04:05:00","F467M"),
    ("1995-11-22","10:27:00","F410M"),
    ("1995-11-22","10:31:00","F467M"),

    ("1996-03-08","01:29:00","F410M"),
    ("1996-03-08","01:33:00","F467M"),
    ("1996-03-08","09:31:00","F410M"),
    ("1996-03-08","09:35:00","F467M"),

    ("1996-08-13","08:33:00","F467M"),
    ("1996-08-13","10:18:00","F467M"),
    ("1996-08-13","11:38:00","F467M"),
    ("1996-08-13","12:15:00","F410M"),
    ("1996-08-13","15:00:00","F467M"),
    ("1996-08-13","16:45:00","F467M"),
    ("1996-08-13","18:04:00","F467M"),
    ("1996-08-13","18:42:00","F410M"),

    ("1997-07-03","09:14:00","F439W"),
    ("1997-07-03","09:17:00","F439W"),
    ("1997-07-03","15:41:00","F439W"),
    ("1997-07-03","15:44:00","F439W"),

    ("1997-07-05","04:48:00","F410M"),
    ("1997-07-05","04:52:00","F467M"),
    ("1997-07-05","11:15:00","F410M"),
    ("1997-07-05","11:19:00","F467M"),

    ("1998-08-11","02:07:00","F467M"),
    ("1998-08-11","05:20:00","F467M"),
    ("1998-08-11","06:56:00","F467M"),
    ("1998-08-11","10:11:00","F467M"),

    ("2000-08-22","03:36:00","F467M"),
    ("2001-06-06","02:10:00","F467M"),
    ("2001-06-25","21:22:00","F467M"),

    ("2002-08-09","12:22:00","F467M"),
    ("2002-08-09","13:34:00","F467M"),
    ("2002-08-09","14:17:00","F467M"),

    ("2004-11-06","09:02:00","F475W"),
    ("2005-04-29","19:35:00","F475W"),
    ("2005-04-29","19:59:00","F435W"),

    ("2006-07-10","10:38:00","F475W"),
    ("2007-08-11","06:00:00","F467M"),

    ("2009-08-19","13:16:00","F475W"),
    ("2010-08-28","13:07:00","F467M"),

    ("2011-06-25","19:37:00","F467M"),
    ("2011-06-25","22:59:00","F467M"),

    ("2015-09-18","07:17:00","F467M"),
    ("2015-09-18","12:03:00","F467M"),

    ("2016-05-15","00:52:00","F467M"),
    ("2016-05-15","01:07:00","F410M"),
    ("2016-05-16","00:42:00","F467M"),
]

add_opal_dates = [
    "2014-11-15", #k2
    "2015-01-18", #k2
    "2016-10-04",
    "2017-10-06", "2017-10-07",
    "2018-09-09", "2018-09-10",
    "2018-11-05", "2018-11-06",
    "2019-09-28", "2019-09-29",
    "2020-08-19", "2020-08-20",
    "2021-09-06", "2021-09-07", "2021-09-08",
    "2022-09-18", "2022-09-19", "2022-09-20",
    "2023-09-22", "2023-09-23",
    "2024-08-24", "2024-08-25",
    "2025-06-28",
    "2025-08-24", "2025-08-25"
]

add_filters = ["F845M", "F467M", "F547M", "F657M", "F763M", "FQ727N", "FQ619N"]

add_VLT_dates = [
    "2019-10-18",
    "2019-11-13"
]
add_VLT_filters= ["551nm"]

from datetime import datetime, timezone

def add_new_observations(data, dates, filters, obs_name="HST"):
    for date_str in dates:
        # assume midnight UTC
        dt = datetime.fromisoformat(f"{date_str}T00:00:00").replace(tzinfo=timezone.utc)
        jd = datetime_to_jd(dt)

        for filt in filters:
            if filt not in data[obs_name]:
                data[obs_name][filt] = []

            # if already numpy array → convert to list temporarily
            if isinstance(data[obs_name][filt], np.ndarray):
                data[obs_name][filt] = data[obs_name][filt].tolist()

            data[obs_name][filt].append(jd)

    # convert back to sorted numpy arrays
    for filt in data[obs_name]:
        data[obs_name][filt] = np.sort(np.array(data[obs_name][filt]))

    return data

data = build_data_dict(hst_dark)

data = add_new_observations(data, add_opal_dates, add_filters, obs_name="HST")
data = add_new_observations(data, add_VLT_dates, add_VLT_filters, obs_name = 'VLT' )

def plot_observation_times(data):
    plt.figure()

    y_offset = 0
    yticks = []
    ylabels = []

    for obs, filters in data.items():
        for filt, jd in filters.items():
            jd_sorted = np.sort(jd)
            
            y = np.ones_like(jd_sorted) * y_offset
            plt.scatter(jd_sorted, y, label=f"{obs} - {filt}", zorder = 10)

            yticks.append(y_offset)
            ylabels.append(f"{obs}-{filt}")
            y_offset += 1

    plt.xlabel("Julian Date")
    plt.yticks(yticks, ylabels)
    plt.title("VIS Observation Timeline")
    plt.grid()
    plt.legend()
    plt.show()

def plot_time_gaps(data):
    plt.figure()

    for obs, filters in data.items():
        for filt, jd in filters.items():
            jd_sorted = np.sort(jd)
            dt = np.diff(jd_sorted)  # time gaps in days

            label = f"{obs} - {filt}"
            plt.plot(dt/365.25, marker='o', linestyle='-', label=label)

    plt.xlabel("Observation Index")
    plt.ylabel("Δt (years)")
    plt.title("Time Between Consecutive Observations")
    plt.legend()
    plt.grid()
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def plot_time_gaps_consolidated(data):
    plt.figure(figsize=(10, 6))

    for obs in data:
        # --- collect all times across filters ---
        all_times = []
        for filt in data[obs]:
            all_times.extend(data[obs][filt])

        if len(all_times) < 2:
            continue

        # --- sort ---
        all_times = np.sort(np.array(all_times))
        all_times = np.unique(np.sort(np.array(all_times)))
        # --- compute gaps ---
        dt = np.diff(all_times)

         # --- convert to datetime ---
        times_dt = np.array([jd_to_datetime(jd) for jd in all_times])
        
        # --- plot (aligned with time of second point) ---
        plt.plot(times_dt[1:], dt/365.25, marker='o', linestyle = 'None', label=obs)
         # highlight the "missing" cadence region (~6 months)
    #import pdb; pdb.set_trace()
    #plt.plot(times_dt[1:][-30], dt[-30]/365.25, 'D', color = 'indigo', label = 'K2', markersize = 7)#, xerr = np.array([[datetime(2014, 12, 15, 0, 0) - datetime(2014, 11, 15, 0, 0)],[datetime(2015, 1, 18, 0, 0) - datetime(2014, 12, 15, 0, 0)]]), elinewidth=20)
    plt.axhspan(0.25, 0.75, color='red', alpha=0.2, label=r'$\Delta t$ = 3-9 months')
    plt.axhspan(0.020833, 0.0833, color='orange', alpha=0.2, label=r'$\Delta t$ = 1-4 weeks')


    # optional: also show 1 year reference
    #plt.axhline(365/365, color='gray', linestyle='--', alpha=0.5, label='1 year')

    plt.xlabel("Year")
    plt.ylabel("time since last observation (years)")
    plt.title("Time Gaps Between Observations (VIS)")
    plt.legend()
    #plt.yscale("log")  # optional but very helpful
    plt.grid(True)

    plt.show()


plot_time_gaps_consolidated(data)
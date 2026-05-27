import datetime
from xml.parsers.expat import model
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy.taup import TauPyModel
import os


REFERENCE_COORDINATE_MMS2 = (35.4134,
        138.7738)
REFERENCE_COORDINATE_FUJ = (

    35.31356 ,138.67578
)


def calculate_P_S_travel_time(lon, lat, depth,REFERENCE_COORDINATE):
    """Calculate P and S wave tt with taup"""

    model = TauPyModel(model="iasp91")
    epicentral_distance = obspy.geodetics.locations2degrees(
        lat1=lat,
        long1=lon,
        lat2=REFERENCE_COORDINATE[0],
        long2=REFERENCE_COORDINATE[1],
    )
    try:
        arrivals = model.get_travel_times(
            source_depth_in_km=depth,
            distance_in_degree=epicentral_distance,
        )
        p_time = next(
            (arr.time for arr in arrivals if "p" in arr.phase.name.lower()),
            None,
        )
        s_time = next(
            (arr.time for arr in arrivals if "s" in arr.phase.name.lower()),
            None,
        )
        r_time = (111 * epicentral_distance) / 3.5
        return p_time, s_time, r_time
    except Exception as e:
        print(
            f"Error calculating travel times for ({lon}, {lat}, {depth}): {e}"
        )
        return None, None, None


def subsidiary_info_to_string(info_code):
    """Convert subsidiary information code to human-readable string."""
    if pd.isna(info_code):
        return None

    info_code = int(info_code)

    if info_code == "":
        return None

    mapping = {
        "1": "Tectonic earthquake",
        "2": "Lack of information",
        "3": "Artificial event",
        "4": "Eruption earthquake and others",
        "5": "Low-frequency earthquake",
    }

    return mapping.get(info_code, info_code)


def decode_magnitude(mag_str):
    """Decode JMA magnitude encoding to numeric value.

    Format F2.1 means implied decimal point - divide by 10.
    Special encoding:
    - -0.1: -1, -0.9: -9
    - -1.0: A0, -1.9: A9
    - -2.0: B0, -2.9: B9
    - -3.0: C0, -3.9: C9
    """
    if pd.isna(mag_str) or mag_str == "":
        return None

    mag_str = str(mag_str).strip()

    if mag_str == "":
        return None

    # Check for letter codes (A/B/C format)
    if len(mag_str) >= 2 and mag_str[0] in ["A", "B", "C"]:
        letter = mag_str[0]
        try:
            digit = int(mag_str[1])
            if letter == "A":
                return -1.0 - (digit / 10.0)
            elif letter == "B":
                return -2.0 - (digit / 10.0)
            elif letter == "C":
                return -3.0 - (digit / 10.0)
        except (ValueError, IndexError):
            return None

    # Try to convert to numeric - all values have implied decimal (F2.1 format)
    try:
        val = int(mag_str)
        # All numeric values (positive or negative) are in F2.1 format with implied decimal
        return val / 10.0
    except ValueError:
        return None


def degrees_minutes_to_decimal(degrees, minutes):
    """Convert degrees and minutes to decimal degrees.

    Args:
        degrees: Degrees component
        minutes: Minutes component

    Returns:
        Decimal degrees
    """
    if pd.isna(degrees) or pd.isna(minutes):
        return None

    return degrees + (minutes / 60.0 / 100.0)


def convert_depth(depth):
    return depth / 100


def arrival_time(q_time, p_time, s_time, r_time):
    if pd.isna(p_time) or pd.isna(s_time) or pd.isna(r_time):
        return None, None, None
    arr_p_time = q_time + datetime.timedelta(seconds=p_time)
    arr_s_time = q_time + datetime.timedelta(seconds=s_time)
    arr_r_time = q_time + datetime.timedelta(seconds=r_time)
    return arr_p_time, arr_s_time, arr_r_time


def seconds(second, millisecond):
    if pd.isna(second) or pd.isna(millisecond):
        return None
    return second + (millisecond / 100.0)


# def Rayleigh_time(epicentrale_distance):
#     """Calculate Rayleigh wave travel time using a simple empirical formula."""

#     dist_km = 111.1 * epicentrale_distance
#     travel_time = dist_km / 3.5
#     return travel_time


def read_catalog(
    start=None,
    end=None,
    max_distance=None,
    max_depth=None,
    min_magnitude=None,
    LFE=False,
    all_phases=False,
    station_name=None,
):
    """Load the JMA catalog and optionally filter by time range.

    Args:
        start: Optional start time (datetime-like). If provided, only records
            with index >= start are kept.
        end: Optional end time (datetime-like). If provided, only records
            with index <= end are kept.
    """

    # Read file with FWF
    print("Reading JMA catalog... ", end="")
    df = pd.read_fwf(
        "/gpfs/scratch/borde/fuji/h2008",
        colspecs=[
            (0, 1),  # Record type identifier
            (1, 5),  # Year
            (5, 7),  # Month
            (7, 9),  # Day
            (9, 11),  # Hour
            (11, 13),  # Minute
            (13, 15),  # Second
            (15, 17),  # Millisecond
            (17, 21),  # Standard error (seconds)
            (21, 24),  # Latitude (degrees)
            (24, 28),  # Latitude (minutes)
            (28, 32),  # Standard error (minutes) for latitude
            (32, 36),  # Longitude (degrees)
            (36, 40),  # Longitude (minutes)
            (40, 44),  # Standard error (minutes) for longitude
            (44, 47),  # Depth (kilometers)
            (47, 49),  # Depth decimal
            (49, 52),  # Standard error (kilometers) for depth
            (52, 54),  # Magnitude 1
            (54, 55),  # Magnitude type 1
            (55, 57),  # Magnitude 2
            (57, 58),  # Magnitude type 2
            (58, 59),  # Travel time table
            (59, 60),  # Hypocenter location precision
            (60, 61),  # Subsidiary information
            (61, 62),  # Maximum intensity
            (62, 63),  # Damage class
            (63, 64),  # Tsunami class
            (64, 65),  # District number
            (65, 68),  # Region number
            (68, 92),  # Region name
            (92, 95),  # Number of stations
            (95, 96),  # Hypocenter determination flag
        ],
        names=[
            "record_type",
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "millisecond",
            "second_error",
            "latitude_degrees",
            "latitude_minutes",
            "latitude_error",
            "longitude_degrees",
            "longitude_minutes",
            "longitude_error",
            "depth_km",
            "depth_decimal_km",
            "depth_error",
            "magnitude_1",
            "magnitude_type_1",
            "magnitude_2",
            "magnitude_type_2",
            "travel_time_table",
            "location_precision",
            "subsidiary_info",
            "max_intensity",
            "damage_class",
            "tsunami_class",
            "district_number",
            "region_number",
            "region_name",
            "num_stations",
            "determination_flag",
        ],
        dtype={
            "record_type": str,
            "year": int,
            "month": int,
            "day": int,
            "hour": int,
            "minute": int,
            "second": int,
            "millisecond": int,
            "second_error": float,
            "latitude_degrees": str,
            "latitude_minutes": float,
            "latitude_error": float,
            "longitude_degrees": str,
            "longitude_minutes": float,
            "longitude_error": float,
            "depth_km": float,
            "depth_decimal_km": float,
            "depth_error": float,
            "magnitude_1": str,  # Keep as string for decoding
            "magnitude_type_1": str,
            "magnitude_2": str,  # Keep as string for decoding
            "magnitude_type_2": str,
            "travel_time_table": str,
            "location_precision": str,
            "subsidiary_info": str,
            "max_intensity": str,
            "damage_class": str,
            "tsunami_class": str,
            "district_number": str,
            "region_number": str,
            "region_name": str,
            "num_stations": str,
            "determination_flag": str,
        },
    )
    
    # Choosing station reference coordinate
    if station_name == "FUJ":
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_FUJ
    elif station_name == "MMS2": 
        
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_MMS2
   
    print("done")

    print("Converting quantities... ", end="")
    # Decode magnitude columns first (before numeric conversion)
    df["magnitude_1"] = df["magnitude_1"].apply(decode_magnitude)
    df["magnitude_2"] = df["magnitude_2"].apply(decode_magnitude)

    # Combine second and millisecond into a single column with decimal seconds
    df["second"] = df.apply(
        lambda row: seconds(row["second"], row["millisecond"]), axis=1
    )
    # Divides depth by 100 (some depth are bigger than earth radius)
    # df["depth_km"] = df["depth_km"].apply(convert_depth)

    df["depth_decimal_km"] = df["depth_decimal_km"].fillna(0)
    df["depth_km"] = df["depth_km"] + (df["depth_decimal_km"] / 100)

    # Turn numeric columns into numbers
    for col in [
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "second_error",
        "latitude_degrees",
        "latitude_minutes",
        "latitude_error",
        "longitude_degrees",
        "longitude_minutes",
        "longitude_error",
        "depth_km",
        "depth_error",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Convert latitude and longitude to decimal degrees
    df["latitude"] = df.apply(
        lambda row: degrees_minutes_to_decimal(
            row["latitude_degrees"], row["latitude_minutes"]
        ),
        axis=1,
    )
    df["longitude"] = df.apply(
        lambda row: degrees_minutes_to_decimal(
            row["longitude_degrees"], row["longitude_minutes"]
        ),
        axis=1,
    )

    # Turn the year, month, day, hour, minute, second into a datetime

    df["datetime"] = pd.to_datetime(
        df[["year", "month", "day", "hour", "minute", "second"]]
    )
    print("done")

    # Drop the original columns
    df.drop(
        [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "millisecond",
            "latitude_degrees",
            "latitude_minutes",
            "longitude_degrees",
            "longitude_minutes",
            "longitude_error",
            "latitude_error",
            "depth_error",
            "second_error",
            "magnitude_2",
            "magnitude_type_2",
            "magnitude_type_1",
            "travel_time_table",
            "location_precision",
            "damage_class",
            "tsunami_class",
            "district_number",
            "region_number",
            "region_name",
            "determination_flag",
            "num_stations",
            "max_intensity",
        ],
        axis=1,
        inplace=True,
    )

    # Use it as the index
    df.set_index("datetime", inplace=True)

    # Convert subsidiary information codes to strings
    df["subsidiary_info"] = df["subsidiary_info"].apply(
        subsidiary_info_to_string
    )
    # df.drop("subsidiary_info", axis=1, inplace=True)

    # Drop na
    df.dropna(inplace=True)

    if start is not None:
        start = pd.to_datetime(start)
        df = df[df.index >= start]
    if end is not None:
        end = pd.to_datetime(end)
        df = df[df.index <= end]


    if station_name is not None:


        df["epicentral_distance"] = df.apply(
            lambda row: obspy.geodetics.locations2degrees(
                lat1=row["latitude"],
                long1=row["longitude"],
                lat2=REFERENCE_COORDINATE[0],
                long2=REFERENCE_COORDINATE[1],
            ),
            axis=1,
            result_type="expand",
        )


        if max_distance is not None:
            df = df[df.epicentral_distance.between(0, max_distance)]

        if max_depth is not None:
            df = df[df.depth_km.between(0, max_depth)]

        if min_magnitude is not None:
            df = df[df.magnitude_1.between(min_magnitude, df.magnitude_1.max())]

        print("Calculating travel times... ", end="")
        if all_phases == False:
            df[["p_time", "s_time", "r_time"]] = df.apply(
                lambda row: calculate_P_S_travel_time(
                    row["longitude"], row["latitude"], row["depth_km"],REFERENCE_COORDINATE
                ),
                axis=1,
                result_type="expand",
            )

            # #Calculate arrival times of P and S waves
            df[["arrival_p_time", "arrival_s_time", "arrival_r_time"]] = df.apply(
                lambda row: arrival_time(
                    row.name, row["p_time"], row["s_time"], row["r_time"]
                ),
                axis=1,
                result_type="expand",
            )
        print("done")

        if LFE == True:
            df = df[df["subsidiary_info"] == 5.0]



    return df

    # def best_contenders(data,max_distance=0.5, max_depth=50.0, min_magnitude=3.0):
    """Select earthquakes that are within a certain distance, depth, and magnitude range."""
    if data is not None:
        df = data

    df.dropna(inplace=True)
    # Calculate distance using obspy geodetics
    df["epicentral_distance"] = df.apply(
        lambda row: obspy.geodetics.locations2degrees(
            lat1=row["latitude"],
            long1=row["longitude"],
            lat2=REFERENCE_COORDINATE[0],
            long2=REFERENCE_COORDINATE[1],
        ),
        axis=1,
        result_type="expand",
    )

    # Create a mask to select earthquakes based on location and magnitude
    mask = df[
        (df.epicentral_distance.between(0, max_distance))
        & (df.depth_km.between(0, max_depth))
        & (df.magnitude_1.between(min_magnitude, df.magnitude_1.max()))
    ]

    # Calculate arrival times of P and S waves for the selected earthquakes
    mask[["p_time", "s_time"]] = mask.apply(
        lambda row: calculate_P_S_travel_time(
            row["longitude"], row["latitude"], row["depth_km"]
        ),
        axis=1,
        result_type="expand",
    )

    # Calculate arrival times of P and S waves
    mask[["arrival_p_time", "arrival_s_time"]] = mask.apply(
        lambda row: arrival_time(row.name, row["p_time"], row["s_time"]),
        axis=1,
        result_type="expand",
    )
    # mask.dropna(inplace=True)

    return mask

    # def picking_LFE(max_distance=1.0, max_depth=50.0, min_magnitude=1.0):
    """Select earthquakes that are within a certain distance, depth, and magnitude range."""
    df = read_catalog()
    df.dropna(inplace=True)
    # Calculate distance using obspy geodetics
    df["epicentral_distance"] = df.apply(
        lambda row: obspy.geodetics.locations2degrees(
            lat1=row["latitude"],
            long1=row["longitude"],
            lat2=REFERENCE_COORDINATE[0],
            long2=REFERENCE_COORDINATE[1],
        ),
        axis=1,
        result_type="expand",
    )

    # Create a mask to select earthquakes based on location and magnitude
    mask = df[
        (df.epicentral_distance.between(0, max_distance))
        & (df.depth_km.between(0, max_depth))
        & (df.magnitude_1.between(min_magnitude, df.magnitude_1.max()))
    ]
    mask = mask[mask["subsidiary_info"] == 5.0]

    # Calculate arrival times of P and S waves for the selected earthquakes
    mask[["p_time", "s_time"]] = mask.apply(
        lambda row: calculate_P_S_travel_time(
            row["longitude"], row["latitude"], row["depth_km"]
        ),
        axis=1,
        result_type="expand",
    )

    # Calculate arrival times of P and S waves
    mask[["arrival_p_time", "arrival_s_time"]] = mask.apply(
        lambda row: arrival_time(row.name, row["p_time"], row["s_time"]),
        axis=1,
        result_type="expand",
    )
    # mask.dropna(inplace=True)

    return mask

    # def between_dates(start, end):
    df = read_catalog()
    df.dropna(inplace=True)
    # Select a time window (if provided)

    if start is not None:
        start = pd.to_datetime(start)
        df = df[df.index >= start]
    if end is not None:
        end = pd.to_datetime(end)
        df = df[df.index <= end]

    df[["p_time", "s_time"]] = df.apply(
        lambda row: calculate_P_S_travel_time(
            row["longitude"], row["latitude"], row["depth_km"]
        ),
        axis=1,
        result_type="expand",
    )

    # #Calculate arrival times of P and S waves
    df[["arrival_p_time", "arrival_s_time"]] = df.apply(
        lambda row: arrival_time(row.name, row["p_time"], row["s_time"]),
        axis=1,
        result_type="expand",
    )

    # print(mask.describe())

    return df


def plot_catalog(
    df,
    ax=None,
    projection=cartopy.crs.PlateCarree(),
    colors={"J": "blue", "U": "orange"},
    **kwargs,
):

    # Create figure if not provided
    if ax is None:
        _, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Scatter plot default parameters
    kwargs.setdefault("s", 1)
    kwargs.setdefault("alpha", 0.6)
    kwargs.setdefault("c", df.record_type.map(colors))
    ax.scatter(df.longitude, df.latitude, **kwargs)

    # Add coastlines
    ax.coastlines()
    return ax


def plot_section(df):
    # See info
    mask = df[
        df.longitude.between(130, 150) & df.latitude.between(35, 35.5)
    ]  # & df.depth_km.between(0,4000)
    print(mask.describe())

    fig = plt.figure()
    plt.plot(mask.longitude, -mask.depth_km, ".")
    plt.plot(138.725263766, 0, "^r")
    plt.xlabel("Longitude")
    plt.ylabel("Profondeur en km")
    return fig, plt.gca()


def map_waveform(df, freq,station_name,LFE = False):
    catalog = df
    # Putting the station name in the path to read the seismic data
    # if station_name == "FUJ":
    #     PATH_SEISMIC_DATA = "/gpfs/scratch/doucet/seismograms/FUJI/mseed/2008/{month}/{day}/EV.FUJ..w*"
    # elif station_name == "MMS":
    #     PATH_SEISMIC_DATA = "/gpfs/scratch/doucet/seismograms/FUJI/mseed/2008/{month}/{day}/EV.MMS2..w*"
    # else:
    #     raise ValueError("Station name must be either 'FUJ' or 'MMS'")
    
    PATH_SEISMIC_DATA = "/gpfs/scratch/doucet/seismograms/FUJI/mseed/2008/{month}/{day}/EV.{station_name}..w*"

    # Choosing station reference coordinate
    if station_name == "FUJ":
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_FUJ  # FUJ coordinates
    elif station_name == "MMS2":
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_MMS2  # MMS coordinates

    
    if LFE == True: 
        save_path = f"figures/LFE/{freq}Hz/{station_name}/"

    else :
        save_path = f"figures/best_contenders/{freq}Hz/{station_name}/"

    for i in range(len(catalog)):
        # Plot catalog where each subplot is a different channel from the stream
        fig, ax = plt.subplot_mosaic(
            [
                ["top"],
                ["{stream[0].stats.channel}"],
                ["{stream[0].spectrogram}"],
                ["{stream[1].stats.channel}"],
                ["{stream[1].spectrogram}"],
                ["{stream[2].stats.channel}"],
                ["{stream[2].spectrogram}"],
            ],
            figsize=(10, 20),
            gridspec_kw={"height_ratios": [3, 1, 2, 1, 2, 1, 2]},
            per_subplot_kw={"top": {"projection": cartopy.crs.PlateCarree()}},
            constrained_layout=True,
        )
        plot_catalog(
            catalog,
            ax=ax["top"],
            s=catalog.magnitude_1.clip(1, 5) * 100,
            c=catalog.p_time,
            linewidth=0,
            alpha=1,
            cmap="viridis",
            vmin=0,
        )

        # Add station
        station_coordinates = REFERENCE_COORDINATE
        ax["top"].plot(
            station_coordinates[1],
            station_coordinates[0],
            marker="v",
            color="w",
            mec="k",
            markersize=6,
            label="Station",
        )

        # Labels
        mappable = ax["top"].collections[0]
        cbar = plt.colorbar(mappable, ax=ax["top"], shrink=0.5)
        cbar.set_label("$P$ travel time (s)")

        # Sort events as a function of p_time
        catalog = catalog.sort_values("p_time")
        event = catalog.iloc[i]
        # if event.magnitude_1 < 1.5:
        #     print(f"Skipping event {event.name} with magnitude {event.magnitude_1}")
        #     continue

        print("Reading seismic data... ", end="")
        p_time = obspy.UTCDateTime(event.arrival_p_time)
        s_time = obspy.UTCDateTime(event.arrival_s_time)
        stream = obspy.read(
            PATH_SEISMIC_DATA.format(
                month=str(event.name.month).zfill(2),
                day=str(event.name.day).zfill(2),
                station_name=station_name,
            ),
            starttime=obspy.UTCDateTime(event.name) - 20,
            endtime=obspy.UTCDateTime(event.name) + 100,
        )

        spectro = stream.copy()

        # The three magic steps to make the seismogram look nice
        stream.detrend("demean")
        stream.taper(max_percentage=0.05)
        stream.filter("highpass", freq=freq)

        # Plot event on map
        ax["top"].plot(
            event.longitude,
            event.latitude,
            marker="*",
            color="r",
            label="Largest event",
        )
        ax["top"].set_extent(
            [
                event.longitude.min() - 1,
                event.longitude.max() + 1,
                event.latitude.min() - 1,
                event.latitude.max() + 1,
            ]
        )

        # Plot seismograms
        for i in range(len(stream)):
            ax[f"{{stream[{i}].stats.channel}}"].plot(
                stream[i].times(),
                stream[i].data,
                color="k",
                linewidth=0.3,
            )
            ax[f"{{stream[{i}].stats.channel}}"].axvline(
                p_time - stream[i].stats.starttime,
                color="r",
                label="P arrival",
            )
            ax[f"{{stream[{i}].stats.channel}}"].axvline(
                s_time - stream[i].stats.starttime,
                color="b",
                label="S arrival",
            )
            spectro[i].spectrogram(
                axes=ax[f"{{stream[{i}].spectrogram}}"],
                log=True,
                wlen=2,
                per_lap=0.75,
                dbscale=True,
                clip=[0.2, 0.9],
            )
            ax[f"{{stream[{i}].stats.channel}}"].text(
                p_time - stream[i].stats.starttime - 5,
                max(stream[i].data) * 0.8,
                "P",
                color="r",
            )
            ax[f"{{stream[{i}].stats.channel}}"].text(
                s_time - stream[i].stats.starttime + 2,
                max(stream[i].data) * 0.8,
                "S",
                color="b",
            )

            ax[f"{{stream[{i}].stats.channel}}"].set_ylim(
                -max(abs(stream[i].data)), max(abs(stream[i].data))
            )
            ax[f"{{stream[{i}].stats.channel}}"].grid()
            ax[f"{{stream[{i}].stats.channel}}"].set_xlabel(
                f"Distance = {event.epicentral_distance:.2f} (°), M={event.magnitude_1} (s) ,channel {stream[i].stats.channel}"
            )
            ax[f"{{stream[{i}].stats.channel}}"].set_ylabel("Amplitude")

            for key in (
                f"{{stream[{i}].stats.channel}}",
                f"{{stream[{i}].spectrogram}}",
            ):
                ax[key].set_xlim(0, 120)

        # Save figure

        os.makedirs(save_path, exist_ok=True)
        fig.savefig(save_path + f"{event.name}.png", dpi=300)
        plt.close(fig)


# if __name__ == "read_jma":
#     map_waveform()


def calculate_all_phases(lon, lat, depth,REFERENCE_COORDINATE):
    """Calculate all phases tt with taup"""

    model = TauPyModel(model="iasp91")
    epicentral_distance = obspy.geodetics.locations2degrees(
        lat1=lat,
        long1=lon,
        lat2=REFERENCE_COORDINATE[0],
        long2=REFERENCE_COORDINATE[1],
    )
    try:
        arrivals = model.get_travel_times(
            source_depth_in_km=depth,
            distance_in_degree=epicentral_distance,
        )
        return arrivals
    except Exception as e:
        print(
            f"Error calculating travel times for ({lon}, {lat}, {depth}): {e}"
        )

        return None


def all_phases(df, save_path, station_name):
    # Applying the function to each row of the DataFrame and returning the arrivals object for each earthquake

    arr = df.apply(
        lambda row: calculate_all_phases(
            row["longitude"], row["latitude"], row["depth_km"]
        ),
        axis=1,
    )

    arr = arr.explode().reset_index(drop=True)

    catalog = df
    PATH_SEISMIC_DATA = (
        "/gpfs/scratch/doucet/seismograms/FUJI/mseed/2008/*/*/EV.MMS2..w*"
    )
    # save_path = f"figures/all_phases/"
    for i in range(len(catalog)):

        fig, ax = plt.subplot_mosaic(
            [
                ["top"],
                ["{stream[0].stats.channel}"],
                ["{stream[0].spectrogram}"],
                ["{stream[1].stats.channel}"],
                ["{stream[1].spectrogram}"],
                ["{stream[2].stats.channel}"],
                ["{stream[2].spectrogram}"],
            ],
            figsize=(10, 20),
            gridspec_kw={"height_ratios": [3, 1, 2, 1, 2, 1, 2]},
            per_subplot_kw={"top": {"projection": cartopy.crs.PlateCarree()}},
            constrained_layout=True,
        )
    plot_catalog(
        catalog,
        ax=ax["top"],
        s=catalog.magnitude_1.clip(1, 5) * 100,
        c=catalog.p_time,
        linewidth=0,
        alpha=1,
        cmap="viridis",
        vmin=0,
    )

    if station_name == "FUJ":
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_FUJ
    elif station_name == "MMS2": 
        
        REFERENCE_COORDINATE = REFERENCE_COORDINATE_MMS2

    # Add station
    station_coordinates = REFERENCE_COORDINATE
    ax["top"].plot(
        station_coordinates[1],
        station_coordinates[0],
        marker="v",
        color="w",
        mec="k",
        markersize=6,
        label="Station",
    )

    mappable = ax["top"].collections[0]
    cbar = plt.colorbar(mappable, ax=ax["top"], shrink=0.5)
    cbar.set_label("$P$ travel time (s)")

    # Sort events as a function of p_time
    catalog = catalog.sort_values("p_time")
    event = catalog.iloc[i]
    # if event.magnitude_1 < 1.5:
    #     print(f"Skipping event {event.name} with magnitude {event.magnitude_1}")
    #     continue
    starttime = obspy.UTCDateTime(event.name) - 20
    endtime = obspy.UTCDateTime(event.name) + 100
    p_time = obspy.UTCDateTime(event.arrival_p_time)
    s_time = obspy.UTCDateTime(event.arrival_s_time)
    stream = obspy.read(
        PATH_SEISMIC_DATA, starttime=starttime, endtime=endtime
    )

    diff = obspy.UTCDateTime(event.name) - starttime

    spectro = stream.copy()

    stream.detrend("demean")
    stream.taper(max_percentage=0.05)
    stream.filter("highpass", freq=1.0)

    # Plot event on map
    ax["top"].plot(
        event.longitude,
        event.latitude,
        marker="*",
        color="r",
        label="Largest event",
    )
    ax["top"].set_extent(
        [
            event.longitude.min() - 1,
            event.longitude.max() + 1,
            event.latitude.min() - 1,
            event.latitude.max() + 1,
        ]
    )

    # Plot seismograms
    for i in range(len(stream)):
        ax[f"{{stream[{i}].stats.channel}}"].plot(
            stream[i].times(),
            stream[i].data,
            color="k",
            linewidth=0.3,
        )

        for j, phase in enumerate(arr):
            if phase.time + diff < endtime - starttime:
                ax[f"{{stream[{i}].stats.channel}}"].axvline(
                    diff + phase.time,
                    color=f"C{j}",
                    label="{phase.name} arrival",
                    alpha=0.5,
                )
                ax[f"{{stream[{i}].stats.channel}}"].text(
                    diff + phase.time + 2,
                    max(stream[i].data) * 0.8,
                    phase.name,
                    color=f"C{j}",
                    alpha=0.5,
                )

            else:
                continue

        spectro[i].spectrogram(
            axes=ax[f"{{stream[{i}].spectrogram}}"],
            log=True,
            wlen=10,
            per_lap=0.5,
        )
        # ax[f"{{stream[{i}].stats.channel}}"].text(p_time - stream[i].stats.starttime - 5, max(stream[i].data) * 0.8, "P", color="r")

    os.makedirs(save_path, exist_ok=True)
    fig.savefig(save_path + f"{event.name}.png", dpi=300)
    plt.close(fig)

    return arr

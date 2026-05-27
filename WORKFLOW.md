# Fuji Continuous Data Workflow

`SonifSismo` provides an initial Julia interface to the Fuji daily miniSEED
archive and JMA fixed-width event catalog.

## Time Convention

The JMA `h2008` catalog is interpreted as Japan Standard Time (`UTC+09:00`) by
default. `JMAEvent.catalog_time` preserves that printed clock time, while
`JMAEvent.origin_time_utc` is the canonical time for miniSEED loading and
trimming. The miniSEED files inspected in this archive begin at `00:00:00`
UTC on their filename date.

If a future catalog is UTC, construct the archive with:

```julia
archive = DataArchive("/path/to/archive"; catalog_timezone=:UTC)
```

## First Exploration

```julia
using Dates
using SonifSismo

archive = DataArchive("/Users/nobuaki/Documents/mtFujiContinuous")
available_stations(archive)
available_channels(archive)

events = read_catalog(
    archive;
    start=DateTime(2008, 5, 25),
    stop=DateTime(2008, 5, 31, 23, 59, 59),
    min_magnitude=2.0,
    max_depth=50.0,
)

traces = read_event(
    archive,
    first(events);
    stations="FUJ",
    channels=["wE", "wN", "wU"],
    processed=true,
    bandpass=(0.5, 15.0),
)
```

Use an arbitrary UTC interval when you do not want to select from the catalog:

```julia
traces = read_window(
    archive,
    DateTime(2008, 5, 27, 12),
    DateTime(2008, 5, 27, 12, 1);
    stations="FUJ",
    channels="wU",
)
```

## Filtering And Time-Frequency Data

`suggest_filters(traces)` supplies conservative exploratory starting bands
bounded by each trace's Nyquist frequency. `process_traces` returns copies, so
the raw waveforms remain available for comparison and eventual sonification.

```julia
recipes = suggest_filters(traces)
filtered = process_traces(traces; bandpass=(0.5, 15.0))
spec = time_frequency(first(filtered); window=2.0, overlap=0.75)
```

## Makie Plots

Plotting is optional. Load the CairoMakie backend in the notebook:

```julia
using CairoMakie, SonifSismo

plot_events = read_catalog(
    archive;
    start=first(events).origin_time_utc - Second(20),
    stop=first(events).origin_time_utc + Second(100),
)

fig = plot_waveforms(traces; events=plot_events)
save("waveforms.png", fig)

specfig = plot_time_frequency(traces; events=plot_events)
save("spectrogram.png", specfig)
```

`plot_waveforms` and `plot_time_frequency` draw catalog markers for the
events within the loaded trace window, with different colors for earthquake,
low-frequency earthquake (LFE), artificial, eruption/other, and unknown
records. Their horizontal axes show seconds after the displayed window's UTC
start time. To display rare LFEs in context, load a longer UTC waveform window
and pass a less restrictive `plot_events` query than the one used to select a
featured earthquake.

The waveform-loading and processing functions intentionally return `Seis`
traces. They can later feed resampling and WAV export without changing the
archive or event-selection layer.

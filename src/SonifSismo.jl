module SonifSismo

using Dates
using Seis
using Statistics

export DataArchive,
       FujiArchive,
       JMAEvent,
       ChannelAvailability,
       ProcessingRecipe,
       event_type,
       event_type_label,
       catalog_time_to_utc,
       parse_jma_event,
       read_catalog,
       available_channels,
       available_stations,
       waveform_files,
       read_window,
       read_event,
       suggest_filters,
       process_traces,
       time_frequency,
       plot_waveforms,
       plot_time_frequency

const JST_UTC_OFFSET = Hour(9)
const MSEED_FILENAME = r"^(?<network>[^.]+)\.(?<station>[^.]+)\.(?<location>[^.]*)\.(?<channel>[^.]+)\.(?<date>\d{8})\.mseed$"

"""
    DataArchive(root; waveform_dir="mseed", catalog="catalog/h2008",
                catalog_timezone=:JST, catalog_utc_offset_hours=nothing)

Description of a Fuji continuous-data archive. Times passed to waveform reading
functions are UTC `DateTime`s. By default JMA catalog origin times are treated
as Japan Standard Time (UTC+09:00) and converted to UTC when read.

Pass `catalog_timezone=:UTC` if a catalog is already in UTC, or provide an
explicit `catalog_utc_offset_hours` for another catalog convention.
"""
struct DataArchive
    root::String
    waveform_root::String
    catalog_path::String
    catalog_timezone::Symbol
    catalog_utc_offset::Hour
end

function DataArchive(
    root::AbstractString;
    waveform_dir::AbstractString="mseed",
    catalog::AbstractString=joinpath("catalog", "h2008"),
    catalog_timezone::Symbol=:JST,
    catalog_utc_offset_hours::Union{Nothing,Integer}=nothing,
)
    offset = if !isnothing(catalog_utc_offset_hours)
        Hour(catalog_utc_offset_hours)
    elseif catalog_timezone === :JST
        JST_UTC_OFFSET
    elseif catalog_timezone === :UTC
        Hour(0)
    else
        throw(ArgumentError(
            "unknown catalog_timezone=$catalog_timezone; use :JST, :UTC, " *
            "or provide catalog_utc_offset_hours",
        ))
    end
    return DataArchive(
        abspath(root),
        abspath(joinpath(root, waveform_dir)),
        abspath(joinpath(root, catalog)),
        catalog_timezone,
        offset,
    )
end

const FujiArchive = DataArchive

"""
A parsed JMA hypocentre record. `catalog_time` retains the clock printed in
the input file and `origin_time_utc` is the UTC time used for waveform reads.
"""
struct JMAEvent
    record_type::Char
    catalog_time::DateTime
    origin_time_utc::DateTime
    latitude::Union{Missing,Float64}
    longitude::Union{Missing,Float64}
    depth_km::Union{Missing,Float64}
    magnitude::Union{Missing,Float64}
    magnitude_type::Union{Missing,String}
    subsidiary_code::Union{Missing,Int}
    region_name::String
    raw::String
end

"""
Availability for one network/station/location/channel combination.
`days` and `files` refer to daily miniSEED files, whose dates are UTC.
"""
struct ChannelAvailability
    network::String
    station::String
    location::String
    channel::String
    days::Vector{Date}
    files::Vector{String}
end

struct ProcessingRecipe
    name::Symbol
    description::String
    bandpass::Union{Nothing,Tuple{Float64,Float64}}
    highpass::Union{Nothing,Float64}
    lowpass::Union{Nothing,Float64}
end

const EVENT_TYPE_LABELS = Dict(
    :earthquake => "Earthquake",
    :insufficient_information => "Insufficient information",
    :artificial => "Artificial event",
    :eruption_other => "Eruption/other event",
    :lfe => "Low-frequency earthquake (LFE)",
    :unknown => "Unknown event type",
)

"""
    event_type(event) -> Symbol

Classify a JMA event from its subsidiary information code. The JMA categories
used here are `:earthquake`, `:insufficient_information`, `:artificial`,
`:eruption_other`, `:lfe`, and `:unknown`.
"""
function event_type(event::JMAEvent)
    ismissing(event.subsidiary_code) && return :unknown
    return get(
        Dict(1 => :earthquake, 2 => :insufficient_information, 3 => :artificial,
             4 => :eruption_other, 5 => :lfe),
        event.subsidiary_code,
        :unknown,
    )
end

event_type_label(event::JMAEvent) = EVENT_TYPE_LABELS[event_type(event)]

catalog_time_to_utc(archive::DataArchive, time::DateTime) =
    time - archive.catalog_utc_offset

function _slice(line::AbstractString, first::Int, last::Int)
    ncodeunits(line) < first && return ""
    return strip(SubString(line, first, min(last, ncodeunits(line))))
end

function _parse_number(::Type{T}, value::AbstractString) where {T}
    isempty(value) && return missing
    return tryparse(T, value) === nothing ? missing : parse(T, value)
end

function _decode_magnitude(value::AbstractString)
    value = strip(value)
    isempty(value) && return missing
    if length(value) >= 2 && first(value) in ('A', 'B', 'C')
        digit = tryparse(Int, string(value[2]))
        digit === nothing && return missing
        base = Dict('A' => -1.0, 'B' => -2.0, 'C' => -3.0)[first(value)]
        return base - digit / 10
    end
    raw = tryparse(Int, value)
    return raw === nothing ? missing : raw / 10
end

function _coordinate(degrees::Union{Missing,Int}, minutes_hundredths::Union{Missing,Float64})
    (ismissing(degrees) || ismissing(minutes_hundredths)) && return missing
    return degrees + minutes_hundredths / 6000
end

"""
    parse_jma_event(line, archive) -> Union{JMAEvent,Nothing}

Parse one fixed-width JMA hypocentre line. Non-event or malformed lines return
`nothing`, allowing a mixed catalog file to be consumed robustly.
"""
function parse_jma_event(line::AbstractString, archive::DataArchive)
    isempty(line) && return nothing
    record_type = first(line)
    record_type == 'J' || return nothing

    components = (
        _parse_number(Int, _slice(line, 2, 5)),
        _parse_number(Int, _slice(line, 6, 7)),
        _parse_number(Int, _slice(line, 8, 9)),
        _parse_number(Int, _slice(line, 10, 11)),
        _parse_number(Int, _slice(line, 12, 13)),
        _parse_number(Int, _slice(line, 14, 15)),
        _parse_number(Int, _slice(line, 16, 17)),
    )
    any(ismissing, components) && return nothing
    year, month, day, hour, minute, second, centisecond = components
    catalog_time = try
        DateTime(year, month, day, hour, minute, second, 10centisecond)
    catch
        return nothing
    end

    lat = _coordinate(
        _parse_number(Int, _slice(line, 22, 24)),
        _parse_number(Float64, _slice(line, 25, 28)),
    )
    lon = _coordinate(
        _parse_number(Int, _slice(line, 33, 36)),
        _parse_number(Float64, _slice(line, 37, 40)),
    )
    depth_whole = _parse_number(Float64, _slice(line, 45, 47))
    depth_decimal = _parse_number(Float64, _slice(line, 48, 49))
    depth = ismissing(depth_whole) ? missing :
        depth_whole + (ismissing(depth_decimal) ? 0.0 : depth_decimal / 100)
    magnitude = _decode_magnitude(_slice(line, 53, 54))
    magnitude_type = let value = _slice(line, 55, 55)
        isempty(value) ? missing : value
    end
    subsidiary_code = _parse_number(Int, _slice(line, 61, 61))
    region_name = _slice(line, 69, 92)

    return JMAEvent(
        record_type,
        catalog_time,
        catalog_time_to_utc(archive, catalog_time),
        lat,
        lon,
        depth,
        magnitude,
        magnitude_type,
        subsidiary_code,
        region_name,
        String(line),
    )
end

function _query_time_utc(archive::DataArchive, value::Union{Nothing,DateTime}, input_timezone::Symbol)
    isnothing(value) && return nothing
    if input_timezone === :UTC
        return value
    elseif input_timezone === :catalog || input_timezone === archive.catalog_timezone
        return catalog_time_to_utc(archive, value)
    end
    throw(ArgumentError("input_timezone must be :UTC, :catalog, or $(archive.catalog_timezone)"))
end

"""
    read_catalog(archive; start=nothing, stop=nothing, input_timezone=:UTC,
                 min_magnitude=nothing, max_depth=nothing, predicate=nothing)

Return JMA events matching the requested interval and filters. `start` and
`stop` are interpreted as UTC by default; pass `input_timezone=:catalog` when
selecting in the catalog's printed clock.
"""
function read_catalog(
    archive::DataArchive;
    start::Union{Nothing,DateTime}=nothing,
    stop::Union{Nothing,DateTime}=nothing,
    input_timezone::Symbol=:UTC,
    min_magnitude::Union{Nothing,Real}=nothing,
    max_depth::Union{Nothing,Real}=nothing,
    predicate::Union{Nothing,Function}=nothing,
)
    isfile(archive.catalog_path) ||
        throw(ArgumentError("catalog file does not exist: $(archive.catalog_path)"))
    start_utc = _query_time_utc(archive, start, input_timezone)
    stop_utc = _query_time_utc(archive, stop, input_timezone)
    (!isnothing(start_utc) && !isnothing(stop_utc) && stop_utc < start_utc) &&
        throw(ArgumentError("stop time is before start time"))

    events = JMAEvent[]
    open(archive.catalog_path, "r") do io
        for line in eachline(io)
            event = parse_jma_event(line, archive)
            isnothing(event) && continue
            !isnothing(start_utc) && event.origin_time_utc < start_utc && continue
            !isnothing(stop_utc) && event.origin_time_utc > stop_utc && continue
            !isnothing(min_magnitude) &&
                (ismissing(event.magnitude) || event.magnitude < min_magnitude) && continue
            !isnothing(max_depth) &&
                (ismissing(event.depth_km) || event.depth_km > max_depth) && continue
            !isnothing(predicate) && !predicate(event) && continue
            push!(events, event)
        end
    end
    return events
end

function _mseed_metadata(file::AbstractString)
    match_result = match(MSEED_FILENAME, basename(file))
    isnothing(match_result) && return nothing
    date = try
        Date(match_result[:date], dateformat"yyyymmdd")
    catch
        return nothing
    end
    return (
        network=String(match_result[:network]),
        station=String(match_result[:station]),
        location=String(match_result[:location]),
        channel=String(match_result[:channel]),
        date=date,
        file=String(file),
    )
end

function _all_mseed_files(archive::DataArchive)
    isdir(archive.waveform_root) ||
        throw(ArgumentError("waveform directory does not exist: $(archive.waveform_root)"))
    files = String[]
    for (root, _, names) in walkdir(archive.waveform_root)
        append!(files, joinpath.(Ref(root), filter(name -> endswith(name, ".mseed"), names)))
    end
    return sort!(files)
end

"""
    available_channels(archive; start=nothing, stop=nothing)

Inspect filenames to return the available station/channel combinations quickly,
without loading full daily waveforms. Date filters refer to UTC waveform dates.
"""
function available_channels(
    archive::DataArchive;
    start::Union{Nothing,Date}=nothing,
    stop::Union{Nothing,Date}=nothing,
)
    grouped = Dict{NTuple{4,String},Vector{Tuple{Date,String}}}()
    for file in _all_mseed_files(archive)
        meta = _mseed_metadata(file)
        isnothing(meta) && continue
        !isnothing(start) && meta.date < start && continue
        !isnothing(stop) && meta.date > stop && continue
        key = (meta.network, meta.station, meta.location, meta.channel)
        push!(get!(grouped, key, Tuple{Date,String}[]), (meta.date, meta.file))
    end
    out = ChannelAvailability[]
    for (key, items) in grouped
        sort!(items)
        network, station, location, channel = key
        push!(out, ChannelAvailability(
            network, station, location, channel, first.(items), last.(items),
        ))
    end
    sort!(out; by=item -> (item.station, item.channel, item.network, item.location))
    return out
end

available_stations(archive::DataArchive; kwargs...) =
    sort!(unique(item.station for item in available_channels(archive; kwargs...)))

_selected(value::AbstractString, choice::Nothing) = true
_selected(value::AbstractString, choice::AbstractString) = value == choice
_selected(value::AbstractString, choices) = value in choices

"""
    waveform_files(archive, start_utc, stop_utc; stations=nothing, channels=nothing)

Find daily miniSEED files overlapping a UTC interval, with optional exact
station and channel selections.
"""
function waveform_files(
    archive::DataArchive,
    start_utc::DateTime,
    stop_utc::DateTime;
    stations=nothing,
    channels=nothing,
)
    stop_utc >= start_utc || throw(ArgumentError("stop time is before start time"))
    files = String[]
    for item in available_channels(archive)
        _selected(item.station, stations) || continue
        _selected(item.channel, channels) || continue
        for (date, file) in zip(item.days, item.files)
            day_start = DateTime(date)
            day_stop = day_start + Day(1)
            day_start <= stop_utc && start_utc < day_stop && push!(files, file)
        end
    end
    return sort!(files)
end

"""
    read_window(archive, start_utc, stop_utc; stations=nothing, channels=nothing,
                processed=false, process_kwargs...)

Read daily miniSEED files overlapping a UTC interval and trim each returned
continuous trace to that interval. Set `processed=true` to apply
`process_traces` with the supplied processing keyword arguments.
"""
function read_window(
    archive::DataArchive,
    start_utc::DateTime,
    stop_utc::DateTime;
    stations=nothing,
    channels=nothing,
    maximum_gap=nothing,
    processed::Bool=false,
    process_kwargs...,
)
    files = waveform_files(archive, start_utc, stop_utc; stations=stations, channels=channels)
    traces = Seis.AbstractTrace[]
    for file in files
        input = isnothing(maximum_gap) ? Seis.read_mseed(file) :
            Seis.read_mseed(file; maximum_gap=maximum_gap)
        for tr in input
            Seis.enddate(tr) < start_utc && continue
            Seis.startdate(tr) > stop_utc && continue
            push!(traces, Seis.cut(tr, start_utc, stop_utc; allowempty=false, warn=false))
        end
    end
    sort!(traces; by=tr -> (
        string(tr.sta.sta), string(tr.sta.cha), Seis.startdate(tr),
    ))
    return processed ? process_traces(traces; process_kwargs...) : traces
end

"""
    read_event(archive, event; pre=Second(20), post=Second(100), kwargs...)

Read a waveform window around an event's UTC origin time.
"""
function read_event(
    archive::DataArchive,
    event::JMAEvent;
    pre::Period=Second(20),
    post::Period=Second(100),
    kwargs...,
)
    return read_window(
        archive,
        event.origin_time_utc - pre,
        event.origin_time_utc + post;
        kwargs...,
    )
end

"""
    suggest_filters(traces) -> Vector{ProcessingRecipe}

Return conservative starting recipes based on the sampling frequency. These
are exploration aids, not event-dependent scientific choices.
"""
function suggest_filters(traces::AbstractVector{<:Seis.AbstractTrace})
    isempty(traces) && return ProcessingRecipe[]
    nyquist = minimum(1 / (2tr.delta) for tr in traces)
    broad_high = min(20.0, 0.8nyquist)
    local_high = min(10.0, 0.8nyquist)
    recipes = ProcessingRecipe[
        ProcessingRecipe(:raw, "demean and taper only", nothing, nothing, nothing),
    ]
    broad_high > 0.5 &&
        push!(recipes, ProcessingRecipe(:broadband, "broad local-event view", (0.5, broad_high), nothing, nothing))
    local_high > 1.0 &&
        push!(recipes, ProcessingRecipe(:low_frequency, "emphasise longer-period motion", (0.2, local_high), nothing, nothing))
    return recipes
end

"""
    process_traces(traces; demean=true, taper_width=0.05, bandpass=nothing,
                   highpass=nothing, lowpass=nothing, poles=4, twopass=true)

Return processed copies of traces. Choose at most one of `bandpass`,
`highpass`, and `lowpass`, with frequencies expressed in Hz.
"""
function process_traces(
    traces::AbstractVector{<:Seis.AbstractTrace};
    demean::Bool=true,
    taper_width::Union{Nothing,Real}=0.05,
    bandpass::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    highpass::Union{Nothing,Real}=nothing,
    lowpass::Union{Nothing,Real}=nothing,
    poles::Integer=4,
    twopass::Bool=true,
)
    count(!isnothing, (bandpass, highpass, lowpass)) <= 1 ||
        throw(ArgumentError("choose at most one filter: bandpass, highpass, or lowpass"))
    output = deepcopy(traces)
    for tr in output
        demean && (Seis.trace(tr) .-= mean(Seis.trace(tr)))
        !isnothing(taper_width) && Seis.taper!(tr, taper_width)
        if !isnothing(bandpass)
            Seis.bandpass!(tr, bandpass[1], bandpass[2]; poles=poles, twopass=twopass)
        elseif !isnothing(highpass)
            Seis.highpass!(tr, highpass; poles=poles, twopass=twopass)
        elseif !isnothing(lowpass)
            Seis.lowpass!(tr, lowpass; poles=poles, twopass=twopass)
        end
    end
    return output
end

"""
    time_frequency(trace; window=2.0, overlap=0.75, kwargs...)

Calculate a Seis/DSP spectrogram, ready for plotting or later sonification
feature extraction.
"""
time_frequency(trace::Seis.AbstractTrace; window::Real=2.0, overlap::Real=0.75, kwargs...) =
    Seis.spectrogram(trace; length=window, overlap=overlap, kwargs...)

function plot_waveforms(args...; kwargs...)
    error("plot_waveforms requires a Makie backend; run `using CairoMakie` before plotting")
end

function plot_time_frequency(args...; kwargs...)
    error("plot_time_frequency requires a Makie backend; run `using CairoMakie` before plotting")
end

end

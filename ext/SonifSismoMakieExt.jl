module SonifSismoMakieExt

using Dates
using Makie
using Seis
using SonifSismo

const EVENT_COLORS = Dict(
    :earthquake => :dodgerblue,
    :lfe => :darkorange,
    :artificial => :purple,
    :eruption_other => :crimson,
    :insufficient_information => :grey50,
    :unknown => :black,
)

function _label(trace)
    return join((trace.sta.net, trace.sta.sta, trace.sta.loc, trace.sta.cha), ".")
end

function _window_start(traces)
    return minimum(Seis.startdate(trace) for trace in traces)
end

_seconds_after(start::DateTime, time::DateTime) = Dates.value(time - start) / 1000

function _sample_times(trace, window_start)
    first_time = _seconds_after(window_start, Seis.startdate(trace))
    return range(first_time; step=trace.delta, length=Seis.nsamples(trace))
end

function _spectrogram_times(trace, seconds, window_start)
    ismissing(trace.evt.time) && error("trace needs an absolute origin time for plotting")
    trace_offset = _seconds_after(window_start, trace.evt.time)
    return trace_offset .+ seconds
end

function _visible_events(traces, events)
    isempty(events) && return SonifSismo.JMAEvent[]
    start = minimum(Seis.startdate(trace) for trace in traces)
    stop = maximum(Seis.enddate(trace) for trace in traces)
    return filter(event -> start <= event.origin_time_utc <= stop, events)
end

function _event_groups(events)
    kinds = (:earthquake, :lfe, :artificial, :eruption_other,
             :insufficient_information, :unknown)
    return [(kind, filter(event -> SonifSismo.event_type(event) === kind, events))
            for kind in kinds if any(event -> SonifSismo.event_type(event) === kind, events)]
end

function _mark_events!(axis, events, window_start; label=false)
    for (_, group) in _event_groups(events)
        kind = SonifSismo.event_type(first(group))
        times = [_seconds_after(window_start, event.origin_time_utc) for event in group]
        if label
            vlines!(
                axis,
                times;
                color=EVENT_COLORS[kind],
                linewidth=1.5,
                label=SonifSismo.event_type_label(first(group)),
            )
        else
            vlines!(axis, times; color=EVENT_COLORS[kind], linewidth=1.5)
        end
    end
end

function SonifSismo.plot_waveforms(
    traces::AbstractVector{<:Seis.AbstractTrace};
    events::AbstractVector{<:SonifSismo.JMAEvent}=SonifSismo.JMAEvent[],
    title::AbstractString="Waveforms",
)
    isempty(traces) && throw(ArgumentError("no traces to plot"))
    visible_events = _visible_events(traces, events)
    window_start = _window_start(traces)
    figure = Figure(size=(1100, 250length(traces)))
    axes = Axis[]
    for (index, trace) in enumerate(traces)
        axis = Axis(
            figure[index, 1],
            ylabel=_label(trace),
            xlabel=index == length(traces) ? "Seconds after $(window_start) UTC" : "",
            title=index == 1 ? title : "",
        )
        lines!(axis, _sample_times(trace, window_start), Seis.trace(trace); color=:black)
        _mark_events!(axis, visible_events, window_start; label=index == 1)
        push!(axes, axis)
    end
    length(axes) > 1 && linkxaxes!(axes...)
    !isempty(visible_events) && axislegend(first(axes); position=:rt)
    return figure
end

function SonifSismo.plot_time_frequency(
    trace::Seis.AbstractTrace;
    kwargs...,
)
    return SonifSismo.plot_time_frequency([trace]; kwargs...)
end

function SonifSismo.plot_time_frequency(
    traces::AbstractVector{<:Seis.AbstractTrace};
    events::AbstractVector{<:SonifSismo.JMAEvent}=SonifSismo.JMAEvent[],
    window::Real=2.0,
    overlap::Real=0.75,
    title::AbstractString="Time-frequency power",
)
    isempty(traces) && throw(ArgumentError("no traces to plot"))
    specs = [SonifSismo.time_frequency(trace; window=window, overlap=overlap) for trace in traces]
    powers = [10 .* log10.(spec.power .+ eps(eltype(spec.power))) for spec in specs]
    colorrange = extrema(vcat(vec.(powers)...))
    visible_events = _visible_events(traces, events)
    window_start = _window_start(traces)
    figure = Figure(size=(1100, 260length(traces)))
    axes = Axis[]
    heatmaps = Any[]
    for (index, (trace, spec, power_db)) in enumerate(zip(traces, specs, powers))
        axis = Axis(
            figure[index, 1],
            ylabel="$(_label(trace))\nFrequency (Hz)",
            xlabel=index == length(traces) ? "Seconds after $(window_start) UTC" : "",
            title=index == 1 ? title : "",
        )
        heatmap = heatmap!(
            axis,
            _spectrogram_times(trace, spec.time, window_start),
            spec.freq,
            permutedims(power_db);
            colorrange=colorrange,
            colormap=:magma,
        )
        _mark_events!(axis, visible_events, window_start; label=index == 1)
        push!(axes, axis)
        push!(heatmaps, heatmap)
    end
    length(axes) > 1 && linkxaxes!(axes...)
    !isempty(visible_events) && axislegend(first(axes); position=:rt)
    Colorbar(figure[1:length(traces), 2], first(heatmaps); label="Power (dB)")
    return figure
end

end

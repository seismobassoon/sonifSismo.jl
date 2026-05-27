using Dates
using SonifSismo

archive = DataArchive("/Users/nobuaki/Documents/mtFujiContinuous")

println("Available stations: ", join(available_stations(archive), ", "))
for availability in available_channels(archive)
    println(
        availability.station,
        " ",
        availability.channel,
        ": ",
        first(availability.days),
        " to ",
        last(availability.days),
    )
end

# Catalog times are converted from JST to UTC by DataArchive.
events = read_catalog(
    archive;
    start=DateTime(2008, 5, 25),
    stop=DateTime(2008, 5, 31, 23, 59, 59),
    min_magnitude=2.0,
    max_depth=50.0,
)
println("Matching May events: ", length(events))

if !isempty(events)
    event = first(events)
    println("First event catalog time: ", event.catalog_time, " JST")
    println("First event waveform time: ", event.origin_time_utc, " UTC")
    traces = read_event(
        archive,
        event;
        stations="FUJ",
        channels=["wE", "wN", "wU"],
        processed=true,
        bandpass=(0.5, 15.0),
    )
    println("Read ", length(traces), " event traces")

    # Optional plotting:
    # using CairoMakie
    # plot_events = read_catalog(
    #     archive;
    #     start=event.origin_time_utc - Second(20),
    #     stop=event.origin_time_utc + Second(100),
    # )
    # display(plot_waveforms(traces; events=plot_events, title="FUJ event waveforms"))
    # display(plot_time_frequency(traces; events=plot_events))
end

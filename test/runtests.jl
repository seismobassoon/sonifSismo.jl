using Dates
using Seis
using SonifSismo
using Test

const DATA_ROOT = "/Users/nobuaki/Documents/mtFujiContinuous"

@testset "JMA time conversion and parsing" begin
    archive = DataArchive(DATA_ROOT)
    line = "J2008010100085321 019 325190 070 1351810 083 334036218V   511   6248S OFF KII PENINSULA      21K"
    event = parse_jma_event(line, archive)
    @test !isnothing(event)
    @test event.catalog_time == DateTime(2008, 1, 1, 0, 8, 53, 210)
    @test event.origin_time_utc == DateTime(2007, 12, 31, 15, 8, 53, 210)
    @test event.latitude ≈ 32.865
    @test event.longitude ≈ 135.30166666666668
    @test event.depth_km ≈ 33.40
    @test event.magnitude ≈ 1.8
    @test event_type(event) == :earthquake
    @test event_type_label(event) == "Earthquake"
    malformed_magnitude = line[1:52] * "AZ" * line[55:end]
    @test ismissing(parse_jma_event(malformed_magnitude, archive).magnitude)
end

@testset "May 2008 archive integration" begin
    if isdir(DATA_ROOT)
        archive = DataArchive(DATA_ROOT)
        availability = available_channels(archive)
        @test length(availability) == 33
        @test "FUJ" in available_stations(archive)
        @test all(first(item.days) == Date(2008, 5, 25) for item in availability)
        @test all(last(item.days) == Date(2008, 5, 31) for item in availability)

        start = DateTime(2008, 5, 27)
        stop = start + Second(1)
        traces = read_window(archive, start, stop; stations="FUJ", channels="wU")
        @test length(traces) == 1
        @test Seis.startdate(only(traces)) == start
        @test Seis.enddate(only(traces)) == stop
        @test Seis.nsamples(only(traces)) == 101

        recipes = suggest_filters(traces)
        @test any(recipe -> recipe.name == :broadband, recipes)
        processed = process_traces(traces; bandpass=(0.5, 10.0))
        @test length(processed) == length(traces)
        @test Seis.trace(only(processed)) != Seis.trace(only(traces))

        events = read_catalog(
            archive;
            start=DateTime(2008, 5, 25),
            stop=DateTime(2008, 5, 26),
        )
        @test !isempty(events)
        @test all(event -> DateTime(2008, 5, 25) <= event.origin_time_utc <= DateTime(2008, 5, 26), events)
        @test any(event -> event_type(event) == :lfe, events)
        @test event_type_label(first(filter(event -> event_type(event) == :lfe, events))) ==
            "Low-frequency earthquake (LFE)"
    else
        @test_skip "local Fuji archive is not available"
    end
end

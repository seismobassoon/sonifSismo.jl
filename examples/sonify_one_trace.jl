using Dates
using DSP: conv
using Seis
using SonifSismo
using Statistics
using WAV

const AUDIO_RATE = 44_100
const AUDIO_OUTPUT = joinpath(@__DIR__, "..", "audio")

"""
    resample_linear(signal, input_rate, output_rate)

Simple dependency-light interpolation for converting a sonified signal to a
standard audio sampling rate. It is adequate for exploratory sonification;
for release-quality rendering, use a band-limited resampler later.
"""
function resample_linear(signal::AbstractVector, input_rate::Real, output_rate::Real)
    input_rate > 0 || throw(ArgumentError("input_rate must be positive"))
    output_rate > 0 || throw(ArgumentError("output_rate must be positive"))
    isempty(signal) && return Float64[]
    length(signal) == 1 && return [Float64(first(signal))]

    output_length = max(1, round(Int, length(signal) * output_rate / input_rate))
    locations = range(1.0, stop=Float64(length(signal)), length=output_length)
    output = Vector{Float64}(undef, output_length)
    for (i, location) in pairs(locations)
        left = floor(Int, location)
        right = min(left + 1, length(signal))
        weight = location - left
        output[i] = (1 - weight) * signal[left] + weight * signal[right]
    end
    return output
end

function normalize_audio(signal::AbstractVector; peak::Real=0.92)
    0 < peak <= 1 || throw(ArgumentError("peak must be between 0 and 1"))
    centered = Float64.(signal) .- mean(signal)
    amplitude = maximum(abs, centered; init=0.0)
    amplitude == 0 && return centered
    return centered .* (peak / amplitude)
end

function fade_edges!(signal::AbstractVector, rate::Real; seconds::Real=0.03)
    seconds <= 0 && return signal
    nfade = min(round(Int, seconds * rate), length(signal) ÷ 2)
    nfade == 0 && return signal
    ramp = range(0.0, stop=1.0, length=nfade)
    signal[1:nfade] .*= ramp
    signal[end-nfade+1:end] .*= reverse(ramp)
    return signal
end

function add_harmonics(signal::AbstractVector; even::Real=0.0, odd::Real=0.0)
    x = normalize_audio(signal; peak=1.0)
    second_harmonic = x .^ 2 .- mean(x .^ 2)
    third_harmonic = x .^ 3
    return normalize_audio(x .+ even .* second_harmonic .+ odd .* third_harmonic)
end

"""
    sonify_trace(trace; acceleration=100, octave=0, harmonics=(0.0, 0.0))

Convert one seismic trace into mono audio. `acceleration` compresses time and
raises frequency. In this minimal resampling approach, `octave=1` doubles
that playback factor and `octave=-1` halves it; octave and duration are not
independent. `harmonics=(even, odd)` adds nonlinear timbral color without
changing the rendered duration.
"""
function sonify_trace(
    trace::Seis.AbstractTrace;
    acceleration::Real=100.0,
    octave::Real=0.0,
    harmonics::Tuple{<:Real,<:Real}=(0.0, 0.0),
    audio_rate::Integer=AUDIO_RATE,
)
    acceleration > 0 || throw(ArgumentError("acceleration must be positive"))
    playback_factor = acceleration * 2.0^octave
    seismic_rate = inv(trace.delta)
    virtual_rate = seismic_rate * playback_factor
    audio = resample_linear(Seis.trace(trace), virtual_rate, audio_rate)
    fade_edges!(audio, audio_rate)
    audio = add_harmonics(audio; even=harmonics[1], odd=harmonics[2])
    return (
        signal=audio,
        sample_rate=audio_rate,
        playback_factor=playback_factor,
        seismic_duration=length(Seis.trace(trace)) / seismic_rate,
        audio_duration=length(audio) / audio_rate,
    )
end

function write_sonification(path::AbstractString, rendered)
    mkpath(dirname(path))
    wavwrite(
        rendered.signal,
        path;
        Fs=rendered.sample_rate,
        nbits=16,
        compression=WAVE_FORMAT_PCM,
    )
    return path
end

function mono(samples::AbstractArray)
    ndims(samples) == 1 && return vec(Float64.(samples))
    return vec(mean(Float64.(samples); dims=2))
end

"""
    convolve_with_audio(rendered, music_path; kernel_seconds=2, mix=0.35)

Use the opening part of a WAV file as a convolution kernel. This is an
experimental texture operation: a true room impulse-response file produces
conventional reverberation, while a music excerpt creates a smeared hybrid.
"""
function convolve_with_audio(
    rendered,
    music_path::AbstractString;
    kernel_seconds::Real=2.0,
    mix::Real=0.35,
)
    0 <= mix <= 1 || throw(ArgumentError("mix must be between 0 and 1"))
    kernel_seconds > 0 || throw(ArgumentError("kernel_seconds must be positive"))
    source, source_rate, _, _ = wavread(music_path)
    kernel = mono(source)
    kernel = resample_linear(kernel, source_rate, rendered.sample_rate)
    kernel = kernel[1:min(length(kernel), round(Int, kernel_seconds * rendered.sample_rate))]
    fade_edges!(kernel, rendered.sample_rate; seconds=min(0.02, kernel_seconds / 4))
    kernel = normalize_audio(kernel; peak=1.0)

    wet = normalize_audio(conv(rendered.signal, kernel))
    dry = vcat(rendered.signal, zeros(length(wet) - length(rendered.signal)))
    mixed = normalize_audio((1 - mix) .* dry .+ mix .* wet)
    return merge(rendered, (signal=mixed, audio_duration=length(mixed) / rendered.sample_rate,))
end

# A reproducible MWE based on the LFE identified at 21:31 JST on 29 May 2008.
function run_lfe_mwe(; output_dir::AbstractString=AUDIO_OUTPUT)
    archive = DataArchive("/Users/nobuaki/Documents/mtFujiContinuous")
    candidate_events = read_catalog(
        archive;
        start=DateTime(2008, 5, 29, 21),
        stop=DateTime(2008, 5, 29, 23),
        input_timezone=:JST,
    )
    lfe = first(filter(event -> event_type(event) == :lfe, candidate_events))
    traces = read_window(
        archive,
        lfe.origin_time_utc - Minute(5),
        lfe.origin_time_utc + Minute(10);
        stations="FUJ",
        channels="wU",
        processed=true,
        bandpass=(0.5, 15.0),
    )
    trace = only(traces)

    plain = sonify_trace(trace; acceleration=100)
    bright = sonify_trace(trace; acceleration=100, harmonics=(0.18, 0.08))
    higher = sonify_trace(trace; acceleration=100, octave=1)

    println("Seismic window: ", round(plain.seismic_duration; digits=1), " s")
    println("Plain audio: ", round(plain.audio_duration; digits=2), " s at x", plain.playback_factor)
    write_sonification(joinpath(output_dir, "lfe_wU_x100.wav"), plain)
    write_sonification(joinpath(output_dir, "lfe_wU_x100_harmonics.wav"), bright)
    write_sonification(joinpath(output_dir, "lfe_wU_x100_octave_up.wav"), higher)
    return (; lfe, trace, plain, bright, higher)
end

(abspath(PROGRAM_FILE) == @__FILE__) && run_lfe_mwe()

# Optional use from a notebook, with a WAV file you choose:
# mwe = run_lfe_mwe()
# hybrid = convolve_with_audio(mwe.bright, "/absolute/path/to/music.wav"; kernel_seconds=1.5, mix=0.3)
# write_sonification(joinpath(AUDIO_OUTPUT, "lfe_wU_music_convolution.wav"), hybrid)

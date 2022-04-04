begin
    using GLMakie
    using FileIO
    using MusicProcessing, FLAC, MP3
    using PortAudio, SampledSignals
    using RollingFunctions
    using Statistics

    GLMakie.set_theme!(theme_black(); resolution=(2048, 1536))
    GLMakie.set_window_config!(; framerate=60)
end

function load_song(song_name)
    loaded = load(song_name)
    MusicProcessing.SampleBuf(Float64.(loaded.data), loaded.samplerate)
end

function rollingwindow(f::Function, data::Vector{T}, window::Int)::Vector{T} where {T}
    if length(data) < window
        error("ya fucked up")
    end

    window -= 1
    output::Vector{T} = similar(@view(data[begin:end-window]))
    idxes = eachindex(data)[begin:end-window]

    @Threads.threads for idx ∈ idxes
        @inbounds output[idx] = @views f(data[idx:idx+window])
    end

    output
end

begin
    audio = load_song("sounds/83_its_good_to_be_d.mp3")
    audio_mono = MusicProcessing.mono(audio)
    audio_rolling = rollingwindow(mean, Float64.(audio_mono.data), audio.samplerate |> Int)
end


function play_audio(soundtrack)
    fig, ax = lines(soundtrack.data);
    
    vline = Observable(0.0)
    vlines!(ax, vline; color=:red)
    display(fig)
    sleep(1)

    audio_stream_play = @async PortAudioStream(0, 2; samplerate=soundtrack.samplerate |> Int) do stream
        write(stream, soundtrack.data)
    end

    try
        start_time = time()
        while !istaskdone(audio_stream_play)
            vline[] = (time() - start_time) * soundtrack.samplerate
            @sync @async sleep(1/60)
        end
    catch e
        if !istaskdone(audio_stream_play)
            Base.throwto(audio_stream_play, e)
            wait(audio_stream_play)
        end
    end

    GLMakie.destroy!(GLMakie.global_gl_screen())
end


begin
    mono_rolling = copy(audio_mono)
    mono_rolling.data = audio_rolling
    #audio_mono.data .+ rand(length(audio_mono.data)) * 0.01
end

play_audio(mono_rolling)
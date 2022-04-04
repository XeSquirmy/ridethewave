begin
    using GLMakie
    using FileIO
    using MusicProcessing, FLAC, MP3
    using PortAudio
    using RollingFunctions

    GLMakie.set_theme!(theme_black(); resolution=(2048, 1536))
    GLMakie.set_window_config!(; framerate=60)
end

function something()
    loaded = load("sounds/83_its_good_to_be_d.mp3")
    MusicProcessing.SampleBuf(Float64.(loaded.data), loaded.samplerate)
end

audio = something()
audio_mono = mono(audio)
audio_rolling = rollmean(audio_mono, audio.samplerate |> Int)


scatter(audio_rolling)#cumsum(audio[:, 1], dims=1); color=1:length(audio[:, 1]))
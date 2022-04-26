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

#begin
#    audio = load_song("sounds/83_its_good_to_be_d.mp3")
#    audio_mono = MusicProcessing.mono(audio)
#    audio_rolling = rollingwindow(mean, Float64.(audio_mono.data), 60000)
#end


function play_audio(soundtrack, samplerate)
    fig, ax = lines(soundtrack);
    
    vline = Observable(0.0)
    vlines!(ax, vline; color=:red)
    display(fig)
    sleep(1)

    audio_stream_play = @async PortAudioStream(0, 2; samplerate=samplerate |> Int) do stream
        write(stream, soundtrack)
    end

    try
        start_time = time()
        while !istaskdone(audio_stream_play)
            vline[] = (time() - start_time) * samplerate
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


#begin
#    mono_rolling = copy(audio_mono)
#    mono_rolling.data = audio_rolling .* 400
#    #audio_mono.data .+ rand(length(audio_mono.data)) * 0.01
#end
#
#play_audio(mono_rolling)

begin
    using Dates
    using GLMakie
    using Primes
    using Distributions
    using Base.MathConstants
    using SpecialFunctions
    using StatsBase
    using DataFrames
    using DataFramesMeta
    using CSV
    using SHA

    GLMakie.activate!()
    GLMakie.set_theme!(
        theme_black();
        resolution = (2048, 1536),
        aspect = DataAspect())

    GLMakie.set_window_config!(; framerate=60)

    γ = MathConstants.eulergamma
    ϕ = MathConstants.golden
end

begin
    isprime_r(p::Float64) = isprime(abs(Int64(round(p))))
    isprime_c(p::Float64) = isprime(abs(Int64(ceil(p))))
    isprime_f(p::Float64) = isprime(abs(Int64(floor(p))))
    Primes.isprime(p::Float64) = isprime_f(p)

    function gen_range(iters, resolution, xs, ys)
        iters = convert(Int64, iters)
        zros = zeros(resolution, resolution)
    
        cs_real = zros .+ LinRange(xs[1], xs[2], resolution)
        cs_imag = zros .+ LinRange(ys[1], ys[2], resolution)
    
        cs::Matrix{Complex{Float64}} = cs_real .- (cs_imag .* 1im)'
        zs::Matrix{Complex{Float64}} = copy(cs)
    
        return zs, cs, zros, iters, resolution
    end

    function mandelbrotp(;
        fig = nothing,
        ax = nothing,
        iters = 10, 
        resolution = 500, 
        xs = (-4, 2),
        ys = (-4, 2),
        escape_criteria = abs2(2 + 2im))

        zs, cs, zi, iters, resolution = gen_range(iters, resolution, xs, ys)

        display(fig)
        mandelbrot_plot(zs, cs, zi, iters, resolution, escape_criteria)
    end

    function mandelbrot_plot(
        zs,
        cs,
        zi,
        iters, 
        resolution,
        escape_criteria)

        zxs = real.(cs[:])
        zys = imag.(cs[:])

        zss = Observable(zs)
        zsr = @lift(real.($zss)[:])
        zsi = @lift(imag.($zss)[:])


        zis = Observable(zi)
        ziz = @lift($zis[:])

        fig = Figure()
        display(fig)

        sleep(1)

        ax = Axis(fig[1, 1])

        scatter!(ax, zxs, zys; color=ziz, markersize=8, colormap=:inferno)
        ax.aspect = DataAspect()

        normtrunc = TruncatedNormal(0, 1.5, -0.0005, 0.0005)
        mandelbrot_iter!(
            zss,
            cs, 
            zis,
            iters, 
            escape_criteria, 
            axes(zs),
            fig,
            ax,
            (z, c, r=0) -> sum(range(;start=min(abs(z), abs(c)), stop=max(abs(z), abs(c)), length=300)),
            normtrunc,
            zxs,
            zys,
            zsr,
            zsi
        )

        zss[], cs, zis[], ziz[]
    end

    function mandelbrot_iter!(
        zs, 
        cs, 
        zi, 
        iters, 
        escape_criteria, 
        dims,
        fig,
        ax,
        formula,
        randtrunc,
        zxs,
        zys,
        zsr,
        zsi)

        known_exceptions = Bool.(zeros(axes(cs)))
        for _ ∈ 1:iters
            fractal_iter!(zs.val, cs, zi.val, dims)#, known_exceptions, escape_criteria, dims, formula, randtrunc)

            notify(zs)
            notify(zi)
        end

        notify(zs)
        notify(zi)
    end

    function fractal_iter!(zs, cs, zi, dims)#, ke, escape_criteria, dims, formula, randtrunc)
        rows, cols = dims
        #i64max = typemax(Int64)

        #rand_id = string(datetime2unix(Dates.now()))
        #rands = rand(randtrunc, (length(rows), length(cols)))
        #prevshas = [[sha256(string(round(Int, ((i + j)*(i + j + 1)/2) + j))) for j in cols] for i in rows]
        #prevshas = [[sha256(replace(string(cs[i, j]), " " => "", "im" => "")) for j in cols] for i in rows]

        #zv = zs.val
        #ziv = zi.val

        @Threads.threads for i ∈ rows
            for j ∈ cols
                @inbounds c = cs[i, j]
                @inbounds z = zs[i, j]
                @inbounds zs[i, j] = zn = fractal_calc(z, c)
                @inbounds zi[i, j] += fractal_satisfies(zn, z)# ⊻ 1
            end
        end
    end

    function fractal_calc(z::ComplexF64, c::ComplexF64)
        zeta(z^ℯ) * zeta(log(sin(z^ℯ) + cos(z^ℯ)))
    end

    function fractal_satisfies(zn::ComplexF64, zo::ComplexF64)
        #abs2(z) <= 8
        abs2(zn) <= abs2(zo)
    end
    
    #xs = (4.722, 4.73)
    #ys = (1.024, 1.032)
    xs = (-8, 8)
    ys = (-8, 8)

    zs, cs, zi, ziz = mandelbrotp(; 
        resolution = 1000,
        iters = 10,
        xs=xs,
        ys=ys);
end

#begin
#    c = -0.5514+0.626im
#    c = -0.03+0.78im
#    z = copy(c)
#    S = 1024
#    Sm = 2S
#
#    zs, cs, _, _, _ = gen_range(
#        0,
#        5,
#        (c.re - 0.006, c.re + 0.006),
#        (c.im - 0.006, c.im + 0.006))
#
#    result = [[] for i in 1:length(cs)]
#    zs = zs[:]
#    cs = cs[:]
#
#    for i in eachindex(zs)
#        for _ in 1:Sm
#            z = zs[i]
#            c = cs[i]
#
#            z = z^2 + c
#
#            zs[i] = z
#            push!(result[i], zs[i])
#        end
#    end
#
#    signal = sum([abs.(map(x -> isnan(x) || isinf(x) ? 0 : x, i)) for i in result]) / length(result)
#    play_audio(signal, S)
#end
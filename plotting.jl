begin
    using Dates
    using GLMakie
    using Primes
    using Distributions
    using Base.MathConstants
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
        #zsa = @lift(abs2.($zss)[:])
        zsr = @lift(real.($zss)[:])
        zsi = @lift(imag.($zss)[:])


        zis = Observable(zi)
        ziz = @lift($zis[:])

        fig = Figure()
        display(fig)

        sleep(1)

        #ax = LScene(fig[1, 1])
        #ax = Axis3(fig[1, 1])#; aspect=DataAspect())
        ax = Axis(fig[1, 1])#; aspect=DataAspect())

        #fig, ax, = scatter(zsr, zsi, zxs, axis=(type=LScene,); color = ziz, colormap=:inferno, markersize=100)
        #surface!(ax, zxs, zys, ziz; color = ziz, markersize=100)
        #figscatter = scatter!(ax, zxs, zys, zsr; color = zxs, colormap=:viridis, markersize=8)
        #scatter!(ax, zsr, zsi, zys; color = zxs, colormap=:viridis, markersize=100)
        #scatter!(ax, zsr, zsi, ziz; color = ziz, colormap=:viridis, markersize=15)
        #scatter!(ax, zxs, zys, ziz; color = ziz, colormap=:viridis, markersize=20)
        #scatter!(ax, zxs, zys, zsa; color=ziz, colormap=:viridis, markersize=16)
        #scatter!(ax, zxs, zys, @lift(map(x -> x > 10 ? 10 : x, $zsa)); color=ziz, colormap=:viridis, markersize=16)


        # 2d
        #scatter!(ax, zxs, zys; color=@lift(map(x -> x < 0.01 ? 5 : x > 5 ? 10 : 0, $zsa)), colormap=:viridis, markersize=8)
        #scatter!(ax, zxs, zys; color=zsa, colormap=:viridis, markersize=8)
        #scatter!(ax, zxs, zys; color=@lift(1 ./ $ziz), colormap=:inferno, markersize=8)
        #arrows!(ax, zxs, zys, zsr, zsi; arrowsize = 32, lengthscale = 0.3, arrowcolor=1:length(cs[:]), linecolor=1:length(cs[:]), colormap=:viridis)#, markersize=8)
        scatter!(ax, zxs, zys; color=ziz, markersize=8, colormap=:inferno)
        #scatter!(ax, zsr, zsi, ziz; color=ziz, markersize=8)
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
            fractal_iter!(zs, cs, zi, known_exceptions, escape_criteria, dims, formula, randtrunc)

            notify(zs)
            notify(zi)
        end
    end

    function fractal_iter!(zs, cs, zi, ke, escape_criteria, dims, formula, randtrunc)
        rows, cols = dims
        i64max = typemax(Int64)

        #rand_id = string(datetime2unix(Dates.now()))
        #rands = rand(randtrunc, (length(rows), length(cols)))
        #prevshas = [[sha256(string(round(Int, ((i + j)*(i + j + 1)/2) + j))) for j in cols] for i in rows]
        #prevshas = [[sha256(replace(string(cs[i, j]), " " => "", "im" => "")) for j in cols] for i in rows]

        @Threads.threads for i ∈ rows
            for j ∈ cols
                @inbounds z = zs.val[i, j]
                @inbounds c = cs[i, j]
                #@inbounds r = rands[i, j]
       
                #if isinf(z.re) || isinf(z.im) || isnan(z.re) || isnan(z.im) || abs2(z) > i64max
                #    continue
                #end

                #shasum = sha256(prevshas[i][j])
                zn = z^2 + c
                satisfies = abs(zn) <= 8

                zo = z
                z = zn

                zs.val[i, j] = z
                zi.val[i, j] += satisfies ⊻ 1
                #prevshas[i][j] = shasum
            end
        end
    end
    
    zs, cs, zi, ziz = mandelbrotp(; 
        resolution = 1600,
        iters = 500,
        xs=(-0.753, -0.751),
        ys=(0.037, 0.039));
end
using LinearAlgebra, Plots, Graphs, GraphRecipes, ColorSchemes, Colors, Pandoc, Random, Distributions;

"""
Aanalysis(A::Array{T}, ver::Matrix{String}; graph = false) where T <:Real

The function allows the accumulation analysis of a graph, provided the Adiacency matrix (A) and list of all vertices (ver).
"""
function Aanalysis(A, ver; graph = false) #ver vector of strings naming all vertex
    length(ver) != length(A[1,:]) && return println("Vertices vector must be consistent with Adiacency matrix")
    length(A[1,:]) != length(A[:,1]) && return println("Adiacency matrix needs to be a square one")

    function _ndv(ver) #find all type of vertex from ver list
        vecount = []
        for i ∈ ver
            if isempty(findall(x -> x == i, vecount))
                push!(vecount, i)
            end
        end
        return vecount
    end

    function _membership(ver)
       memb = zeros(Int, length(ver))
        ndv = _ndv(ver)
        for i ∈ 1:length(ndv)
            for j ∈ findall(x -> x == ndv[i], ver)
                memb[j[2]] = i
            end
        end
    return memb
    end

    function _legenda!(av, colorset)
       println("Legend color - type")
        for i ∈ 1:length(av)
            #display(colorset[i])
            println("Color: ", colorset[i],", type: ", av[i])
        end
    end

    function _accanalysis!(A)
        x₀ = zeros(Float64, length(A[1,:])) .+ 1
        x = A' * x₀
        t = findall(y -> y > 1, x)
        acc = x - x₀
        println("There is accumulation in: ", ver[t], " with net accumulation of ", acc[t], " respectively.")
    end

    function _massdistr!(A)
        V = eigvecs(A')
        λ = eigvals(A')
        π = abs.(Float64.(V[:, end])) ./ sum(abs.(Float64.(V[:, end])))
        return println("Mass distribution on vertices: π = ", π)
    end

    av = _ndv(ver) #vector containing all classes
    nav = length(av) #number of classes
    mem = _membership(ver) #type are transformed in 1, 2, 3.. (depending on number of classes) number
    colorset = ColorScheme(distinguishable_colors(nav+1))[2:end] #first color, black, is skipped
    mark = colorset[mem] #color assignment

    if graph
        _legenda!(av, colorset) #show legend
    end

    _accanalysis!(A) #find accumulation vertices in the graph

    _massdistr!(A) #find mass distribution

    if graph
        graphplot(A, #graph sketch
              markercolor = mark,
              names = ver,
              markersize = 0.2,
              nodeshape = :circle,
              nodesize = 5,
              linecolor = :darkgrey,
              linealpha = 0.5,
              )
    end

end

"""
MCanalysis(A::Array{T}, ver::Matrix{String}; iter = 1000, var = 0.1, plot = true, deb = false, mdvar = true) where T <: Real

The function allows to do a single Monte Carlo analysis of iter steps, provided the Adiacency matrix (A) and list of all vertices (ver).

"""
function MCanalysis(A::Array{T}, ver::Matrix{String}; iter::Int = 1000, var::T = 0.1, plot::Bool = true, deb::Bool = false, mdvar::Bool= true) where T <: Real
    B = copy(A)
    data = zeros(iter, length(ver))
    Ctrl = copy(A)
    Ctrl[findall(x -> x > 0, Ctrl)] .= 1 #all non-zero entries are set equal to 1

    function _chooseoneentry(Ctrl) #this choose a random non-zero entry
        c = rand(findall(x -> x > 0, Ctrl))
        return c
    end

    function _onestepmdp!(B, Ctrl, var) #MC step with mean dependent variance
        c = _chooseoneentry(Ctrl)
        d = truncated(Normal(B[c], var), 0.0, Inf) #it defines a Normal distribution with μ = B[c] and σ = sqrt(var * B[c])
        B[c] = rand(d, 1)[1] #it samples a value from that Normal and assign it to B[c]
        B[c[1], :] ./= sum(B[c[1], :]) #it normalises all the line where the change took place
        Ctrl[c] += 1 #it adds 1 in the c position of the control matrix
    end

    function _onestep!(B, Ctrl, var) #MC step with mean independent variance
        c = _chooseoneentry(Ctrl)
        d = truncated(Normal(B[c], var), 0.0, Inf) #it defines a Normal distribution with μ = B[c] and σ = sqrt(var)
        B[c] = rand(d, 1)[1] #it samples a value from that Normal and assign it to B[c]
        B[c[1], :] ./= sum(B[c[1], :]) #it normalises all the line where the change took place
        Ctrl[c] += 1 #it adds 1 in the c position of the control matrix
    end

    function _accanalysis!(A) #it performs one step of the recursive relation
        x₀ = zeros(Float64, length(A[1,:])) .+ 1
        #x₀ ./= sum(x₀)
        x = A' * x₀
        return x - x₀
    end

    function _myplot!(data, ver)
        for j ∈ 1:length(ver)
            display(plot(data[:,j], title = "Accumulation in " *  ver[j], legend = false))
            name = "Accumulation_plot_" * ver[j] * ".svg"
            savefig(name)
            deb && println("Plot done")
        end
    end

    function _MC!(B, var, data, Ctrl, iter, deb, mdvar) #this modifies data structure, line by line, introducing accumulation results for each simulation step
        if mdvar
            for i ∈ 1:iter
                _onestepmdp!(B, Ctrl, var)
                data[i, :] = _accanalysis!(B)
                deb && println("Internal cycle _MC mdv, step:", i)
            end
        else
            for i ∈ 1:iter
                _onestep!(B, Ctrl, var)
                data[i, :] = _accanalysis!(B)
                deb && println("Internal cycle _MC muv, step:", i)
            end
        end
    end

    function _ctrlpost(Ctrl)
        c = findall(x -> x > 0, Ctrl)
        Ctrl[c] = Ctrl[c] .- 1
        return Ctrl
    end

    deb && println("Declaration done")

    _MC!(B, var, data, Ctrl, iter, deb, mdvar)

    deb && println("Data acquired")

    if plot
        _myplot!(data, ver)
    end

    return data, Int.(_ctrlpost(Ctrl))
end


"""
MC(A::Array{T}, ver::Matrix{String}; itmc::Int=10000, itr::Int = 1000, var::T = 0.1, deb::Bool = false) where T <:Real

The function allows to do itmc Monte Carlo analysis each one of iter steps, provided the Adiacency matrix (A) and list of all vertices (ver).
Final result is then averaged on all realisations of each MC analysis.
"""
function MC(A::Array{T}, ver::Matrix{String}; itmc::Int=10000, itr::Int = 1000, var::T = 0.1, deb::Bool = false, mdvar::Bool = true, savef::Bool = true) where T <:Real
    data = []
    μ = zeros(itr, length(ver))
    σ = zeros(itr, length(ver))
    B = copy(A)
    #Bm = copy(A)
    Ctrl = copy(A)
    Ctrl[findall(x -> x > 0, Ctrl)] .= 1

    function _myplot!(μ, σ, ver, itr, savef, mdvar)
        thr = zeros(itr)
        for j ∈ 1:length(ver)
            pl = plot(μ[:,j], title = "Net accumulation measure in " *  ver[j], legend = false, lw = 2, label = "μ", xlabel = "Time steps", ylabel = "Net accumulation measure", linecolor = :black)
            up = μ + σ
            down = μ - σ
            plot!(up[:, j], label = "μ + σ", lw = 1.5, linecolor = :black, linestyle = :dot)
            plot!(down[:, j], label = "μ - σ", lw = 1.5, linecolor = :black, linestyle = :dot)
            plot!(thr, label = "Acc. threshold", lw = 1.5, linestyle = :dash, linecolor = :black)
            display(pl)
            if mdvar
                name = "hm_mdv_" *  ver[j] * ".svg"
            else
                name = "hm_nmdv_" *  ver[j] * ".svg"
            end
            savef && savefig(name)
            deb && println("Plot done")
        end
    end

    function _ctrlpost(Ctrl)
        c = findall(x -> x > 0, Ctrl)
        Ctrl[c] = Ctrl[c] .- 1
        return Ctrl
    end

    deb && println("Declaration done")

    for i ∈ 1 : itmc
        D = MCanalysis(B, ver, iter = itr, var = 0.1, plot = false, deb = deb, mdvar = mdvar)
        push!(data, D[1])
        Ctrl += D[2]
        #Bm += D[3]
        deb && println(i, " simulation done")
    end

    for i ∈ 1:itmc
        μ += data[i]
    end

    μ ./= itmc

    for i ∈ 1:itmc
       σ += (((data[i] - μ) .^ 2) ./ (itmc - 1))
    end

    std = sqrt.(σ)

    nCt = _ctrlpost(Ctrl) ./ itmc
    #Bm ./= itmc

    deb && println("Data and Ctrl normalised")

    _myplot!(μ, std, ver, itr, savef, mdvar)

    return nCt
end

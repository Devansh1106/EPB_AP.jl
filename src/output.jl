# CSV output (SPEC: one header block of "# key = value" lines, then a column-name row, then data).

"""
    write_csv(path, meta, names, columns)

Write `columns` (equal-length vectors) under `names`, preceded by `# key = value` lines from `meta`.
"""
function write_csv(path, meta, names, columns)
    mkpath(dirname(path))
    open(path, "w") do io
        for (k, v) in meta
            println(io, "# ", k, " = ", v)
        end
        println(io, join(names, ","))
        for i in eachindex(columns[1])
            println(io, join((@sprintf("%.16e", c[i]) for c in columns), ","))
        end
    end
    return path
end

"""
    read_csv(path) -> (meta::Dict{String,String}, names::Vector{String}, data::Matrix{Float64})
"""
function read_csv(path)
    meta, names, rows = Dict{String,String}(), String[], Vector{Float64}[]
    for line in eachline(path)
        if startswith(line, "#")
            k, v = split(line[2:end], "="; limit = 2)
            meta[strip(k)] = strip(v)
        elseif isempty(names)
            names = String.(split(line, ","))
        elseif !isempty(line)
            push!(rows, parse.(Float64, split(line, ",")))
        end
    end
    return meta, names, reduce(vcat, permutedims.(rows); init = zeros(0, length(names)))
end

function case_meta(case::TestCase, res::Result, t)
    s = case.scheme
    m = Pair{String,Any}["case" => case.name, "scheme" => s.kind, "eta" => s.eta]
    s.kind == :dissipation &&
        append!(m, ["theta" => s.theta, "kappa" => "$(s.kappa_c)*dx*(max|u|+$(s.kappa_s))"])
    return [m; Pair{String,Any}["N" => case.N, "domain" => "[$(case.domain[1]), $(case.domain[2])]", "lambda" => case.lambda,
        "T" => case.T, "t" => t, "steps" => res.steps, "status" => res.status,
        "timestep" => "mode=$(case.timestep.mode), cfl=$(case.timestep.cfl), dtmax=$(case.timestep.dtmax)"]]
end

function solution_columns(case, x, t, ρ, u, φ)
    names, cols = ["x", "rho", "u", "phi"], Any[x, ρ, u, φ]
    if case.exact !== nothing
        ρe, ue, φe = case.exact(x, t)
        append!(names, ["rho_exact", "u_exact", "phi_exact"]); append!(cols, [ρe, ue, φe])
    end
    return names, cols
end

"""
    run_case(name_or_file; dir = "data", kwargs...) -> (result, paths)

Solve a test case, given by name (`params/<name>.toml`) or by parameter file, with keywords as in
[`testcase`](@ref), and write CSV files to `dir`:
the solution at the final time, one file per snapshot, and the per-step history.
"""
function run_case(src::Union{Symbol,AbstractString}; dir = "data", kw...)
    case = testcase(src; kw...)
    res = solve(case)
    base = joinpath(dir, "$(case.name)_$(case.scheme.kind)_N$(case.N)_lambda$(case.lambda)_T$(case.T)")
    paths = String[]
    for (t, ρ, u, φ) in res.snapshots
        meta = case_meta(case, res, t)
        push!(paths, write_csv("$(base)_t$(t).csv", meta, solution_columns(case, res.grid.x, t, ρ, u, φ)...))
    end
    meta = [case_meta(case, res, res.t); case.diagnostics(case, res)]
    push!(paths, write_csv("$(base).csv", meta,
                           solution_columns(case, res.grid.x, res.t, res.rho, res.u, res.phi)...))
    h = res.history
    push!(paths, write_csv("$(base)_history.csv", meta, collect(String.(keys(h))), collect(values(h))))
    return res, paths
end

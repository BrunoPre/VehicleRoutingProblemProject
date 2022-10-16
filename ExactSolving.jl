#= Operational Research Course Project 2021
    Exact solving
    PREKA Bruno
    ZELLE Yannick
=#

using JuMP, GLPK
using TravelingSalesmanExact
include("DataVRP.jl")
using .DataVRP

### Modeling
function build_model_exact(selectedSolver::DataType, AllS_i::Vector{Vector{Int64}}, nbClients::Int64, resTSP::Vector{Tuple{Vector{Int64},Int64}}, nbSubsets::Int64)
    #= Glossary:
        * AllS_i:    for all i, AllS_i[i] is the vector containing the routes in which client `i` is visited   
        * resTSP:    contains all the tuples (route, min_distance) given by TSP solving
        * nbSubsets: number of possible subsets of clients
    =#

    # Init the model
    m::Model = Model(selectedSolver)

    # Decision variables
    # if the j-th route is picked, then x[j]=1. Else, x[j]=0 
    @variable(m, x[1:nbSubsets]>=0, binary = true)

    # Objective function to be optimised
    # `resTSP` is a vector of tuples (route, distance), hence only distance is kept, i.e.: resTSP[j][2]
    @objective(m, Min, sum(resTSP[j][2]*x[j] for j in 1:nbSubsets))

    # Constraint: every client is visited only once
    @constraint(m, VisitOnlyOnceClient[i=2:nbClients], sum(x[j] for j in AllS_i[i-1]) == 1)

    return m
end

#= Computes the possible subsets of clients ("constrained" by the capacity of the delivery vehicle),
    given the following:
    * @param P: set of the pending clients to add
    * @param S: subsets of clients
    * @param capacity: the drone capacity
    * @param demand: demand per client
    * @param index: current client
    * @param d: demand accumulator that has to be less or equal than `capacity`
    @returns the set of all the possible subsets of clients
=#
function getSubsets_recursive(P::Vector{Int64}, S::Vector{Vector{Int64}}, capacity::Int64, demand::Vector{Int64}, clientIndex::Int64, d::Int64)

    if clientIndex < 0
        error("clientIndex is out of bounds")
    end

    # Stop condition: every client has been treated
    if clientIndex > length(demand)
        return S                        
    end
    
    while clientIndex <= length(demand) 
        if d + demand[clientIndex] <= capacity    # constraint
            toAdd::Vector{Int64} = copy(P)               # shallow copy of the pending clients to add
            toAdd = append!(toAdd,clientIndex+1)         # with current `clientIndex`, constraint still holds ==> add next client to the pending to-add set
            S = vcat(S,[toAdd])                          # `toAdd` fits the constraint ==> keep it

            Snew::Vector{Vector{Int64}} = getSubsets_recursive(toAdd, S, capacity, demand, clientIndex + 1, demand[clientIndex] + d)

            # add new subsets to S
            if length(S) < length(Snew)
                S = vcat(S, Snew[length(S) + 1 : length(Snew)])        
            end
        end
        clientIndex += 1
    end
    return S
end

### Wrapper function of getSubsets_recursive
function getSubsets(capacity::Int64,demand::Vector{Int64})
    P::Vector{Int64} = []
    S::Vector{Vector{Int64}} = []
    return getSubsets_recursive(P,S,capacity,demand,1,0)
end

### Compute the index of each subset/route where client `cli` is visited
## @returns vector of indexes in S
function getSetofCyclesClient(S::Vector{Vector{Int64}}, cli::Int64)
    res::Vector{Int64} = []
    for i in 1:length(S)
        if cli in S[i]
            res = push!(res,i)
        end
    end
    return res
end

### Computes the shortest cycle/route in a given set of to-be-visited clients
## @param Si: set of clients (as indexes)
## @param d: distances between the clients
function determineShortestCycle(Si::Vector{Int64}, d::Matrix{Int64})
    Si = append!([1],Si) # add depot for TSP solving
    _d::Matrix{Int64} = d[Si,Si] # restrict the matrix of distances to the targeted clients

    _cycle::Vector{Int64} = []
    distMin::Int64 = 0

    # TSP solving
    _cycle, distMin = solveTSPExact(_d)

    # re-arrange the indexes
    cycle_new::Vector{Int64} = []
    k::Int64 = 0
    
    for i in eachindex(_cycle)
        k = _cycle[i] # index provided by TSP is the client (as index) in `Si`
        cycle_new = push!(cycle_new,Si[k]) # "new" cycle has the right indexes
    end

    return cycle_new, distMin
end

#= Returns the tuples (route, min_dist) for each subset of S
    * @param S: subsets of clients
    * @param d: distances for TSP solving
    @returns (route, min_dist) tuples resulting from TSP solving for each subset of clients
=#
function getAllShortestCycles(S::Vector{Vector{Int64}}, d::Matrix{Int64})
    res::Vector{Tuple{Vector{Int64},Int64}} = []
    tmin::Int64 = length(S[1])
    tmax::Int64 = length(S[1])
    for i in eachindex(S)
        if (length(S[i]) < tmin) tmin = length(S[i]) end
        if (length(S[i]) > tmax) tmax = length(S[i]) end

        res = push!(res, determineShortestCycle(S[i], d))
    end
    println("       (Length (number of nodes) of the shortest cycle (without the depot!) : ", tmin, ")")
    println("       (Length (number of nodes) of the longest cycle (without the depot!) : ", tmax, ")")
    return res
end

#= Solves the Traveling Salesman Problem (aka TSP)
    @param d: Matrix of distances between each client to visit
    @returns a tuple `(cycles, distance)`, with `cycles` as the sequence of clients to visit (including the return trip to the depot)
=#
function solveTSPExact(d::Matrix{Int64})
    cycle::Vector{Int64} = []
    totalDist::Int64 = 0
    nbPlacesToVisit::Int64 = size(d,1)

    # Trivial case ==> TSP solving function is no longer needed
    if nbPlacesToVisit <= 3
        cycle = [i for i in 1:nbPlacesToVisit]
        totalDist = d[nbPlacesToVisit, 1]
        for i in 1:nbPlacesToVisit-1
            totalDist += d[i,i+1]
        end
        
    # General case
    else
        cycle, totalDist = TravelingSalesmanExact.get_optimal_tour(d, GLPK.Optimizer)
    end

    return cycle, round(Int64,totalDist)
end

### Reads and parse data prior to solving
function data_then_solve_exact(filename::String)

    data::DataVRP.DataObj = DataVRP.read_data(filename)
    nbClients::Int64 = data.nbClients
    distances::Matrix{Int64} = data.distance
    capa::Int64 = data.capacity
    demand::Vector{Int64} = data.demand

    # Compute all the possible subsets of clients, satisfying the constraint of vehicle capacity
    println("CPU time for building the instance of set partitionning:")
    S::Vector{Vector{Int64}} = @time getSubsets(capa,demand)
    println()
    println("Total of possible subsets: ", length(S))
    println()

    # Tuples (cycle, dist_min)
    println("CPU time to get the shortest cycles (cycle, dist_min) for each subset using TSP:")
    l::Vector{Tuple{Vector{Int64},Int64}} = @time getAllShortestCycles(S,distances)
    println()

    # For each client `i`, AllS_i[i] contains the indexes to the cycles visting client `i`
    AllS_i::Vector{Vector{Int64}} = []
    for i in 2:nbClients
        AllS_i = push!(AllS_i, getSetofCyclesClient(S,i))
    end

    # Number of possible subsets of clients
    nbSubsets::Int64 = length(S)

    # Building the model from the given data
    m::Model = build_model_exact(GLPK.Optimizer, AllS_i, nbClients, l, nbSubsets)

    # SOLVING
    println("CPU time for the linear program solving: ")
    @time optimize!(m)
    println()

    # Print results (handling solver status code)
    status::MOI.TerminationStatusCode = termination_status(m)

    if status == MOI.OPTIMAL
        println("Problem exactly solved with optimal solution")

        # Print the optimal cycles
        println("Cycles found: ")
        for j in 1:nbSubsets
            if (value(m[:x][j]) == 1.0)
                println("no. ", j, ": ", l[j][1], " with distance ", l[j][2])
            end
        end

        # Print the optimal distance
        println("Minimal total distance: z = ", round(Int64, objective_value(m)))

        println()

    elseif status == MOI.INFEASIBLE
        println("Unbounded problem")

    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Infeasible problem")
    end

    # Last print before that of @time used in `timer_res_exact` function
    println("Total CPU time:")

end

### Main/wrapper function to run from REPL
function timer_res_exact(filename::String)
    @time data_then_solve_exact(filename)
end




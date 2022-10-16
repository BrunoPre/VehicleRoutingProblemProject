#= Operational Research Course Project 2021
    Solving approach based on Clark and Wright's algorithm
    PREKA Bruno
    ZELLE Yannick
=#

using JuMP, GLPK
using TravelingSalesmanExact
include("DataVRP.jl")
using .DataVRP


#= Computes the savings for the Clark and Wright's algorithm, given the distances between the clients
    @param d: distances between the clients
    @returns the tuples `((i,j), s_ij)`, such that `s_ij` is the distance saving resulting from merging clients `i` and `j` into a single tour
=#
function getSavings(d::Matrix{Int64})
    #= Matrix `d` is symmetric, hence the upper triangular part is iterated (diagonal is zero-ed).
        Savings are only computed between clients, not the depot. That is why iteration starts from 2nd row
        Matrix looks like:
        [ 0 a12 a13 ... a1m  <-- excluded because we're not interested in savings between any i>1 and depot (indexed 1)
             0  a23 ... a2m  <-- iteration starts at row i=2, column j=3
                 0  ... a3m  <-- an offset is required to shift `j` index 
          ................
                         0 ]

    =#
    n,m::Int64 = size(d)
    savings::Vector{Tuple{Tuple{Int64,Int64},Int64}} = []
    offset::Int64 = 3        # offset from the zero-ed diagonal
    s_ij::Int64 = 0
    
    for i in 2:n    # rows
        for j in offset:m   # columns after the zero-ed diagonal
            # compute the associated saving
            s_ij = d[i,1] + d[1,j] - d[i,j]    
            savings = push!(savings, ((i,j), s_ij))
        end
        offset += 1 # adjust the offset for the next line 
    end
    return savings
end

#= Merges elementary cycles/tours, as part of the Clark and Wright's algorithm.
    * @param ordered_savings: all the eligible candidates ((i,j), s_ij) for merger, 
                            such that s_ij is the saving distance between i and j,
                            and descending-sorted by savings s_ij
    * @param init_sol: initial tuples (tour, 0),
                    such that `tour` is an elementary cycle, i.e. a vector of the form [1,i,1].
    * @param demand: client demands
    * @param capacity: vehicle capacity
    @returns admissible solution tuples (tour, idx_merged_tour), where `tour` is the result of a merger, 
                                                            and `idx_merged_tour` is the index of the cycle that has been merged with `tour`.
                                                            If `idx_merged_tour == 0`, then `tour` was not merged with any other cycle.
=# 
function getAllMergedCycles(ordered_savings::Vector{Tuple{Tuple{Int64,Int64},Int64}}, init_sol::Vector{Tuple{Vector{Int64},Int64}}, demand::Vector{Int64}, capacity::Int64)
    sol::Vector{Tuple{Vector{Int64},Int64}} = init_sol
    iterator_sol::Int64 = 1

    for saving::Tuple{Tuple{Int64,Int64},Int64} in ordered_savings
        firstClient::Int64 = saving[1][1]
        secondClient::Int64 = saving[1][2]

        # find the index of the merged tour where the client has to be considered (if it's not the same index)
        indexFirstClient::Int64 = getIndex(firstClient-1, sol)   
        indexSecondClient::Int64 = getIndex(secondClient-1, sol)

        # tours containing respectively 1st & 2nd client
        firstClientTour::Vector{Int64} = sol[indexFirstClient][1]
        secondClientTour::Vector{Int64} = sol[indexSecondClient][1]

        # merging happens if the capacity constraint holds
        if(demand[indexFirstClient] + demand[indexSecondClient] <= capacity) 
            
            merger::Vector{Int64} = mergeTours(firstClientTour, secondClientTour, firstClient, secondClient)

            # firstClientTour is elementary [1,firstClient,1], and merger[2] includes `firstClient` ==> solution is kept
            if(length(firstClientTour) == 3)
                firstClient=-1
            end

            # same for `secondClient`
            if(length(secondClientTour) ==3)
                secondClient=-1
            end

            # if the merger becomes more complex and of the form [1,a,...i,j,...,b,1],
            # client candidate tuples (a,b), (_,i), (i,_), (j,_), (_,j) are removed from the ordered_savings since they have been used in the merging process and are no longer necessary 
            if (length(merger) > 4 && iterator_sol+1 != length(ordered_savings))
                ordered_savings = removeRedundantCandidates(firstClient,     # i
                                                            secondClient,    # j
                                                            (merger[2], merger[length(merger)-1]), # (a,b)
                                                            ordered_savings, 
                                                            iterator_sol + 1)
            end

            # overwrite at indexFirstClient; fresh merger is ready to be merged again ==> set its `idx_merged_tour` to 0 in the set of merged paths `sol`
            sol[indexFirstClient] = (merger, 0)
            demand[indexFirstClient] += demand[indexSecondClient] # do not forget to sum the demands up

            # 2nd client having been used for merger, indexSecondClient can be emptied and point to the index of the merger
            sol[indexSecondClient] = ([], indexFirstClient)
            
        end

        iterator_sol += 1
        
    end
    return sol
end

#= Finds the index of the merged tour where a client is visited
    * @param client
    * @param sol: tuples (tour, idx_merged_tour),  where `tour` is the result of a merger, 
                                                and `idx_merged_tour` is the index of the cycle that has been merged with `tour`.
                                                If `idx_merged_tour == 0`, then `tour` was not merged with any other cycle.
    @returns the expected index
=#
function getIndex(client::Int64, sol::Vector{Tuple{Vector{Int64},Int64}})
    # try: the client has not moved from its initial elementary cycle 1 -> client -> 1, hence the tour is at the same index
    idx_merged_tour::Int64 = sol[client][2] 
    if (idx_merged_tour == 0)
        return client
    end

    # else: the tour was merged with another one, hence the client has to be considered at the merger
    return getIndex(idx_merged_tour, sol)
end

#= Merges two given paths/tours
    * @params ai, bj: two paths
    * @params i, j: two clients contained respectively in ai and bj
=#
function mergeTours(ai::Vector{Int64}, bj::Vector{Int64}, i::Int64, j::Int64)   
    dist_ai::Int64 = length(ai)
    dist_bj::Int64 = length(bj)

    ## in the code comments, concatenation operation is abbreviated "ai+bj"    

    ## trivial case: ai is an elementary cycle, i.e.: ai = [1,i,1]
    if (dist_ai == 3)
        # if bj = [1,j,...,1]
        if (bj[2] == j)
            return append!(ai[1:dist_ai-1], bj[2:dist_bj])  # ai+bj = [1,i,j,...,1]
        end
        # else, i.e. bj = [1,...,j,1]
        return append!(bj[1:dist_bj-1], ai[2:dist_ai]) # bj+ai = [1,...,j,i,1]
    end

    ## ai = [1,...,i,1]
    if (ai[2] != i)
        # bj = [1,j,...,1]
        if (bj[2] == j)                                                     
            return append!(ai[1:dist_ai-1], bj[2:dist_bj]) # ai+bj = [1,...,i,j,...,1]
        end
        # bj = [1,...,j,1]
        bj=reverse(bj)  # force bj = [1,j,...,1]                                               
        return append!(ai[1:dist_ai-1], bj[2:dist_bj]) # ai+reverse(bj) = [1,...,i,j,...,1]
    end    

    ## ai = [1,i,...,1] and bj = [1,...,j,1]
    if (bj[dist_bj-1] == j)
        return append!(bj[1:dist_bj-1], ai[2:dist_ai]) # bj+ai = [1,...,j,i,...,1]
    end

    ## ai = [1,i,...,1] and bj = [1,j,...,1]
    ai = reverse(ai)                              
    return append!(ai[1:dist_ai-1], bj[2:dist_bj]) # reverse(ai)+bj= [1,...,i,j,...,1]
end

#= Removes redundant client candidates according to the given rules.
    * @params removeA, removeB: clients s.t. there are paths [1,..., removeA, removeB, ..., 1] --> removeA, removeB have to be removed from the set of the saving candidates ((i,j), s_ij)
                If their value is 0, then they're not intended to be removed.
    * @param removeC: client tuple (a,b) s.t. there are paths [1,a,...,b,1] --> a,b have to be removed as well
    * @param ordered_savings: set of the saving candidates ((i,j), s_ij)
    * @param start_index: iterator from where to start looping through `ordered_savings` (looking before `start_index` is useless because prior solutions are already treated)
    @returns the altered set of the saving candidates
=#
function removeRedundantCandidates(removeA::Int64, removeB::Int64, removeC::Tuple{Int64,Int64}, ordered_savings::Vector{Tuple{Tuple{Int64,Int64},Int64}}, start_index::Int64)

    # set an order in the tuple 
    if (removeC[1] > removeC[2])
        removeC = (removeC[2], removeC[1])         
    end

    _idx::Int64 = start_index
                                                
    while (_idx < length(ordered_savings))
        _client12::Tuple{Int64,Int64} = ordered_savings[_idx][1]
        client1::Int64 = _client12[1]
        client2::Int64 = _client12[2]
        
        if (client1 == removeA || client2 == removeB || (client1, client2) == removeC)           
            ordered_savings = deleteat!(ordered_savings, _idx)
        else
            _idx += 1
        end
    end

    return ordered_savings
end

## Sets the tuples (tour, idx_merged_tour), s.t. `tour` is elementary and idx_merged_tour=0 (because there's no merger initially)
function getElementaryTours(nbClients::Int64)
    return [([1,client,1], 0) for client in 2:nbClients]
end

## Removes unused cycles (either unmerged or merged), i.e. their idx_merged_tour != 0, as well as getting rid of the `idx_merged_tour` attribute
function removeUnusedCyclesAndFormat(allCycles::Vector{Tuple{Vector{Int64},Int64}})
    return [cycle[1] for cycle in allCycles if cycle[2]==0] # Vector{Vector{Int64}}
end

## Sorts the tuples ((i,j), s_ij) in descending order of s_ij
function sortVector(toSort::Vector{Tuple{Tuple{Int64,Int64},Int64}})
    return sort!(toSort, by = x::Tuple{Tuple{Int64,Int64},Int64} -> x[2], rev=true)
end

## Computes the total distance of a cycle/tour, given a matrix of distances
function getCycleLength(cycle::Vector{Int64}, d::Matrix{Int64})
    sum::Int64 = 0
    i::Int64 = 1
    while i+1 <= length(cycle)
        sum += d[cycle[i], cycle[i+1]]
        i += 1
    end
    return sum
end

## Main function to call in the REPL; wraps the data initialisation and solving

function data_then_solve_appr(filename::String)
    # Read and parse data from file, then init to data structure
    data::DataVRP.DataObj = DataVRP.read_data(filename)
    nbClients::Int64 = data.nbClients
    distance::Matrix{Int64} = data.distance
    capa::Int64 = data.capacity
    dmd::Vector{Int64} = data.demand

    savings::Vector{Tuple{Tuple{Int64,Int64},Int64}} = getSavings(distance)

    savings = sortVector(savings)
    
    elementaryCycles::Vector{Tuple{Vector{Int64},Int64}} = getElementaryTours(nbClients)

    fusionnedCycles::Vector{Tuple{Vector{Int64},Int64}} = getAllMergedCycles(savings, elementaryCycles, dmd, capa)

    formattedCycles::Vector{Vector{Int64}} = removeUnusedCyclesAndFormat(fusionnedCycles)

    # Print every cycle & total distance
    sumtot::Int64 = 0
    distcycle::Int64 = 0

    println("Solution tours obtained by Clark & Wright's algorithm are:")

    for cycle::Vector{Int64} in formattedCycles
        distcycle = getCycleLength(cycle, distance)
        sumtot += distcycle
        println("* ", cycle, ", distance = ", distcycle)
    end

    println("Total distance: ", sumtot)
end
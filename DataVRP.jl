module DataVRP

mutable struct DataObj
    nbClients::Int64 # Number of clients (including the depot)
    capacity::Int64 # Fixed capacity of the delivery vehicle 
    demand::Vector{Int64} # Demand by client
    distance::Matrix{Int64} # distance[i,j] is the distance between clients i and j (Matrix{Int64} is equivalent to Array{Int64,2})
end

### Reads an instance's data
function read_data(filename::String)
    # Open the file in read-only mode
    f::IOStream = open(filename,"r")

    # String storing each line before parsing
    s::String = ""

    # First line gives the number of clients
    s = readline(f)
    nbClients::Int64 = parse(Int64,s)

    # Second line gives the capacity of the delivery vehicle
	s = readline(f)
	capacity::Int64 = parse(Int64,s);

	# Get demands
    s = readline(f)
    demand::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

	# Get distances and store them in a square matrix
    distance::Matrix{Int64} = Matrix{Int64}(undef,nbClients,nbClients)
    for i in 1:nbClients
        s = readline(f)
        distance[i,:] = parse.(Int64,split(s," ",keepempty = false))
    end

	close(f)

    return DataObj(nbClients, capacity, demand, distance)
end

end
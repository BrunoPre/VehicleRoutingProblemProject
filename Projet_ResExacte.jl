#= PROJET RO 2021
    Résolution exacte
    PREKA Bruno
    ZELLE Yannick
    681C
=#

using JuMP, GLPK
using TravelingSalesmanExact # pour des fonctions du problème du voyageur de commerce

# Structure contenant les données du problème
mutable struct donnees
    nbClients::Int64 # Nombre de clients (y compris le dépôt)
    capacite::Int64 # Capacité du véhicule de livraison
    demande::Vector{Int64} # Demande de chaque client
    distance::Matrix{Int64} # Distancier (Matrix{Int64} est équivalent à Array{Int64,2})
end

# Fonction de lecture des données du problème
function lecture_donnees(nom_fichier::String)
    # Ouverture d'un fichier en lecture
    f::IOStream = open(nom_fichier,"r")

    # Lecture de la première ligne pour connaître le nombre de clients
    s::String = readline(f) # lecture d'une ligne et stockage dans une chaîne de caractères
    nbClients::Int64 = parse(Int64,s) # Transformation de la chaîne de caractère en entier

    # Lecture de la deuxième ligne pour connaître la capacité des véhicules
	s = readline(f)
	capacite::Int64 = parse(Int64,s);

	# Lecture des demandes
    s = readline(f)
    demande::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

	# Lecture du distancier
    distance::Matrix{Int64} = Matrix{Int64}(undef,nbClients,nbClients) # Allocation mémoire pour une matrice de taille nbClients * nbClients
    for i in 1:nbClients
        s = readline(f)
        distance[i,:] = parse.(Int64,split(s," ",keepempty = false))
    end

    # Fermeture du fichier
	close(f)

    # Retour
    return donnees(nbClients,capacite,demande,distance)
end

# Fonction de modèle
function model_exact(solverSelected::DataType, AllS_i::Vector{Vector{Int64}}, nbClient::Int64, l::Vector{Tuple{Vector{Int64},Int64}}, nbRegroup::Int64)
    # Vocabulaire :
    #= AllS_i contient les ensembles des indices des tournées dans lesquels un client i est desservi,
        i.e. pour tout i, AllS_i[i] est le vecteur des indices des tournées dans lesquels le client i est visité =#
    # l est le vecteur de tous les couples (tournée,distance_min) déterminés par TSP
    # nbRegroup est le nombre de regroupements 

    # Déclaration d'un modèle (initialement vide)
    m::Model = Model(solverSelected)

    # Déclaration des variables de décision
    # si on choisit la tournée indicée j dans l'ensemble \mathcal(S), alors x[j]=1. Sinon x[j]=0
    @variable(m, x[1:nbRegroup]>=0, binary = true)

    # Déclaration de la fonction objectif (avec le sens d'optimisation)
    # "l" étant un vecteur de couples (tournée, distance), on ne s'intéresse qu'à la distance, soit : l[j][2]
    @objective(m, Min, sum(l[j][2]*x[j] for j in 1:nbRegroup))

    # Déclaration de la contrainte garantissant que chaque client soit visité une seule fois
    @constraint(m, VisitOnlyOnceClient[i=2:nbClient], sum(x[j] for j in AllS_i[i-1]) == 1)

    return m
end

#= Fonction déterminant l'ensemble des regroupements, étant donnés :
* deux ensembles initialement vides,
* la capacité du drone, 
* un vecteur contenant la demande des clients,
* deux entiers :
    * index est un curseur initialisé à 0 au premier appel
    * d est un accumulateur de la demande (cf condition sur \mathcal(S) : "somme des d_j <= Ca")
La fonction retourne alors l'ensemble S qui contient tous les regroupements possibles =#
function getSubsets_recursive(P::Vector{Int64}, S::Vector{Vector{Int64}},capacite::Int64, demande::Vector{Int64}, index::Int64, d::Int64, distances::Matrix{Int64})
    # Condition d'arrêt : le curseur "index" a passé en revue toutes les demandes / tous les clients
    if index > length(demande)
        return S                        
    end
    # Itération sur toutes les demandes des clients
    while index <= length(demande) 
        if d+demande[index]<=capacite   # condition de \mathcal(S)                                                                                       
            toadd::Vector{Int64} = copy(P)                   # copie profonde de P pour que les changements effectués dans toadd ne soit pas effectués dans P
            toadd = append!(toadd,index+1)                     # d'après la conditionnelle, la demande des totale des clients dans P, sommée avec demande[index], est inférieure à la capacité. Donc on construit l'union des candidats dans P et du candidat "index"
            S = vcat(S,[toadd])                                 # On concatène cette union à S
            Snew::Vector{Vector{Int64}}=getSubsets_recursive(toadd,S,capacite,demande,index+1,demande[index]+d,distances)  #appel récursif sur cette union, la nouvelle demande et la demande de tous les clients dans "toadd". Le curseur "index" est incrémenté parce qu'on veut pas considérer le même client
            if length(S)<length(Snew)                         #On veut seulement ajouter le fin du résultat de l'appel récursif si les deux vecteurs sont différents
                S=vcat(S,Snew[length(S)+1:length(Snew)])        #dans ce cas, on ne concatène que la partie nouvelle du vecteur résultant de l'appel récursif
            end
        end
        index+=1
    end
    return S
end

# fonction d'encapsulation de la fonction récursive getSubsets_recursive
function getSubsets(capacite::Int64,demande::Vector{Int64}, distances::Matrix{Int64})
    P::Vector{Int64}=[]
    S::Vector{Vector{Int64}}=[]
    return getSubsets_recursive(P,S,capacite,demande,1,0,distances)
end

# fonction déterminant l'ensemble de numéros de regroupements dans lesquels le client "cli" est livré
function getSetofCyclesClient(S::Vector{Vector{Int64}}, cli::Int64)
    res::Vector{Int64} = []
    for i in 1:length(S) # pour chaque regroupement...
        if cli in S[i]  # ...si le client "cli" y appartient...
            res = push!(res,i)    # ...alors on ajoute l'indice de ce regroupement dans notre vecteur final
        end
    end
    return res
end


# fonction déterminant la tournée et sa distance la plus courte dans un regroupement "Si", étant donné un distancier "d"
function determineShortestCycle(Si::Vector{Int64}, d::Matrix{Int64}) # Si est l'ensemble des clients à visiter et d est le distancier
    Si = append!([1],Si) # on n'oublie pas d'inclure le dépôt pour TSP
    newd::Matrix{Int64} = d[Si,Si] # restriction du distancier sur les clients à visiter et le dépôt

    # initialisations
    tournee::Vector{Int64} = []
    distmin::Int64 = 0

    # résolution du problème TSP
    tournee, distmin = solveTSPExact(newd)

    # réajustement des indices
    newtournee::Vector{Int64} = []
    k::Int64 = 0 # initialisation
    for i in 1:length(tournee)
        k = tournee[i] # l'indice donné par TSP est le numéro du client de Si
        newtournee = push!(newtournee,Si[k]) # insertion du numéro correspondant dans la "nouvelle" tournée
    end

    return newtournee, distmin
end

# fonction fournissant le vecteur des couples (tournée,distance_min) de chaque regroupement de S
function getAllShortestCycles(S::Vector{Vector{Int64}}, d::Matrix{Int64})
    res::Vector{Tuple{Vector{Int64},Int64}} = [] # initialisation
    tmin::Int64 = length(S[1])
    tmax::Int64 = length(S[1])
    for i in 1:length(S)    # pour chaque regroupement...
        # déterminer les tailles min et max
        if (length(S[i]) < tmin) tmin = length(S[i]) end
        if (length(S[i]) > tmax) tmax = length(S[i]) end

        res = push!(res,determineShortestCycle(S[i],d)) # ...ajouter son couple correspondant au vecteur à retourner
    end
    println("       (Taille du plus petit regroupement (sans compter le dépôt !) : ", tmin, ")")
    println("       (Taille du plus grand regroupement (sans compter le dépôt !) : ", tmax, ")")
    return res
end

# Fonction de résolution du problème du voyageur de commerce
# Entrée : une matrice de distances
# Sortie : un couple composé d'une séquence de visites (ne pas oublier le retour au "premier" lieu (dépôt)) et de sa longueur
function solveTSPExact(d::Matrix{Int64})
    # Déclaration de variables (Julia oblige de mettre une valeur initiale même si elle n'a aucun sens)
    cycle::Vector{Int64} = []
    l::Int64 = 0
    taille::Int64 = size(d,1)
    # Cas triviaux
    if taille <= 3
        cycle = [i for i in 1:taille]
        l = d[taille,1]
        for i in 1:taille-1
            l += d[i,i+1]
        end
    else
        cycle, l = TravelingSalesmanExact._get_optimal_tour(d, GLPK.Optimizer, true, false, true)
    end

    # Retour du résultat
    return cycle, round(Int64,l)
end

# fonction de prise des données et de résolution
function data_then_solve_exact(filename::String)
    # Conversion fichier -> structures de données
    data::donnees = lecture_donnees(filename)
    nbClients::Int64 = data.nbClients
    distancier::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande

    # déterminer l'ensemble des regroupements possibles
    println("Temps CPU pour la construction de l'instance de partitionnement d'ensemble :")
    S::Vector{Vector{Int64}} = @time getSubsets(capa,dmd,distancier)
    println()
    println("Nombre de regroupements : ", length(S))
    println()

    # vecteur des couples ordres de tournées / longueurs
    println("Temps CPU pour fournir le vecteur des couples (tournee,dist_min) pour chaque regroupement via TSP :")
    l::Vector{Tuple{Vector{Int64},Int64}} = @time getAllShortestCycles(S,distancier)
    println()

    # vecteur des vecteurs des indices de tournées contenant le client i
    AllS_i::Vector{Vector{Int64}} = []
    for i in 2:nbClients
        AllS_i = push!(AllS_i,getSetofCyclesClient(S,i))
    end

    # nb de regroupements
    nbRegroup::Int64 = length(S)

    # création du modèle à partir des données
    m::Model = model_exact(GLPK.Optimizer,AllS_i,nbClients,l,nbRegroup)

    # résolution
    println("Temps CPU pour la résolution du programme linéaire :")
    @time optimize!(m)
    println()

    # Affichage des résultats (ici assez complet pour gérer certains cas d'"erreur")
    status = termination_status(m)

    if status == MOI.OPTIMAL
        println("Problème résolu à l'optimalité")

        # Affichage de la tournée retenue
        println("Tournées retenues : ")
        for j in 1:nbRegroup
            if (value(m[:x][j]) == 1.0)
                println("n°",j," : ", l[j][1], " de distance ", l[j][2])
            end
        end

        # affichage de la valeur optimale
        println("Distance totale minimale : z = ", round(Int64,objective_value(m)))

        println()

    elseif status == MOI.INFEASIBLE
        println("Problème non-borné")

    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Problème impossible")
    end

    # dernier affichage avant celui de la macro @time dans timer_res_exact
    println("Temps CPU total :")

end

# Fonction à exécuter dans le REPL : exécution de data_then_solve, avec le timer
function timer_res_exact(filename::String)
    @time data_then_solve_exact(filename)
end




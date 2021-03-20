#= méthode exacte =#

using JuMP, GLPK
using TravelingSalesmanExact # pour des fonctions du problème du voyageur de commerce
using Test #Pour faire les tests

# Structure contenant les données du problème
mutable struct donnees
    nbClients::Int64 # Nombre de clients (y compris le dépôt)
    capacite::Int64 # Capacité du véhicule de livraison
    demande::Vector{Int64} # Demande de chaque client
    distance::Matrix{Int64} # Distancier (Matrix{Int64} est équivalent à Array{Int64,2})
end

# Fonction de lecture des données du problème
# TODO : modifier les commentaires et les variables
function lecture_donnees(nom_fichier::String)
    # Ouverture d'un fichier en lecture
    f::IOStream = open(nom_fichier,"r")

    # Lecture de la première ligne pour connaître le nombre de villes
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


function model_exact(solverSelected::DataType, AllS_i::Vector{Vector{Int64}}, nbClient::Int64, l::Vector{Int64}, nbRegroup::Int64)

    # Déclaration d'un modèle (initialement vide)
    m::Model = Model(solverSelected)

    n::Int64 = nbClient;

    # Déclaration des variables de décision
    # si on choisit la tournée indicée j dans l'ensemble S = \mathcal(S)
    @variable(m, x[1:nbRegroup]>=0, binary = true)

    # Déclaration de la fonction objectif (avec le sens d'optimisation)
    @objective(m, Min, sum(l[j]x[j] for j in 1:nbRegroup))

    # Déclaration des contraintes : A REPRENDRE à cause d'un pb d'index
    @constraint(m, VisitOnlyOnceClient[i=2:n], sum(x[j] for j in AllS_i[i-1]) == 1)

    return m
end

#Fonction qui prend deux ensemble  initialement vide, la capacite de drone, 
#une vecteur avec les demandes des endroitrs et deux integer,
#qui doit être 0 à la debut et retourne l'ensemble S qui contient
#tous les combinaison possible une fois 
function getSubsets_recursive(P::Vector{Int64}, S::Vector{Vector{Int64}},capacite::Int64, demande::Vector{Int64}, index::Int64, d::Int64, distances::Matrix{Int64})
    #Arreter la recursion si on est à la fin d'index et retour S
    if index > length(demande)
        return S                        
    end
    #buckle sur tout les demandes de clients
    while index <= length(demande) 
        if d+demande[index]<=capacite                                                                                       
            toadd::Vector{Int64}= copy(P)                   #Copy P pour avoir un nouveau tableau et les changements dans toadd ne soit pas transfere à Par 
            toadd=append!(toadd,index+1)                     #Comme demande des elements dans P et demande[index] sont ensemble plus petite que la capacité on construit l'union de les deux
            S=vcat(S,[toadd])                                 #On ajoute cette union dans S
            Snew::Vector{Vector{Int64}}=getSubsets_recursive(toadd,S,capacite,demande,index+1,demande[index]+d,distances)  #appelle recursive avec cette union le nouveau demande est la demande de tous elements dans toadd et on augmenter l'index avec 1 parce que on veut pas considere le même element
            if length(S)<length(Snew)                         #On veut seulement ajoute le fin de resultat de l'apelle recursice si les deux arrays sont pas egale
                S=vcat(S,Snew[length(S)+1:length(Snew)])        #Si ils sont pas égale on append le part d'array de l'appelle récursive qui est nouveau 
            end
        end
        index+=1
    end
    return S
end

#wrapper pour la méthod getSubsets_recursive
function getSubsets(capacite::Int64,demande::Vector{Int64}, distances::Matrix{Int64})
    P::Vector{Int64}=[]
    S::Vector{Vector{Int64}}=[]
    return collect(getSubsets_recursive(P,S,capacite,demande,1,0,distances))
end

# déterminer l'ensemble de numéros de regroupements dans lesquels le client "cli" est livré
function getSetofCyclesClient(S::Vector{Vector{Int64}}, cli::Int64)
    res::Vector{Int64} = []
    for i in 1:length(S) # pour chaque regroupement...
        if cli in S[i]  # ...si notre client "cli" y appartient...
            res = push!(res,i)    # ...alors on ajoute l'indice de ce regroupement dans notre vecteur final
        end
    end
    return res
end


# déterminer la distance la plus courte dans un regroupement Si, étant donné un distancier d
function determineShortestCycle(Si::Vector{Int64}, d::Matrix{Int64}) # S est l'ensemble d'indices et d est le distancier
    Si = append!([1],Si)
    newd::Matrix{Int64} = d[Si,Si]
    return solveTSPExact(newd)[2]
end

# fournir le vecteur des longueurs de chaque regroupement dans S
function getAllShortestCycles(S::Vector{Vector{Int64}}, d::Matrix{Int64})
    res::Vector{Int64} = [] # initialisation
    for i in 1:length(S)    # pour chaque regroupement
        res = push!(res,determineShortestCycle(S[i],d)) # ajouter sa distance minimale au vecteur à retourner
    end
    return res
end


function test()
    #First test all subsets and respective shortest distances
    data::donnees = lecture_donnees("exemple.dat") # fichier dans le même dossier (cf ex. du sujet)
    d::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande
    nbClients::Int64 = data.nbClients
    S::Vector{Vector{Int64}} = getSubsets(capa,dmd,d)
    l::Vector{Int64} = getAllShortestCycles(S,d)
    AllS_i::Vector{Vector{Int64}} = []
    for i in 2:nbClients
        AllS_i = push!(AllS_i,getSetofCyclesClient(S,i))
    end
    #=
    println("n=",nbClients)
    println(S)
    println(l)
    println(AllS_i)
    =#
    for i in 2:nbClients
        for j in AllS_i[i-1]
            println(j)
        end
    end

#=
    #seconde test: get getSubsets
    de::Vector{Int64}=[2,4,2,4,2]
    cap::Int64=10
    testset::Set=Set([2])
    Es::Array{Array{Int64}}=getSubsets(cap,de,d)  
    @testset "method tests" begin
       @test determineShortestCycle(Set)
    end;
    print(Es)
    =#
end


function data_then_solve(filename::String)
    # Conversion fichier -> structures de données

    data::donnees = lecture_donnees(filename)
    nbClients::Int64 = data.nbClients
    distancier::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande

    # déterminer l'ensemble des regroupements possibles
    S::Vector{Vector{Int64}} = getSubsets(capa,dmd,distancier)
    
    # vecteur des longueurs
    l::Vector{Int64} = getAllShortestCycles(S,distancier)

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
    optimize!(m)

    # Affichage des résultats (ici assez complet pour gérer certains cas d'"erreur")
    status = termination_status(m)

    if status == MOI.OPTIMAL
        println("Problème résolu à l'optimalité")

        # Affichage de la tournée retenue
        println("Tournées retenues : ")
        for j in 1:nbRegroup
            if (value(m[:x][j])==1.0)
                println("n°",j," : ", S[j], " de distance ", l[j])
            end
        end

        # affichage de la valeur optimale
        println("Distance totale minimale : z = ", objective_value(m))

        println()

    elseif status == MOI.INFEASIBLE
        println("Problème non-borné")

    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Problème impossible")
    end

end



# Fonction de résolution du problème du voyageur de commerce
# Entrée : une matrice de distances
# Sortie : un couple composé d'une séquence de visites (ne pas oublier le retour à la "première" ville) et de sa longueur
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

#=
# Exemple d'application de solveTSPExact() sur le problème de l'exercice 2.5
# Valeur retournée : ([7, 1, 5, 6, 2, 3, 4], 2575.0)
# Cela signifie que la problème est résolu à l'optimalité avec le cycle 7 -> 1 -> 5 -> 6 -> 2 -> 3 -> 4 -> 7
# La longueur de ce cycle est 2575
function solveEx25()
    d::Matrix{Int64} = [
    0 786 549 657 331 559 250;
    786 0 668 979 593 224 905;
    549 668 0 316 607 472 467;
    657 979 316 0 890 769 400;
    331 593 607 890 0 386 559;
    559 224 472 769 386 0 681;
    250 905 467 400 559 681 0
    ]
    return solveTSPExact(d)
end
=#


#= Il existe plusieurs façons (plus ou moins efficaces) de réaliser les implémentations demandées.
Des possibilités offertes par Julia sont présentées dans les quelques lignes qui suivent.
Suivant la façon dont vous voulez aborder les implémentations, vous n'êtes pas tenu de les exploiter.
Tout d'abord, on peut créer un tableau à partir d'un autre tableau.
Par exemple :
    T = [1,2,4,8,16]
    T2 = T[2:4] # On a alors T2 = [2,4,8]
    T3 = T[[1,2,4]] # On a alors T3 = [1,2,8]
On peut de la même façon créer une matrice à partir d'une sous matrice
Par exemple :
    d = [
    0 786 549 657 331 559 250;
    786 0 668 979 593 224 905;
    549 668 0 316 607 472 467;
    657 979 316 0 890 769 400;
    331 593 607 890 0 386 559;
    559 224 472 769 386 0 681;
    250 905 467 400 559 681 0
    ]
    d2 = d[1:4,1:4] # On a alors d2 = [   0  786  549  657;
                                        786    0  668  979;
                                        549  668    0  316;
                                        657  979  316    0]
    d3 = d[[1,4,6],[1,4,6]] # On a alors d2 = [   0  657  559
                                                657    0  769;
                                                559  769    0]
En ce qui concerne la manipulation des tableaux, la fonction vcat pourra être utile pour créer un tableau créé en concaténant deux ou plusieurs tableaux.
La fonction append! qui modifie un tableau en ajoutant à la fin les valeurs se trouvant dans un autre tableau pourra également être utile.
Le fonction reverse! dont le principe est d'inverser l'ordre des éléments d'un tableau pourra aussi être intéressante.
Enfin, la fonction sort! qui sert à faire des tris sera sans doute d'un intérêt évident.
Pour rappel, on peut trouver de l'aide sur ces fonctions en tapant '?', puis le nom de la fonction dans le REPL.

Une dernière chose très importante : la mesure du temps !
En Julia, on peut mesurer le temps pris par une fonction très simplement à l'aide de la macro @time.
Par exemple,
    @time f(...)
permet de mesurer le temps d'exécution d'une fonction f (les points de suspension correspondent aux paramètres de la fonction).
Il faudra bien entendu que cette fonction ait été exécutée au moins une fois au préalable, la compilation s'effectuant lors de la première exécution.
=#

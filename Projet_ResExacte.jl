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
#Structure pour le chemins
struct chemin
    clients::Set
    longeur::Int64
end



# Fonction de lecture des données du problème
# TODO : modifier les commentaires et 3es variables
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

# TODO : remplacer Set par Array pour l'"ensemble" des regroupements

function model_exact(solverSelected::DataType, S::Set, nbClient::Int64)

    # Déclaration d'un modèle (initialement vide)
    m::Model = Model(solverSelected)

    cardS::Int64 = length(S);

    n::Int64 = nbClient;

    l::Vector{Int64} = getAllShortestCycles(S,d)

    # Déclaration des variables de décision
    # si on choisit la tournée indicée j dans l'ensemble S = \mathcal(S)
    @variable(m, x[1:cardS]>=0, binary = true)

    # TODO : structure de données associant S et les longueurs respectives (cf Yannick)

    # Déclaration de la fonction objectif (avec le sens d'optimisation)
    @objective(m, Min, sum(l[j]x[j] for j in 1:cardS))

    # Déclaration des contraintes
    # TODO (deprecated) : déclarer SetIndicesTournees[i] = S_i \subset [cardS] d'indices de tournées contenant le client i
    @constraint(m, VisitOnlyOnceClient[2:n], sum(x[j] for j in getSetofCyclesClient(S,i)) == 1)

end
#Fonction qui prend deux ensemble  initialement vide, la capacite de drone, 
#une vecteur avec les demandes des endroitrs et deux integer,
#qui doit être 0 à la debut et retourne l'ensemble S qui contient
#tous les combinaison possible une fois 
function getSubsets_recursive(P::chemin, S::Set,capacite::Int64, demande::Vector{Int64}, index::Int64, d::Int64, distances::Matrix{Int64})
    if index > length(demande)
        return S
    end
    while index <= length(demande) 
        if d+demande[index]<=capacite
            toadd::Set= union(P.clients,Set([index+1]))
            print(toadd)
            chemin_to_add=chemin(toadd,determineShortestCycle(toadd,distances))
            S=push!(S,chemin_toadd)
            S=union(S,getSubsets_recursive(chemin_toadd,S,capacite,demande,index+1,demande[index]+d))
        end
        index+=1
    end
    return S
end
#wrapper pour la méthod getSubsets_recursive
function getSubsets(capacite::Int64,demande::Vector{Int64}, distances::Matrix{Int64})
    ch::chemin=chemin(Set([]),0)
    return collect(getSubsets_recursive(ch,Set([]),capacite,demande,1,0,distances))
end

# déterminer l'ensemble de numéros de regroupements dans lesquels le client "cli" est livré
function getSetofCyclesClient(S::Set, cli::Int64)
    res::Set = Set([])
    for s in S
        if cli in s
            res = union(res,Set([cli]))
        end
    end
end

# déterminer le minimum (sauf zéro) et son indice dans un array (ligne d'une matrice ou d'un array en dim 2)
function minNonZero(a::Array{Int64,2}, S::Set)
    tmp::Int64 = typemax(Int64) # valeur arbitrairement grande
    ind::Int64 = 1
    for i in S
        if (a[i]<tmp && a[i] !== 0) # a[i] !== 0 est redondante car : a[i] == 0 <=> boucle en i. Or un tel i n'est pas dans S car il est exclu en pré-traitement (cf fonction dSC)
            tmp = a[i]
            ind = i
        end
    end
    return tmp, ind
end

# déterminer la distance la plus courante dans un regroupement S, étant donné un distancier d
function determineShortestCycle(S::Set, d::Matrix{Int64}) # S est l'ensemble d'indices et d est le distancier

    m::Tuple{Int64,Int64} = minNonZero(d[1:1,:],S)
    println(m)
    dtot::Int64 = m[1]      # min de la 1ere ligne pour avoir le trajet minimal de 1 vers un certain lieu i
    k::Int64 = m[2]         # récupération de l'indice concerné (lieu à visiter)
    S = setdiff(S,Set([k])) # enlever le lieu visité
    while (length(S) > 1)   # tant qu'il ne reste pas qu'un seul lieu dans S
        m = minNonZero(d[k:k,:],S)
        println(m)
        dtot = dtot + m[1]
        k = m[2]
        S = setdiff(S,Set([k]))
    end

    # récupérer l'élément restant dans S (pas d'accès direct dans un Set, d'où l'utilisation d'une boucle)
    l::Int64 = 0 # initialisation de l à l'extérieur de la boucle
    for z in S
        l = z
    end

    dtot = dtot + d[k,l] + d[l,1] # sommer les dernières distances

    return dtot
end

# fournir le vecteur des longueurs de chaque regroupement
# TODO : remplacer le type de l'ensemble des regroupements par Array (car Set n'est pas ordonné)
function getAllShortestCycles(S::Set, d::Matrix{Int64})
    res::Vector{Int64} = []
    for s in S
        res = push(res,determineShortestCycle(s,d))
    end
end


function test()
    #First test determineShortestPast
    data::donnees = lecture_donnees("exemple.dat") # fichier dans le même dossier (cf ex. du sujet)
    d::Matrix{Int64} = data.distance
    S::Set = Set([2,3,6])
    dtot::Int64 = determineShortestCycle(S,d)
    println(dtot)
    



    #seconde test: get getSubsets
    de::Vector{Int64}=[2,4,2,4,2]
    cap::Int64=10
    Es::Array{chemin}=getSubsets(cap,de,d)  
    @testset "method tests" begin
        @test typeof(Es)==Array{Set,1}
        @test dtot == 787
        @test 
    end;
    typeof(Es)
    for ch in Es
        typeof(ch)
        #print("Clients:")
        #print(ch.clients)
        #print("longeur:")
        #print(ch.longeur)
    end
end


#=
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

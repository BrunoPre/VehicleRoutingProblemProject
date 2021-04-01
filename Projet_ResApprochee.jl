#= PROJET RO 2021
    Résolution approchée
    PREKA Bruno
    ZELLE Yannick
    681C
=#

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

# A MODIFIER Fonction de modèle
function model_appr(solverSelected::DataType, AllS_i::Vector{Vector{Int64}}, nbClient::Int64, l::Vector{Tuple{Vector{Int64},Int64}}, nbRegroup::Int64)
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

# Fonction calculant le vecteur des gains, étant donné un distancier
function calcGainVector(d::Matrix{Int64})

    # parcours dans la partie triangulaire supérieure (par symétrie du distancier)
    n,m::Int64 = size(d)
    G::Vector{Tuple{Tuple{Int64,Int64},Int64}} = []
    k::Int64 = 3    # indice de décalage du début de la lecture de la ligne, suivant un zéro (parcours triangulaire)
    gij::Int64 = 0 # déclaration de la variable
    
    for i in 2:n # parcours selon les lignes de d (on commence à 2, d'où, par ailleurs, k=3 au début)
        # parcours selon les colonnes à partir de k
        for j in k:m
            gij = d[i,1] + d[1,j] - d[i,j] # gain
            G = push!(G, ((i,j),gij) ) # ajout de la paire (indices,gain) au vecteur de gains
        end
        k += 1 # décaler le début de la lecture du triangle supérieur de d
    end
    return G
end

#=
#fonction prend la solution et réalise la fusion de chemins
function fusioner(chemins::Vector{Vector{Int64}}, sortedGains::Vector{Tuple{Tuple{Int64,Int64},Int64}}, demande::Vector{Int64}, capacite::Int64)
    saveFusions::Vector{Int64} =[]   #Vecteur pour garder l'information si deux chemins sont fusionnés
    index::Int64=1
    for ch in chemins           #Initialisation de ce vecteur. Tout d'abord chaque indice correspond au chemin qui lui est assossié
        append!(saveFusions,index)
        index+=1
    end
    for gain in sortedGains                 # on itère sur tous les gains (toujours dans le sens décroissant)
        fusionA = saveFusions[gain[1][1]-1]       # pour chaque gain, on prend les indices (i,j) associés à ce gain
        fusionB = saveFusions[gain[1][2]-1]
        println(fusionA != fusionB && demande[fusionA]+demande[fusionB]<=capacite)
        if (fusionA != fusionB && demande[fusionA]+demande[fusionB]<=capacite)       #si fusionA=fusionB, alors la fusion a déjà eu lieu. On ne fait rien, même dans le cas où la demande de la fusion dépasse la capacité
            saveFusions[gain[1][2]-1] =fusionA                                        #Autrement, on fait la fusion et on remplace la position "fusionB" de saveFusions par "fusionA"
            demande[fusionA]= demande[fusionA]+demande[fusionB]                     #On actualise aussi la demande,
            append!(chemins[fusionA],chemins[fusionB])                    #puis le chemin qui correspond à cette fonction
            chemins[fusionB]=[-99]
            i::Int64=1
            while i <= size(saveFusions)[1]
                if (saveFusions[i]==fusionB)
                    saveFusions[i]=fusionA
                end
                i+=1
            end
        end    
    end                                            #Le chemin à la position B n'existe plus. Il est maintenant inclus dans le chemin à la position fusionA. Alors on le supprime
    return chemins
end




function formatSolution(solution::Vector{Vector{Int64}})
    formatedSolution::Vector{Vector{Int64}}=[]
    for sol in solution 
        if (sol!= [-99])
            sol=append!([1],sol)
            append!(sol,[1])
            push!(formatedSolution,sol)
        end
    end
    return formatedSolution

end
=#




# fonction de prise des données et de résolution
function data_then_solve_exact(filename::String)
    # Conversion fichier -> structures de données
    data::donnees = lecture_donnees(filename)
    nbClients::Int64 = data.nbClients
    distancier::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande

    # vecteur des gains ((i,j),gij)
    G::Vector{Tuple{Tuple{Int64,Int64},Int64}} = calcGainVector(d)

    # tri décroissant par gij
    G = sortVector(G)
    
    # Initialisation des tournées élémentaires pour la fusion
    elementaryCycles::Vector{Tuple{Vector{Int64},Int64}}= initialiserChemins(dmd)

    # tournées fusionnées
    fusionnedCycles::Vector{Tuple{Vector{Int64},Int64}}= getAllFusionnedCycles(G,elementaryCycles, dmd,capa)

    # formattage des tournées
    fusionnedCycles = formatCycles(fusionnedCycles)


    #=
    # A MODIFIER création du modèle à partir des données
    m::Model = model_appr(GLPK.Optimizer,AllS_i,nbClients,l,nbRegroup)

    # résolution
    optimize!(m)

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
=#
end

# exécution de data_then_solve, avec le timer
function timer_res_exact(filename::String)
    @time data_then_solve_exact(filename)
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

#= deuxième version comme dans le sujet
Les structures de données employéees pour la résolution approchée sont les suivantes : 
* Retournée par la fonction calcGainVector, le vecteur "ordered" contient toutes les fusions possibles ((i,j),gij), 
triées dans l'ordre décroissant des gains ;
* Le vecteur des tuples "sol". Un tuple contient une solution de tournée et un entier, et il décrit si la solution a déjà subi une fusion.
    Le cas échéant, le tuple renseigne avec quelle autre tournée la fusion a eu lieue :
    ** La première composante du tuple est une solution de tournée (vecteur) de la forme [1,...,1]. 
    ** La deuxième composante est un entier qui indique avec quelle autre tournée cette solution a été fusionée. 
    Dans ce cas, la première composante est un vecteur vide. Initialement les entiers qui portent l'information de la fusion sont mis à 0.
=#
function getAllFusionnedCycles(ordered::Vector{Tuple{Tuple{Int64,Int64},Int64}}, sol::Vector{Tuple{Vector{Int64},Int64}}, demande::Vector{Int64}, capacité::Int64)

    for el in ordered                   #itération sur toutes les fusions possibles
        firstEl::Int64=el[1][1]-1       #premier client
        secondEl::Int64=el[1][2]-1      #deuxième client
        indexFirst::Int64 = getIndex(firstEl,sol)   #trouver l'index de la fusion dans laquelle se trouve le premier client (il peut différer à cause de précédentes fusions)
        indexSecond::Int64=getIndex(secondEl,sol)   #Même action pour le deuxième client
        lengthOld::Int64=size(sol[indexFirst][1])[1]   #sauvegarder la taille de la fusion pour savoir après s'il faut enlever des fusions

        if(demande[indexFirst]+demande[indexSecond]<=capacité) #si la demande est plus petite que la capacité, alors on fusionnne
            fusioned::Vector{Int64} = fusion(sol[indexFirst][1], sol[indexSecond][1], firstEl+1, secondEl+1)    # réaliser cette fusion (le résultat est un array qui contient la nouvelle solution) 
            if(fusioned[2]==firstEl+1)      #ceci signifie que fusionned[2] contient le premier client, donc on veut pas enlever cette solution
                firstEl=0
            end
            if(fusioned[lengthOld+1]==1)    #idem pour le deuxième élément
                secondEl=0
            end
            if (size(fusioned)[1]>4)        #Si la taille dépasse 4,alors on enlève les clients (a,b) (a < b) dans la fusion = [1,a,...,b,1]
                ordered=removeFusions(firstEl,secondEl,(fusioned[2],fusioned[size(fusioned)[1]-1]),ordered)
            end
            overwriteFirst=(fusioned,0)                             #la nouvelle solution 
            overwriteSecond=([],indexFirst)                         #on actualise aussi l'élément qui a été fusionné avec l'index indiquant où se trouve l'élement suite à la fusion
            sol[indexFirst]=overwriteFirst                          #on garde ces informations dans l'array des solutions
            sol[indexSecond]=overwriteSecond
            demande[indexFirst]=demande[indexFirst]+demande[indexSecond]        #on actualise la demande
        end
        
    end
    return sol
end

#fonction récursive donnant l'index de la solution dans laquelle un certain élement "el" a été fusionné
function getIndex(el::Int64, sol::Vector{Tuple{Vector{Int64},Int64}})           
    fusion=sol[el][2]
    if (fusion==0)
        return el
    end
    return getIndex(fusion,sol)
end

#fonction réalisant une fusion
function fusion(ai::Vector{Int64},bj::Vector{Int64},i::Int64,j::Int64)   

    if(ai[2]==i)                #alors ai est de la forme [1,i,...,1]
        if(bj[size(bj)[1]-1]==j && size(ai)[1]!=3)    #si bj est de la forme [1,...,j,1] on retourne b+a, sauf dans le cas où ai = [1,i,1] (tournée élémentaire)
                return append!(bj[1:size(bj)[1]-1],ai[2:size(ai)[1]]) 
        end                                     #sinon, bj est forcément de la forme [1,j,...,1]      
        ai=reverse(ai)                           #ainsi on inverse ai et on retourne la concaténation de ai et bi
        return append!(ai[1:size(ai)[1]-1],bj[2:size(bj)[1]])
    end

    if(bj[size(bj)[1]-1]==j)                    #ici on sait que ai est de la forme [1,...,i,1]. Le seul cas qui peut être problématique est si bj est de la forme [1,...,j,1]. Dans ce cas, on inverse b
        bj=reverse(bj)
    end
    return append!(ai[1:size(ai)[1]-1],bj[2:size(bj)[1]])   #après on retourne la concaténation de ai et bi
end
            
#= fonction prenant les fusions de "ordered" et trois arguments et
elle itère une fois sur tous les fusions et elle enlève les fusions données par les arguments "removeA", "removeB", "removeC". 
Pour chque fusion où il y a plus de 5 éléments,
si le résultat est de la forme [1,a,...,b,1], alors on enlève au moins les clients (a,b) donnés par "removeC". 
Si on a fait une fusion de la forme [1,...,i,j,...,1], alors on peut enlever aussi i et j, respectivement donnés par
"removeA" et "removeB". 
Si on veut pas enlever ces élements, "removeA" et "removeB" sont mis à 0 dans la fonction "getAllFusionnedCycles"
=#

function removeFusions(removeA::Int64, removeB::Int64, removeC::Tuple{Int64,Int64}, ordered::Vector{Tuple{Tuple{Int64,Int64},Int64}})

    if (removeC[1]>removeC[2])      #il faut que le tuple soit de la forme (a,b), avec a<b
        removeC=(removeC[2],removeC[1])            # si ce n'est pas le cas, on renverse le couple             
    end

    i::Int64=1                                                    #index pour supprimer la fusion qui correspond

    for el in ordered
        #Est-ce que removeA ou removeB sont dans le couple de la fusion regardée ? ou est-ce que le couple est égale à removeC ?
        if(el[1][1]==removeA || el[1][1]==removeA || el[1][2]==removeA|| el[1][2]==removeB || el[1]==removeC)           
            ordered= deleteat!(ordered, i)                                                                              #si oui, on supprime ce couple
        else
            i+=1
        end
    end
    return ordered
end

# fonction qui initialise les chemins du type 1->X->1 (X \in [|2,n|] ) et donne 0 comme deuxième argument, parce qu'il n'y a aucune fusion au départ
function initialiserChemins(demande::Vector{Int64})
    index::Int64=2
    chemins::Vector{Tuple{Vector{Int64},Int64}}=[]
    for d in demande
        chemin::Vector{Int64}=[1,index,1]
        toPush::Tuple{Vector{Int64},Int64}=(chemin,0)
        push!(chemins,toPush)
        index+=1
    end
    return chemins
end

#fonction de mise en forme des tournées 
function formatCycles(allCycles::Vector{Tuple{Vector{Int64},Int64}})
    newallCycles::Vector{Vector{Int64}} = []
    for cycle in allCycles
        if (cycle[2]==0)
            push!(newallCycles,sol[1])
        end
    end
    return newallCycles
end

# fonction triant le vecteur ((i,j),gij) dans l'ordre décroissant des gij
function sortVector(toSort::Vector{Tuple{Tuple{Int64,Int64},Int64}})
    return sort!(toSort,by = x -> x[2],rev=true)
end

# fonction de test (hors-sujet)
function test()
    #First test all subsets and respective shortest distances
    data::donnees = lecture_donnees("A/VRPA15.dat") # fichier dans le même dossier (cf ex. du sujet)
    d::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande
    nbClients::Int64 = data.nbClients

    G::Vector{Tuple{Tuple{Int64,Int64},Int64}} = calcGainVector(d)
    println(G)
    print(sortVector(G))

    G = sortVector(G)
    println(G)
    println()
    
    sol::Vector{Tuple{Vector{Int64},Int64}}= initialiserChemins(dmd)
    print(sol)
    println()


    solution::Vector{Tuple{Vector{Int64},Int64}}= getAllFusionnedCycles(G,sol, dmd,capa)
    println(solution)
    println("Les resultats sont")
    for sol in solution
        if (sol[2]==0)
            println()
            println(sol[1])
        end
    end
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

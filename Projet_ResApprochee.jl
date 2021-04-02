#= PROJET RO 2021
    Résolution approchée
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


#= fonction déterminant toutes les tournées fusionnées possibles
Les structures de données employées pour la résolution approchée sont les suivantes : 
* Retournée par la fonction calcGainVector, le vecteur "ordered" contient toutes les fusions possibles ((i,j),gij), 
triées dans l'ordre décroissant des gains ;
* Le vecteur des tuples "sol". Un tuple contient une solution de tournée et un entier, et il décrit si la solution a déjà subi une fusion.
    Le cas échéant, le tuple renseigne avec quelle autre tournée la fusion a eu lieue :
    ** La première composante du tuple est une solution de tournée (vecteur) de la forme [1,...,1]. 
    ** La deuxième composante est un entier qui indique avec quelle autre tournée cette solution a été fusionée. 
    Dans ce cas, la première composante est un vecteur vide. Initialement les entiers qui portent l'information de la fusion sont mis à 0.
=#
function getAllFusionnedCycles(ordered::Vector{Tuple{Tuple{Int64,Int64},Int64}}, sol::Vector{Tuple{Vector{Int64},Int64}}, demande::Vector{Int64}, capacité::Int64)
    i::Int64=1                          #index pour donner la position de la fusion actuelle à la fonction "fusion"

    for el in ordered                   #itération sur toutes les fusions possibles
        firstEl::Int64=el[1][1]-1       #premier client
        secondEl::Int64=el[1][2]-1      #deuxième client
        indexFirst::Int64 = getIndex(firstEl,sol)   #trouver l'index de la fusion dans laquelle se trouve le premier client (il peut différer à cause de précédentes fusions)
        indexSecond::Int64=getIndex(secondEl,sol)   #Même action pour le deuxième client
        lengthFirst::Int64=size(sol[indexFirst][1])[1]   #sauvegarder la taille de la fusion pour savoir après s'il faut enlever des fusions
        lengthSecond::Int64=size(sol[indexSecond][1])[1] 

        if(demande[indexFirst]+demande[indexSecond]<=capacité) #si la demande est plus petite que la capacité, alors on fusionnne
            fusioned::Vector{Int64} = fusion(sol[indexFirst][1], sol[indexSecond][1], firstEl+1, secondEl+1)    # réaliser cette fusion (le résultat est un array qui contient la nouvelle solution) 
            if(lengthFirst==3)      #ceci signifie que fusioned[2] contient le premier client, donc on veut pas enlever cette solution
                firstEl=-1
            end
            if(lengthSecond==3)    #idem pour le deuxième élément
                secondEl=-1
            end
            if (size(fusioned)[1]>4 && i+1!= size(ordered)[1])        #Si la taille dépasse 4,alors on enlève les clients (a,b) (avec a < b) dans la fusion = [1,a,...,b,1]
                ordered=removeFusions(firstEl+1,secondEl+1,(fusioned[2],fusioned[size(fusioned)[1]-1]),ordered, i+1)
            end
            overwriteFirst=(fusioned,0)                             #la nouvelle solution 
            overwriteSecond=([],indexFirst)                         #on actualise aussi l'élément qui a été fusionné avec l'index indiquant où se trouve l'élement suite à la fusion
            sol[indexFirst]=overwriteFirst                          #on garde ces informations dans l'array des solutions
            sol[indexSecond]=overwriteSecond
            demande[indexFirst]=demande[indexFirst]+demande[indexSecond]        #on actualise la demande
        end

        i+=1
        
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
    sia= size(ai)[1]
    sib= size(bj)[1]

    # dans les commentaires suivants, on note la concaténation de ai et bj par "ai+bj"    

    # cas trivial où ai = [1,i,1]
    if (sia==3)
        if( bj[2]==j)                                                   #si bj = [1,j,...,1]
            return append!(ai[1:sia-1],bj[2:sib])                       #Alors on retoure ai+bj=[1,i,j,...,1]
        end
        return append!(bj[1:size(bj)[1]-1],ai[2:size(ai)[1]])           # sinon, i.e. bj = [1,...,j,1]. Alors on retourne bj+ai=[1,...,j,i,1]
    end

    # si ai = [1,...,i,1]
    if (ai[2]!=i)                                                       # où "..." n'est pas vide, i.e. ai n'est pas une tournée élémentaire 
        if (bj[2]==j)                                                     
            return append!(ai[1:sia-1],bj[2:sib])                       #Alors si bj= [1,j,...,1] on retourne ai+bj=[1,...,i,j,...,1]
        end
        bj=reverse(bj)                                                  
        return append!(ai[1:sia-1],bj[2:sib])                           #Ici  ai = [1,...,i,1] et bj = [1,...,j,1] alors on retourne ai+inverse(bj)= [1,...,i,j,...1]
    end    

    # sinon ai = [1,i,...,1] n'est pas élémentaire
    # et si bj = [1,...,j,1]
    if (bj[sib-1]==j)
        return append!(bj[1:size(bj)[1]-1],ai[2:size(ai)[1]])           # on retourne bj+ai= [1,...,j,i,...,1]
    end

    # toujours si ai = [1,i,...,1] n'est pas élémentaire, et si bj = [1,j,...,1]
    ai=reverse(ai)                              
    return append!(ai[1:sia-1],bj[2:sib])                     #on retourne inverse(ai)+bj= [1,...,i,j,...,1]                            

end
            
#= fonction prenant les fusions de "ordered" et trois arguments et
elle itère une fois sur tous les fusions et elle enlève les fusions données par les arguments "removeA", "removeB", "removeC". 
Pour chque fusion où il y a plus de 5 éléments,
si le résultat est de la forme [1,a,...,b,1], alors on enlève au moins les clients (a,b) donnés par "removeC". 
Si on a fait une fusion de la forme [1,...,i,j,...,1], alors on peut enlever aussi i et j, respectivement donnés par
"removeA" et "removeB". 
Si on veut pas enlever ces élements, "removeA" et "removeB" sont mis à 0 dans la fonction "getAllFusionnedCycles"
=#

function removeFusions(removeA::Int64, removeB::Int64, removeC::Tuple{Int64,Int64}, ordered::Vector{Tuple{Tuple{Int64,Int64},Int64}}, index::Int64)

    if (removeC[1]>removeC[2])      #il faut que le tuple soit de la forme (a,b), avec a<b
        removeC=(removeC[2],removeC[1])            # si ce n'est pas le cas, on renverse le couple             
    end
                                                
    while index < size(ordered)[1]                          #l'index est donné par la fonction "getAllFusionnedCycles", ainsi on ne regarde pas les solutions qui étaient déjà fusionnées
        #Est-ce que removeA ou removeB sont dans le couple de la fusion regardée ? ou est-ce que le couple est égale à removeC ?
        if(ordered[index][1][1]==removeA || ordered[index][1][1]==removeA || ordered[index][1][2]==removeB|| ordered[index][1][2]==removeB || ordered[index][1]==removeC)           
            ordered= deleteat!(ordered, index)                                                                              #si oui, on supprime ce couple
        else
            index+=1
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
            push!(newallCycles,cycle[1])
        end
    end
    return newallCycles
end

# fonction triant le vecteur ((i,j),gij) dans l'ordre décroissant des gij
function sortVector(toSort::Vector{Tuple{Tuple{Int64,Int64},Int64}})
    return sort!(toSort,by = x -> x[2],rev=true)
end

#calcul de la longueur totale d'une tournée "cycle", étant donné le distancier "d"
function calcLengthOfCycle(cycle::Vector{Int64}, d::Matrix{Int64})
    sum::Int64 = 0
    i::Int64 = 1
    while i+1 <= length(cycle)
        sum = sum + d[cycle[i],cycle[i+1]]
        i = i+1
    end
    return sum
end

# Fonction à exécuter dans le REPL : fonction de prise des données et de résolution
function data_then_solve_appr(filename::String)
    # Conversion fichier -> structures de données
    data::donnees = lecture_donnees(filename)
    nbClients::Int64 = data.nbClients
    distancier::Matrix{Int64} = data.distance
    capa::Int64 = data.capacite
    dmd::Vector{Int64} = data.demande

    # vecteur des gains ((i,j),gij)
    G::Vector{Tuple{Tuple{Int64,Int64},Int64}} = calcGainVector(distancier)

    # tri décroissant par gij
    G = sortVector(G)
    
    # Initialisation des tournées élémentaires pour la fusion
    elementaryCycles::Vector{Tuple{Vector{Int64},Int64}}= initialiserChemins(dmd)

    # tournées fusionnées
    fusionnedCycles::Vector{Tuple{Vector{Int64},Int64}}= getAllFusionnedCycles(G,elementaryCycles, dmd,capa)

    # formattage des tournées
    formattedCycles::Vector{Vector{Int64}} = formatCycles(fusionnedCycles)

    # affichages finaux tournée+distance totale
    sumtot::Int64 = 0
    distcycle::Int64 = 0

    println("Les tournées retenues (heuristique Clark & Wright) sont :")

    for cycle in formattedCycles
        distcycle = calcLengthOfCycle(cycle,distancier)
        sumtot = sumtot + distcycle
        println("* ", cycle, ", longueur = ", distcycle)
    end

    println("Longueur totale : ", sumtot)
end
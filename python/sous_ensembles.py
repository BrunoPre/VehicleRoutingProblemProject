import numpy as np
import copy
class chemins:
    #constructeur
    def __init__(self, cout, endroits):
        self.cout = cout
        self.endroits=endroits
    #ajoute une element dans le chemin
    def ajoute(self, element):
        self.endroits+=element.endroits
        self.cout+=element.cout
#print method
    def print(self):
        for e in self.endroits:
            print(e)
        print("couts:"+str(self.cout))

#fonction pour sous-ensembles
def construitSous(S, cout):
    P=[chemins(0,[1])]
    for p in P:
        for s in S:
            if(p.cout+s.cout<=cout and s.endroits[len(s.endroits)-1] > p.endroits[len(p.endroits)-1]):
                toAdd= copy.deepcopy(p)
                toAdd.ajoute(s)
                P.append(toAdd)
    return P
#Main fonction
if __name__ == "__main__":
    #C= [[0,334,262,248,277,302],
     #   [334,0,118,103,551,105],
      #  [262,118,0,31,517,180],
      #  [248,103,31,0,495,152],
       # [277,551,517,495,0,476],
        #[302,105,180,152,476,0]]
    D=[0,2,4,2,4,2]
    S=[]
    #construction des objets
    index=1
    for d in D:
        toAdd=chemins(d,[index])
        S.append(toAdd)
        index+=1
    #appliquer la fonction
    P=construitSous(S, 10)
    for p in P :
        p.print()





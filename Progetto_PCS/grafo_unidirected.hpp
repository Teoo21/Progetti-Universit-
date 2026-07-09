//Si è considerato che nell'insieme di archi compare solo uno dei due archi 
//come orientazione, ma è percorribile in entrambi i versi

#include <iostream>
#include <list>
#include <set>
#include <map>

using namespace std;
class unidirected_edge {
    int nodo1,nodo2;
    public:
    unidirected_edge() {}
    unidirected_edge(int nodoa,int nodob){
        
        if (nodoa>nodob){
            nodo2=nodoa;
            nodo1=nodob;
        }
        else {
            nodo1=nodoa;
            nodo2=nodob;
        }
    }

    int from () const{
        return nodo1;
    }
    int to() const{
        return nodo2;
    }
    
    bool operator<(const unidirected_edge & arco) const{
        if(nodo1!=arco.nodo1){
            return nodo1<arco.nodo1;
        }
        return nodo2<arco.nodo2;
    }
    bool operator==(const unidirected_edge & arco) const{
        
        return (nodo1==arco.nodo1 && nodo2==arco.nodo2);
    }
    
    //Usiamo friend per poter accedere ai membri privati della classe unidirected_edge e 
    //stampare i nodi dell'arco in modo leggibile. Ci siamo aiutati con l'IA per capire come utilizzare friend
    friend ostream& operator<<(ostream& out,const unidirected_edge& e){
        out<<"("<<e.nodo1<<" , "<<e.nodo2<<")";
        return out;
    }
    
};

class undirected_graph {
    
    public:

set<unidirected_edge> archi;
set<int> nodi;
map<int, set<int>> vicini;

undirected_graph() {}

undirected_graph(const undirected_graph& grafo) {
archi = grafo.archi;
nodi = grafo.nodi;
vicini = grafo.vicini;
}

set<int> neighhours(int nodo) {
return vicini[nodo];
}

void add_edge(int a, int b) {
unidirected_edge nuovo_arco(a, b);
archi.insert(nuovo_arco);
nodi.insert(a);
nodi.insert(b);
vicini[a].insert(b);
vicini[b].insert(a);
}

set<unidirected_edge> all_edges()  {
return archi;
}

set<int> all_nodes()  {
return nodi;
}

};  
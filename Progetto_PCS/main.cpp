#include "grafo_unidirected.hpp"
#include "de_pina.hpp"
#include "gradiente_coniugato.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "Eigen/Dense"
#include "Eigen/SVD"

using namespace std;

//Le variabili vengono passate alle varie funzioni fuori dal main by reference

//La funzione apre e legge il file in input. Ogni riga viene divisa in 4 parti
//(tipo, intensità, nodo di partenza, nodo di arrivo) che vengono suddivise in varie
//strutture utili alla manipolazione.
void lettura_file(const string& filename, undirected_graph& grafo, map<unidirected_edge,
    double>& peso_arco, map<unidirected_edge, string>& tipo_arco) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Impossibile aprire il file: " << filename << endl;
        exit (0);
    }else{
    string type;
    double intensità;
    int partenza;
    int arrivo;

    while(file>>type>>intensità>>partenza>>arrivo){

        grafo.add_edge(partenza,arrivo);
        peso_arco[unidirected_edge(partenza,arrivo)] = intensità;
       
        if(type[0]=='R'||type[0]=='r'){
            tipo_arco[unidirected_edge(partenza,arrivo)]="resistenza";
        }
        else if(type[0]=='V'||type[0]=='v'){
            tipo_arco[unidirected_edge(partenza,arrivo)]="generatore";
            if(partenza>arrivo){
                tipo_arco[unidirected_edge(partenza,arrivo)]="generatore_invertito";
            }
        }
        else{
            cerr<<"Tipo di componente non riconosciuto: "<<type[0]<<endl;
            exit(0);
        }

        cout << "Tipo: " << type << ", Intensità: " << intensità << ", Partenza: " << partenza << ", Arrivo: " << arrivo << endl;
       
    }
}
file.close();

}

//Crea il set di archi a cui è associata l'etichetta "resistenza"
void ricevi_resistori(const undirected_graph& grafo, const map<unidirected_edge, string>& tipo_arco, set<unidirected_edge>& resistori) {
    for (const auto& arco : grafo.archi) {
        if (tipo_arco.at(arco) == "resistenza") {
            resistori.insert(arco);
        }
    }
}

//Crea il set di archi a cui è associata l'etichetta "generatore" o "generatore_invertito"
void ricevi_generatori(const undirected_graph& grafo, const map<unidirected_edge, string>& tipo_arco, set<unidirected_edge>& generatori) {
    for (const auto& arco : grafo.archi) {
        if (tipo_arco.at(arco) == "generatore" || tipo_arco.at(arco) == "generatore_invertito") {
            generatori.insert(arco);
        }
    }
}

//Scorriamo gli elementi di resistori e inseriamo i loro valori nella matrice R,
//che è una matrice diagonale con i valori delle resistenze sui suoi elementi diagonali
void riempi_R(Eigen::MatrixXd& R, const set<unidirected_edge>& resistori, const map<unidirected_edge, double>& peso_arco) {
    int i = 0;
    for (const auto& arco : resistori) {
        R(i, i) = peso_arco.at(arco);
        i++;
    }
}

//Se il resistore appartiene ad una maglia (ciclo minimo), inseriamo 1 o -1 a seconda della direzione del resistore
//rispetto alla maglia, considerando positivo il verso dal nodo più piccolo a quello più grande,
//altrimenti 0 se non appartiene al ciclo
void riempi_B(Eigen::MatrixXd& B, const set<unidirected_edge>& resistori, map<int, vector<int>>& registro_cicli) {
   
    for (size_t i = 0; i < registro_cicli.size(); i++) {
        const auto& ciclo = registro_cicli[i];
       
        for (size_t nodo = 0; nodo < ciclo.size() - 1; nodo++) {
            int u = ciclo[nodo];    
            int v = ciclo[nodo + 1];
            int riga_resistore = 0;
           
            for (const auto& arco : resistori) {
               
                if ((arco.from() == u && arco.to() == v) || (arco.from() == v && arco.to() == u)) {
                   
                   
                    if (u < v) {
                        B(riga_resistore, i) = 1.0;
                    }
                   
                    else {
                        B(riga_resistore, i) = -1.0;
                    }
                   
                   
                    break;
                }
               
                riga_resistore++;
            }
        }
    }
}

//per ogni generatore, se appartiene alla maglia (ciclo minimo), inseriamo il suo peso nella posizione corrispondente
//del vettore v, con segno positivo o negativo a seconda che il generatore sia normale o inverito rispetto alla maglia,
//considerando positivo il verso dal nodo più piccolo a quello più grande
void riempi_v(Eigen::VectorXd& v, const set<unidirected_edge>& generatori,
    const map<unidirected_edge, double>& peso_arco, const map<unidirected_edge,
    string>& tipo_arco, vector<vector<int>>& base_cicli) {
    int i = 0;
    for (const auto& arco : generatori) {
        for (size_t j = 0; j < base_cicli.size(); j++) {
            const auto& ciclo = base_cicli[j];
            for (size_t nodo = 0; nodo < ciclo.size() - 1; nodo++) {
                int u = ciclo[nodo];    
                int w = ciclo[nodo + 1];

                if ((arco.from() == u && arco.to() == w) || (arco.from() == w && arco.to() == u)){

                    if(u < w){
                        if (tipo_arco.at(arco) == "generatore") {
                            v(j) = v(j) - peso_arco.at(arco);
                        } else if (tipo_arco.at(arco) == "generatore_invertito") {
                            v(j) = v(j) + peso_arco.at(arco);
                        }
                       
                        break;
                }
                    if(u > w) {

                        if (tipo_arco.at(arco) == "generatore") {
                            v(j) = v(j) + peso_arco.at(arco);
                        } else if (tipo_arco.at(arco) == "generatore_invertito") {
                            v(j) = v(j) - peso_arco.at(arco);
                        }
                       
                        break;
                    }
                     
                }
            }
        }
    }
}
// funzione per la dfs ricorsiva
void dfs_ricorsiva(int u, int p, undirected_graph& grafo, map<int, bool>& visitato,
                   map<int, int>& padre, vector<vector<int>>& base_cicli,
                   set<unidirected_edge>& archiindietro_visti) {
    visitato[u] = true;
    padre[u] = p;

    for (int v : grafo.neighhours(u)) {
        if (v == p) {
            continue;
        }

        if (visitato[v]) {
            unidirected_edge arco(u, v);
           
           
            if (archiindietro_visti.find(arco) == archiindietro_visti.end()) {
               
               
                archiindietro_visti.insert(arco);
               
                vector<int> ciclo;
                ciclo.push_back(v);
                int curr = u;
                while (curr != v && curr != -1) {
                    ciclo.push_back(curr);
                    curr = padre[curr];
                }
                ciclo.push_back(v);
               
                base_cicli.push_back(ciclo);
            }
        }
        else {
           
            dfs_ricorsiva(v, u, grafo, visitato, padre, base_cicli, archiindietro_visti);
        }
    }
}

// lo usiamo nella dfs ricorsiva per trovare i cicli
vector<vector<int>> trova_cicli_dfs(undirected_graph& grafo) {
    vector<vector<int>> base_cicli;
    map<int, bool> visitato;
    map<int, int> padre;
    set<unidirected_edge> archiindietro_visti;

    for (int nodo : grafo.all_nodes()) {
        visitato[nodo] = false;
        padre[nodo] = -1;
    }

    for (int nodo : grafo.all_nodes()) {
        if (!visitato[nodo]) {
            dfs_ricorsiva(nodo, -1, grafo, visitato, padre, base_cicli, archiindietro_visti);
        }
    }
   
    return base_cicli;
}

int main() {
    string filename;
    cout << "Inserire nome del file da leggere: "<< endl;
    cin>>filename;
   
    undirected_graph circuito;
    map<unidirected_edge, double> peso_arco;
    map<unidirected_edge, string > tipo_arco;

    lettura_file(filename, circuito, peso_arco, tipo_arco);
    circuito.all_edges();
    circuito.all_nodes();

    cout << "" << endl;
    cout<<"mettere 1 per dfs e mettere 2 per de pina"<<endl;
    int scelta;
    cin >> scelta;

    vector<vector<int>> base_cicli;

    if (scelta == 1) {
        base_cicli = trova_cicli_dfs(circuito);
    } else {
        DePina solutore;
        base_cicli = solutore.calcola(circuito, peso_arco);
    }

    cout << endl;
    cout << endl;
    cout << "Le maglie sono: " << endl;
   
    map<int, vector<int>> registro_cicli;; //Dato il ciclo i-esimo associa un vettore con i nodi che lo compongono
    for (size_t i = 0; i < base_cicli.size(); ++i) {
        cout << "C" << i + 1 << " = (";
        for (size_t j = 0; j < base_cicli[i].size(); ++j) {
            registro_cicli[i].push_back(base_cicli[i][j]);
            cout<<registro_cicli[i][j];
            if (j < base_cicli[i].size() - 1) {
                cout << ", ";
            }
        }
        cout << ")" << endl;
    }

    set<unidirected_edge> resistori;
    ricevi_resistori(circuito, tipo_arco, resistori);

    set<unidirected_edge> generatori;
    ricevi_generatori(circuito, tipo_arco, generatori);
   
    Eigen::MatrixXd R(resistori.size(), resistori.size());
    R = Eigen::MatrixXd::Zero(resistori.size(), resistori.size());
    riempi_R(R, resistori, peso_arco);

    Eigen::MatrixXd B(resistori.size(), (circuito.archi.size()-circuito.nodi.size()+1));
    B = Eigen::MatrixXd::Zero(resistori.size(), (circuito.archi.size()-circuito.nodi.size()+1));
    riempi_B(B, resistori, registro_cicli);

    cout << endl;
    cout << endl;

    cout<< "B = "<<B<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"R = "<<R<<endl;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(base_cicli.size());
    riempi_v(v, generatori, peso_arco, tipo_arco, base_cicli);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(base_cicli.size(), base_cicli.size());
    A = B.transpose() * R * B;
    cout<<endl;
    cout<<endl;

    cout<<"v = "<<v<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"A = "<<A<<endl;


    gradiente_coniugato oggetto(A, v, 10e-12);
    Eigen::VectorXd I = Eigen::VectorXd::Zero(base_cicli.size());
    I = oggetto.conti();

    cout<<endl;
    cout<<endl;
    cout<<"I = "<<I<<endl;

    Eigen::VectorXd Vr = R * B * I;
    cout<<endl;
    cout<<endl;

    int i = 0;
    for(auto arco : resistori){
       
        cout<<"La caduta di potenziale sul resistore "<<arco<<" è "<< "V = "<<Vr(i)<<endl;
        i++;
    }

    return 0;
}
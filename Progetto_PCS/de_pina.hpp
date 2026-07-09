#include <vector>
#include <set>
#include <map>
#include <queue>
#include <algorithm>
#include <utility>

using namespace std;

class DePina {
   //Mi dice dati 2 insiemi di archi se essi condividono un numero di archi pari se
   //ritorna 0 e dispari se ritorna 1
    int prodotto_scalare(const set<unidirected_edge>& ciclo, const set<unidirected_edge>& S) {
        int comuni = 0;
        for (const auto& e : ciclo) {
            // entro se e è contenuto in S
            if (S.count(e)){   
            comuni=comuni+1;   
            }
        }
        return comuni % 2;
    }

    set<unidirected_edge> trova_ciclo_minimo(undirected_graph& grafo,
    const map<unidirected_edge, double>& peso, const set<unidirected_edge>& S) {
        
        //Metto un costo enorme cosi sono sicuro che non sarà il minore
        double min_costo = 1e9;
        //Creo l'insieme che conterrà gli archi di un ciclo ottimo per ognuno di quelli presenti
        set<unidirected_edge> ciclo_ottimo;
        
        for (const auto& e_test : grafo.all_edges()) {
            //Per ogni arco nel grafo prendo i suoi nodi e creo dist e padre
            int u = e_test.from(); 
            int v = e_test.to();
            
            // Per il grafo sdoppiato la mappa collega uno stato {nodo, parità} alla sua distanza o al suo padre
            map<pair<int, int>, double> dist;
            map<pair<int, int>, pair<pair<int, int>, unidirected_edge>> padre;
            
            //Questa riga è stata inizializzata con l'ausilio di IA 
            priority_queue<pair<double, pair<int, int>>, vector<pair<double, pair<int, int>>>, greater<pair<double, pair<int, int>>>> pq; //greater permette di considerare i pesi minori come prioritari
            
            //Operazioni di inizializzazione
            for (int n : grafo.all_nodes()){
                dist[{n, 0}] = 1e9;
                dist[{n, 1}] = 1e9;
            }
            
            dist[{u, 0}] = 0.0;
            pq.push({0.0, {u, 0}});

            while (!pq.empty()) {
                auto [d, stato_curr] = pq.top(); //Mi restituisce il nodo e la sua distanza che è la minore
                pq.pop(); //E poi qua li tolgo
                
                auto [curr, par] = stato_curr;
                if (d > dist[{curr, par}]) {
                    continue;
                }

                for (int vicino : grafo.neighhours(curr)) {
                    unidirected_edge e_vicino(curr, vicino); //Arco che contiene quello vicino
                    if (e_vicino == e_test){
                     continue;
                    }
                    
                    int in_S = S.count(e_vicino) ? 1 : 0;
                    int prox_par = par ^ in_S; // Cambio di parità logica tramite XOR, abbiamo cercato come farla su internet

                    if (dist[{curr, par}] + 1.0 < dist[{vicino, prox_par}]) {
                        dist[{vicino, prox_par}] = dist[{curr, par}] + 1.0;
                        padre[{vicino, prox_par}] = {{curr, par}, e_vicino};
                        pq.push({dist[{vicino, prox_par}], {vicino, prox_par}});
                    }
                }
            }
            
            int delta = S.count(e_test) ? 1 : 0;
            int par_obiettivo = 1 - delta;

            if (dist[{v, par_obiettivo}] < 1e9) {
                set<unidirected_edge> ciclo_corr = {e_test};
                pair<int, int> stato_corr = {v, par_obiettivo};
                
                double costo = dist[{v, par_obiettivo}] + 1.0;
                
                //Questo ciclo while è stato implementato con l'ausilio di IA per adattare la ricerca del minimo alla nostra 
                //classe di grafi non orientati, infatti per trovare il nodo precedente al nodo corrente, 
                //devo controllare se l'arco padre va da curr a vicino o da vicino a curr, e in base a questo aggiornare curr
                while (stato_corr != make_pair(u, 0)) {
                    auto [stato_prec, e_p] = padre[stato_corr];
                    ciclo_corr.insert(e_p);
                    stato_corr = stato_prec;
                }
                
                if (costo < min_costo) {
                    min_costo = costo;
                    ciclo_ottimo = ciclo_corr;
                }
            }
        }
        return ciclo_ottimo;
    }

public:
    vector<vector<int>> calcola(undirected_graph& grafo, const map<unidirected_edge, double>& peso) {
        vector<vector<int>> risultato_nodi;
        set<int> visitati; 
        set<unidirected_edge> albero;
        queue<int> q;
        
        if (grafo.all_nodes().empty()){
         return risultato_nodi;
        }

        for (int radice : grafo.all_nodes()) {
            if (visitati.find(radice) == visitati.end()) { 
                q.push(radice);
                visitati.insert(radice);
                while (!q.empty()) {
                    int u = q.front();
                    q.pop();
                    for (int v : grafo.neighhours(u)) {
                        if (visitati.find(v) == visitati.end()) {
                            visitati.insert(v); albero.insert(unidirected_edge(u, v)); q.push(v);
                        }
                    }
                }
            }
        }

        //Siccome non avevamo implementato precedentemente l'operatore '-' per grafi non orientati, per trovare il co-albero 
        //dobbiamo scorrere tutti gli archi del grafo e inserire quelli che non sono presenti nell'albero 
        //in un nuovo vettore di archi, che chiameremo co_albero
        vector<unidirected_edge> co_albero;
        for (const auto& e : grafo.all_edges()) {
            if (albero.find(e) == albero.end()){
             co_albero.push_back(e);
            }
        }

        vector<set<unidirected_edge>> S(co_albero.size());
        for (size_t i = 0; i < co_albero.size(); ++i){
         S[i].insert(co_albero[i]);
        }

        for (size_t i = 0; i < co_albero.size(); ++i) {
            set<unidirected_edge> C_i = trova_ciclo_minimo(grafo, peso, S[i]);
            if (C_i.empty()){
             continue;
            }

            //Ora devo convertire l'insieme di archi in una sequenza ordinata di nodi
            map<int, vector<int>> adj;
            for (const auto& e : C_i) {
                adj[e.from()].push_back(e.to());
                adj[e.to()].push_back(e.from());
            }
            //Da qui fino al prossimo commento siamo stati aiutati dall'IA per capire come convertire 
            //l'insieme di archi in una sequenza ordinata di nodi
            vector<int> nodi_ciclo = { C_i.begin()->from() };
            int curr = nodi_ciclo[0], pred = -1;
            do {
                int prox = (adj[curr][0] != pred) ? adj[curr][0] : (adj[curr].size() > 1 ? adj[curr][1] : -1);
                if (prox == -1) break;
                nodi_ciclo.push_back(prox); pred = curr; curr = prox;
            } while (curr != nodi_ciclo[0] && nodi_ciclo.size() <= C_i.size());
            risultato_nodi.push_back(nodi_ciclo);

            //Modifica dei vettori ortogonali S: scorriamo il co_albero attuale e per ogni ciclo minimo trovato, 
            //se il ciclo condivide un numero dispari di archi con S[j], allora facciamo l'operazione di XOR tra S[j] e S[i]
            for (size_t j = i + 1; j < co_albero.size(); j++) {
                if (prodotto_scalare(C_i, S[j]) == 1) {
                    for (const auto& e : S[i]) {
                        if (S[j].count(e)) S[j].erase(e); else S[j].insert(e);
                    }
                }
            }
        }
        return risultato_nodi;
    }
};

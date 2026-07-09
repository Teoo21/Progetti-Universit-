import numpy as np
import scipy as sp



class MultipleFisherDiscriminantAnalysis:
    def __init__(self, n_dimensions=None, delta_reg=1e-7):
        self.n_dimensions_ = n_dimensions
        self.reg_ = None
        self.delta_reg = delta_reg
        self.within_ = None
        self.between_ = None
        self.eigenvectors_ = None  # Matrice con colonne che sono gli autovettori
        self.eigenvalues_ = None  # i-esimo autovalore corrispondente all'i-esima colonna/autovettore in self.eigenvectors_

    def fit(self, X, y):
        """
        param X: numpy array (matrice) N-by-n di N istanze descritte da n features
        param y: numpy array (vettore) di N elementi tale che y[i] indica la classe di X[i, :]
        """
        N, n = X.shape
        classes = np.unique(y)

        if self.n_dimensions_ is None:
            self.n_dimensions_ = min(len(classes) - 1, n)
        elif self.n_dimensions_ > min(len(classes) - 1, n):
            raise ValueError("n_dimensions_ must be <= min(n_classes - 1, n_features).")
        
        # Inizializzazione matrice avente per righe i centroidi delle classi
        M = np.zeros((np.unique(y).size, n))
        # Inizializzazione lista contenente cardinalità classi
        Ns = np.zeros(np.unique(y).size)
        
        # Vettore medio del dataset
        m = X.mean(axis=0)

        # Inizializzazione matrice "within class scatter matrix"
        Sw = np.zeros((n, n))
        # Inizializzazione matrice avente per righe i centroidi delle classi
        M = np.zeros((len(classes), n))
        # Inizializzazione lista contenente cardinalità classi
        Ns = np.zeros(len(classes))

        for i, c in enumerate(classes):
            Xi = X[y == c, :]
            Ni = Xi.shape[0]
            mi = Xi.mean(axis=0)

            Xi_ = Xi - mi
            Si = Xi_.T @ Xi_

            Sw += Si
            M[i, :] = mi
            Ns[i] = Ni

        # Calcolo della matrice "between class scatter matrix"
        M_ = (M - m) * np.sqrt(np.expand_dims(Ns, axis=1))
        Sb = M_.T @ M_

        # Salvataggio come attributi
        self.within_ = Sw
        self.between_ = Sb

        if np.linalg.matrix_rank(self.within_) ==  n:  # Sw Invertibile, quindi Definita Positiva
            Sw_ = self.within_
        else:  # Sw non invertibile, quindi Sw Semi-Definita Positiva con necessità di regolarizzazione
            # k=1  un solo autovalore
            # which='SA' → Smallest Algebraic
            lambda_min = sp.sparse.linalg.eigsh(self.within_,
                                                k=1,  # un solo autovalore
                                                which='SA',  # Smallest Algebraic
                                                return_eigenvectors=False
                                                )[0]
            self.reg_ = - lambda_min + self.delta_reg
            Sw_ = self.within_ + self.reg_ * np.eye(n)
            
        S_ = np.linalg.solve(Sw_, self.between_)
            
        # Calcolo autovalori e autovettori
        self.eigenvalues_ , self.eigenvectors_ = np.linalg.eig(S_)
            
        # Selezione degli n_dimensions autovalori con val. assoluti maggiori (conserviamo solo quelli ed i risp. autovettori)
        eigen_ii = np.argsort(np.abs(self.eigenvalues_))
        eigen_ii = eigen_ii[-1::-1]

        k = self.n_dimensions_
        self.eigenvalues_ = self.eigenvalues_[eigen_ii[:k]]
        self.eigenvectors_ = self.eigenvectors_[:, eigen_ii[:k]]
            
    def transform(self, X):
        # Trasformo i dati rispetto agli autovettori calcolati con il metodo fit
        return X @ self.eigenvectors_
    

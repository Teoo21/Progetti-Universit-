#include <iostream>
#include "Eigen/Dense"
#include "Eigen/SVD"
using namespace std;

class gradiente_coniugato{
     
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd r;
    Eigen::VectorXd p;
    
    double tol;
    
    public:
    
    gradiente_coniugato(Eigen::MatrixXd& B,Eigen::VectorXd& vettore,
    double tollerenza){
      A=B;
      b=vettore;
      tol=tollerenza;
      
      //Controlliamo se la amtrice è simmetrica valutando se A è simile alla sua trasposta
      if(!A.isApprox(A.transpose())){
        cout<<"c'è un errore poichè non è simmetrica"<<endl;
        exit(1);


      }
    }

    Eigen::VectorXd conti(){
    int n = b.size(); //Prendo la dimensione del vettore b che è anche la dimensione di x
    int cont=0;
    Eigen::VectorXd x=Eigen::VectorXd::Zero(n);
    Eigen::VectorXd r0=b-A*x;
    Eigen::VectorXd p0=r0;
    double alfa;
    double beta;
    double norma=r0.norm();
    while (norma>abs(tol)){
        cont=cont+1;
        alfa=((p0.transpose()*r0)/(p0.transpose()*A*p0)).value();
        x=x+alfa*p0;
        r0=b-A*x;
        norma=r0.norm();
     
        beta=((p0.transpose()*A*r0)/(p0.transpose()*A*p0)).value();
        p0=r0-beta*p0;

        //Controlliamo se il metodo converge, se non è così è perchè la matrice è mal condizionata
        if(cont > n*10){
          cout<<"il metodo non converge"<<endl;
          exit(1);
        }
    };

    cout<<endl;
    cout<<endl;
    cout<<"Il numero di iterazioni del metodo del gradiente coniugato è "<<+cont<<endl;

    return x;

    }    
};


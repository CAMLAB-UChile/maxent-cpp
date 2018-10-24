/********************************************************************************
 *                                 Maxent                                       *
 *------------------------------------------------------------------------------*
 *     C++ class that implements the maximum-entropy basis functions            *
 *------------------------------------------------------------------------------*
 *  Version    : 1.0                                                            *
 *  Date       : 24-OCT-2018                                                  *
 *  Source code: http://camlab.cl/software/maxent                               *
 *  Author     : R. Silva-Valenzuela, MSc student, rsilvavalenzue@ing.uchile.cl *
 *  Supervisor : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro     *
 *                                                                              *
 *            (See Copyright and License notice in "license.txt")               *
 *            (See updates and version details in "version.txt")                *
 *------------------------------------------------------------------------------*
 *                                                                              *
 * References                                                                   *
 * =============                                                                *
 * [1] Sukumar N.(2008), Maxent-F90: Fortran 90 Library for Maximum-Entropy     *
 *     Basis Function. User's Reference Manual Version 1.4. UC Davis            *
 *                                                                              *
 * [2] A. Ortiz-Bernardin. Elementos Finitos Generalizados Universidad de       *
 *     Chile. Material Docente. Semestre Otoño 2016                             *
 *                                                                              *
 * [3] A. Ortiz-Bernardin.(2011), Maximum-Entropy Meshfree Method for Linear    *
 *     and Nonlinear Elasticity. PhD thesis, University of California, Davis.   *
 *                                                                              *
 * Adapted from A.O-B's "Maxent basis functions for MATLAB" Version 3.4:        *
 * http://camlab.cl/software/other-software/                                    *
 *                                                                              *
 *******************************************************************************/

#ifndef MAXENT_HPP_INCLUDED
#define MAXENT_HPP_INCLUDED

#include <iostream>
using namespace std;

#include<cstdlib> //para exit

#include <vector>
#include "Nodes.hpp"

#include <iomanip> //Dar formato a nœmeros para la salida del programa

#include <Eigen/Dense>
using namespace Eigen;

class Maxent
{
private:

  int dim; //spatial dimension
  string scheme; // solution scheme "newton"
  string prior; // prior in the weight function either "uniform", "cubic", "quartic", "gaussian".
  double radious; //is the support size
  int MaxIter; //maximum number of iteration
  double ctol;  //convergence tolerance
  bool maxentprint; // print convergence information
  double gamma;
  double hnode;
  //variables internas
  MatrixXd WeightFuncDer;
  VectorXd WeightFunc;
  VectorXd Phi;
  MatrixXd Phider;
  MatrixXd H; //Hessian Matrix
  MatrixXd MA;
  MatrixXd A;
  MatrixXd MC;
  MatrixXd diCoord;
  int NNodesC;
  MatrixXd NodeCoord;
  VectorXd Point;
  //funciones internas
  void CheckConsistency();
  void MakePhi();
  void MakePhider();
  void PriorWeightFunction();
  void CalFlambda(VectorXd&,VectorXd&, MatrixXd&);
  void PhiderMatrix();
  void ShowBasicFunc();

public:

  Maxent(); //constructor por omision
  Maxent(int, string,  int, double, bool, double, double); //constructor
  Maxent(const Maxent&); //Constructor copia
  void operator=(const Maxent&);//Operador de asignacion
  ~Maxent(); //destructor clase
  void BasicFunc(int, vector <Node>, VectorXd);
  VectorXd GetPhi();
  MatrixXd GetPhider();

};

Maxent::Maxent(int auxdim, string auxprior, int auxMaxIter, double auxctol, bool auxmaxentprint, double auxgamma, double auxhnode)
{
  dim=auxdim;
  scheme="newton";
  prior=auxprior;
  MaxIter=auxMaxIter;
  ctol=auxctol;
  maxentprint=auxmaxentprint;
  gamma=auxgamma;
  hnode=auxhnode;

  if(auxprior=="gaussian")
    radious=hnode*sqrt(-log(auxctol)/auxgamma);
  else
    radious=auxgamma*hnode;
    
  cout <<"*************************************************************"<<endl;
  cout <<"******************** MAXENTC++ CLASS ************************"<< endl;
  cout <<"*************************************************************"<<endl<<endl;

  if(maxentprint)
  {
    cout << "Dimension: " <<auxdim << endl<<endl;
    cout << "Type Prior Weigth function: " << auxprior<< endl<<endl;
    cout << "Maximum permissible iterations: " << auxMaxIter<< endl<<endl;
    cout << "Convergence tolerance: " << auxctol<< endl<<endl;
    cout <<"*************************************************************"<<endl;
  }

}

Maxent::Maxent()
{

}

Maxent::Maxent(const Maxent& Maxentold)
{
  dim=Maxentold.dim;
  scheme=Maxentold.scheme;
  prior=Maxentold.prior;
  radious=Maxentold.radious;
  MaxIter=Maxentold.MaxIter;
  ctol=Maxentold.ctol;
  maxentprint=Maxentold.maxentprint;
  gamma=Maxentold.gamma;
  hnode=Maxentold.hnode;
}

void Maxent::operator=(const Maxent& Maxentnew)
{
  dim=Maxentnew.dim;
  scheme=Maxentnew.scheme;
  prior=Maxentnew.prior;
  radious=Maxentnew.radious;
  MaxIter=Maxentnew.MaxIter;
  ctol=Maxentnew.ctol;
  maxentprint=Maxentnew.maxentprint;
  gamma=Maxentnew.gamma;
  hnode=Maxentnew.hnode;
}

Maxent::~Maxent()
{

}

void Maxent::BasicFunc(int auxNNodesC,vector <Node> auxNodesContribute, VectorXd auxPoint)
{
  Point=auxPoint;
  NNodesC=auxNNodesC;
  if(maxentprint)
  {
    cout <<"Using Sample Point: " << setw(0)<<fixed<<setprecision(3)<<Point.transpose() << endl << endl;
    cout << "Number of nodes contributes: " << auxNNodesC << endl<<endl;
    cout <<"Node" << "            Coord"<<endl;
    for(int i=0;i<auxNNodesC;i++)
    {
     cout << " "<<auxNodesContribute[i].GetTagNode()<<"          ";
     cout <<auxNodesContribute[i].GetCoordx() << " "
          <<auxNodesContribute[i].GetCoordy() << " "
          <<auxNodesContribute[i].GetCoordz() <<endl;;
    }
  }
  else
  {
    cout<<left<<setw(6)<<fixed<<setprecision(3)<< Point(0)<<" "<<setw(6)<<Point(1)<<" "<<setw(6)<<Point(2)<<" ";
  }

  NodeCoord.resize(auxNNodesC,3); // NodeCoord = [x1 y1 z1; x2 y2 z2; ...]
  for (int i=0; i<auxNNodesC; i++)
  {
    NodeCoord(i,0)= auxNodesContribute[i].GetCoordx();
    NodeCoord(i,1)= auxNodesContribute[i].GetCoordy();
    NodeCoord(i,2)= auxNodesContribute[i].GetCoordz();
  }

  diCoord.resize(auxNNodesC,dim); // diCoord = [x1-xp y1-yp z1-zp; x2-xp y2-yp z2-zp; ...]
  for (int i = 0; i < auxNNodesC; i++)
    for (int j = 0; j < dim; j++)
      diCoord(i,j) = (NodeCoord(i,j) - Point(j));


  this->PriorWeightFunction();

  this->MakePhi();

  this->MakePhider();
 
  if(maxentprint)
  this-> ShowBasicFunc();

  this->CheckConsistency();

}

void Maxent::ShowBasicFunc()
{
  cout << "Basic Function Phi" << endl<< setw(0) << fixed << setprecision(5) <<Phi<<endl;
  cout <<endl;
  cout << "Derivative Basic Function Phider" << setw(0) << fixed << setprecision(5) <<endl<<Phider<<endl<<endl;
}

void Maxent::MakePhi()
{
  int nn=0; //number iter

  MatrixXd jacobian(dim,dim);// jacobian (Hessian matrix H in maxent context)

  VectorXd lambda(dim); // lambda, Lagrange multipliers
  for (int i=0; i<dim; i++)
  {
    lambda(i)=0.0; // initial lambda
  }

  VectorXd Flambda(dim);
  for (int i=0; i<dim; i++) // initial  Newton's residual: F(lambda)=-Sum(Phi*(x-xp),1,NNodesC)
  {
    Flambda(i)=0.1; //Flambdam= [0.1;0.1;0.1]
  }

  VectorXd dlam(dim); //declaration Increase of lambda

  if(maxentprint)
    cout << endl <<"iterating...  ";
  else
    cout <<"  iterating...  ";

  while(Flambda.norm()>ctol) // F(lambda) = Grad(log Z) = 0
  {

    CalFlambda(lambda, Flambda, jacobian); //calculate Jacobiano and Flambda for current Lagragian multipliers

    dlam = -jacobian.colPivHouseholderQr().solve(Flambda);

    lambda=lambda+dlam; //Increase of lambda

    nn=nn+1; //increasa number of iterations

    if (MaxIter<nn)
    {
      cout <<"Newton Method Failed, no convergence in " << MaxIter << " iterations" << endl;
      exit(1);
    }

  }


  if(maxentprint)
    cout << "Total iterations: " << nn << endl << endl;
  else
    cout << "total iter " <<setw(3)<<nn <<"    ";

}



void Maxent::MakePhider()
{
  this->PhiderMatrix();

  MatrixXd invH=H.inverse();
  MatrixXd MPhi(NNodesC,dim);

  for(int i=0; i<dim;i++)
  {
    for(int j=0; j<NNodesC; j++)
    {
      MPhi(j,i)=Phi(j);
    }
  }

  Phider=MPhi.cwiseProduct(diCoord*(invH-invH*A)) + MPhi.cwiseProduct(MA)-Phi*MC;

}


VectorXd Maxent::GetPhi()
{
    return Phi;
}

MatrixXd Maxent::GetPhider()
{
    return Phider;
}


void Maxent::CalFlambda(VectorXd &lambda, VectorXd &Flambda, MatrixXd &jacobian)
{
  VectorXd argu = -(diCoord * lambda);// argument partition function components for current lambda

  VectorXd zi(NNodesC);
  for (int i = 0; i < NNodesC; i++)
    zi(i) = WeightFunc(i) * exp(argu(i));

  double Z=0;
  for(int i=0;i<NNodesC;i++)
  {
    Z=Z+zi(i);  // partition function
  }

  Phi = zi / Z;  // maxent basis functions for current lambda

  for (int i=0; i<dim; i++) // Restart Newton's residual F(lambda)
  {
    Flambda(i)=0;
  }


  for (int i = 0; i < dim; i++)  // Newton's residual:   F(lambda)=-Sum(Phi*(x-xp),1,NNodesC)
  {
    for(int j=0; j<NNodesC; j++)
    {
      Flambda(i) = Flambda(i)  - (diCoord(j,i) * Phi(j));
    }
  }


  double Sum=0;
  for(int i=0; i<dim; i++)
  {
    for(int j=0; j<dim;j++)
    {
      Sum=0;
      for(int k=0;k<NNodesC;k++)
      {
        Sum=Sum+Phi(k)*diCoord(k,i)*diCoord(k,j);
      }
      jacobian(i,j)= Sum - Flambda(i)*Flambda(j); //Jacobian (Hessian matrix H in maxent context)
    }
  }

}

void Maxent::PhiderMatrix()
{
  H.resize(dim,dim);
  for (int i=0;i<dim;i++)
  {
    for(int j=0;j<dim;j++)
    {
      H(i,j)=Phi.dot(diCoord.col(i).cwiseProduct(diCoord.col(j)));
    }
  }

  MA.resize(NNodesC,dim);
  for(int i=0; i<NNodesC;i++)
  {
    for(int j=0;j<dim;j++)
    {
      MA(i,j)=(1/WeightFunc(i))*WeightFuncDer(i,j);
    }
  }

  A.resize(dim,dim);
  for (int i=0;i<dim;i++)
  {
    for(int j=0;j<dim;j++)
    {
      A(i,j)=(Phi.cwiseProduct(MA.col(j))).dot( diCoord.col(i) );
    }
  }

  MC.resize(1,dim);
  for(int i=0; i<dim;i++)
  {
    MC(0,i)=Phi.dot(MA.col(i));
  }

}

void Maxent::PriorWeightFunction()
{
  WeightFunc.resize(NNodesC);
  WeightFuncDer.resize(NNodesC,dim);

  VectorXd di(NNodesC);
  for (int i=0; i<NNodesC; i++)
  {
    di(i)= (diCoord.row(i)).norm();
  }

  if(prior=="cubic")
  {
    double Iradious = 1.0/radious;

    VectorXd q = di*Iradious;

    for (int j=0; j<NNodesC;j++)
    {
      if( 0.0<=q(j) && q(j)<=0.5 )
      {
        WeightFunc(j)= (2.0/3.0) - 4.0*q(j)*q(j) + 4.0*q(j)*q(j)*q(j);
        for(int k=0; k<dim; k++)
          WeightFuncDer(j,k)=Iradious*Iradious*(8-12*q(j)) * diCoord(j,k);
      }
      else if(0.5<q(j) && q(j)<=1.0 )
      {
        WeightFunc(j)= (4.0/3.0) - 4.0*q(j) + 4*q(j)*q(j) - (4.0/3.0)*q(j)*q(j)*q(j);
        for(int k=0; k<dim; k++)
          WeightFuncDer(j,k)=Iradious*Iradious*(4.0/q(j) - 8.0 + 4.0*q(j)) * diCoord(j,k);
      }
      else
      {
        cout << "Fatal error!, error calculus cubic PriorWeightFunction " << endl;
        exit(1);
      }
    }
  }
  else if (prior =="uniform")
  {
    for(int i=0;i<NNodesC;i++)
      WeightFunc(i)=1.0;

    for(int j=0;j<NNodesC;j++)
    {
      for(int k=0;k<dim;k++)
      WeightFuncDer(j,k)=0.0;
    }

  }
  else if(prior=="gaussian")
  {
    VectorXd beta(NNodesC);
    for(int i=0; i<NNodesC; i++)
      beta(i)=gamma/(hnode*hnode);

    VectorXd argu(NNodesC);
    for(int i=0;i<NNodesC;i++)
      argu(i)=-beta(i)*di(i)*di(i);

    for (int j=0; j<NNodesC;j++)
      WeightFunc(j)=exp(argu(j));

    for(int i=0;i<dim;i++)
    {
      for(int j=0; j<NNodesC;j++)
      {
        WeightFuncDer(j,i)=2*WeightFunc(j)*diCoord(j,i)*beta(j);
      }
    }

  }
  else if (prior=="quartic")
  {
    double Iradious = 1.0/radious;

    VectorXd q = di*Iradious;

    for (int j=0; j<NNodesC;j++)
    {
      if( 0.0<=q(j) && q(j)<=1.0 )
      {
        WeightFunc(j)= 1.0 - 6.0*q(j)*q(j) + 8.0*q(j)*q(j)*q(j) - 3.0*q(j)*q(j)*q(j)*q(j);
        for(int k=0; k<dim; k++)
          WeightFuncDer(j,k)= ( 12.0*Iradious*Iradious - 24.0*q(j)*Iradious*Iradious + 12.0*q(j)*q(j)*Iradious*Iradious ) * diCoord(j,k);
      }
      else
      {
        cout << "Fatal error!, error calculus quartic PriorWeightFunction " << endl;
        exit(1);
      }
    }

  }

}

void Maxent::CheckConsistency()
{
  int num1 = -log10(ctol);
  int num2 = double(num1)/2.0;
  double num3 = num2;
  double tol = (pow(10.0,-num3));

  if(dim==1)
  {
    VectorXd Phider_x=Phider.col(0);
    VectorXd Coordx=NodeCoord.col(0);

    double SumPhi=0;
    for(int i=0; i<NNodesC;i++)
      SumPhi=SumPhi+Phi(i);

    double SumPhi_xi=Phi.dot(Coordx);
    double SumPhider_x_xi=Phider_x.dot(Coordx);

    if(abs(SumPhi-1.0)<tol &&
      abs(SumPhi_xi-Point(0))<tol &&
      abs(SumPhider_x_xi-1.0)<tol)
    {
        cout << "Consistency check... ok" << endl;
    }
    else
    {
        cout << "Consistency check... Failed" << endl;
    }
  }
  else if (dim==2)
  {
    VectorXd Phider_x=Phider.col(0);
    VectorXd Phider_y=Phider.col(1);
    VectorXd Coordx=NodeCoord.col(0);
    VectorXd Coordy=NodeCoord.col(1);

    double SumPhi=0;
    for(int i=0; i<NNodesC;i++)
      SumPhi=SumPhi+Phi(i);

    double SumPhi_xi=Phi.dot(Coordx);
    double SumPhi_yi=Phi.dot(Coordy);

    double SumPhider_x_xi=Phider_x.dot(Coordx);
    double SumPhider_y_xi=Phider_y.dot(Coordx);

    double SumPhider_x_yi=Phider_x.dot(Coordy);
    double SumPhider_y_yi=Phider_y.dot(Coordy);

    if(abs(SumPhi-1.0)<tol &&
      abs(SumPhi_xi-Point(0))<tol &&
      abs(SumPhi_yi-Point(1))<tol &&
      abs(SumPhider_x_xi-1.0)<tol &&
      abs(SumPhider_y_xi-0.0)<tol &&
      abs(SumPhider_x_yi-0.0)<tol &&
      abs(SumPhider_y_yi-1.0)<tol)
    {
      cout << "Consistency check... ok" << endl;
    }
    else
    {
      cout << "Consistency check... Failed" << endl;
    }
  }
  else if(dim==3)
  {
    VectorXd Phider_x=Phider.col(0);
    VectorXd Phider_y=Phider.col(1);
    VectorXd Phider_z=Phider.col(2);

    VectorXd Coordx=NodeCoord.col(0);
    VectorXd Coordy=NodeCoord.col(1);
    VectorXd Coordz=NodeCoord.col(2);

    double SumPhi=0;
    for(int i=0; i<NNodesC;i++)
        SumPhi=SumPhi+Phi(i);

    double SumPhi_xi=Phi.dot(Coordx);
    double SumPhi_yi=Phi.dot(Coordy);
    double SumPhi_zi=Phi.dot(Coordz);
    double SumPhider_x_xi=Phider_x.dot(Coordx);
    double SumPhider_y_xi=Phider_y.dot(Coordx);
    double SumPhider_z_xi=Phider_z.dot(Coordx);
    double SumPhider_x_yi=Phider_x.dot(Coordy);
    double SumPhider_y_yi=Phider_y.dot(Coordy);
    double SumPhider_z_yi=Phider_z.dot(Coordy);
    double SumPhider_x_zi=Phider_x.dot(Coordz);
    double SumPhider_y_zi=Phider_y.dot(Coordz);
    double SumPhider_z_zi=Phider_z.dot(Coordz);

    if(abs(SumPhi-1.0)<tol &&
       abs(SumPhi_xi-Point(0))<tol &&
       abs(SumPhi_yi-Point(1))<tol &&
       abs(SumPhi_zi-Point(2))<tol &&
       abs(SumPhider_x_xi-1.0)<tol &&
       abs(SumPhider_y_xi-0.0)<tol &&
       abs(SumPhider_z_xi-0.0)<tol &&
       abs(SumPhider_x_yi-0.0)<tol &&
       abs(SumPhider_y_yi-1.0)<tol &&
       abs(SumPhider_z_yi-0.0)<tol &&
       abs(SumPhider_x_zi-0.0)<tol &&
       abs(SumPhider_y_zi-0.0)<tol &&
       abs(SumPhider_z_zi-1.0)<tol)

    {
      cout << "Consistency check... ok" << endl;
    }
    else
    {
      cout << "Consistency check... Failed" << endl;
    }
  }

  if(maxentprint)
  cout <<"*************************************************************"<<endl;

}

#endif // MAXENT_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;

#include "Nodes.hpp"
#include "Maxent.hpp"

int main(int argc, const char * argv[])
{
    /*  
     a 3D cloud of 6 nodes
    */    
    
    int dim=3; //spatial dimension
    
    int NNodesC = 6; //Number Nodes contribute
    
    VectorXd point(3); //coordinates of the sample points p(x,y,z)
    point(0)=0.1;
    point(1)=0.1;
    point(2)=0.1;

    vector <Node> NodesContribute(NNodesC); //Vector class Node contribute
    NodesContribute[0] = Node(1,0.0,0.0,0.0);
    NodesContribute[1] = Node(2,1.0,0.0,0.0);
    NodesContribute[2] = Node(3,0.0,1.0,0.0);
    NodesContribute[3] = Node(4,0.0,0.0,1.0);
    NodesContribute[4] = Node(5,1.0,1.0,0.0);
    NodesContribute[5] = Node(6,1.0,0.0,1.0);
    
    int MaxIter=100; //maximum number of iteration
    
    double ctol= 1E-10; //convergence tolerance
    
    double gamma=4.2426; // parameter that controls the support size of the basic function

    double hnode = 1; //characteristic nodal spacing
  
    string prior = "cubic"; // prior in the weight function either "uniform", "cubic", "quartic", "gaussian".
    
    bool maxentprint = true; // print total information

    Maxent maxent(dim, prior, MaxIter, ctol, maxentprint, gamma, hnode);
    
    maxent.BasicFunc(NNodesC,NodesContribute,point);

    VectorXd Phi = maxent.GetPhi();
    
    MatrixXd Phider =  maxent.GetPhider();
    
    return 0;
}

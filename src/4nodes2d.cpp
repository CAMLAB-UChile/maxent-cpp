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
     a 2D cloud of 4 nodes
    */  

    int dim=2; //spatial dimension
    
    int NNodesC = 4; //Number Nodes contribute
    
    VectorXd point(3); //coordinates of the sample points p(x,y,z)
    point(0)=0.85;
    point(1)=0.572;
    point(2)=0.0;
    
    vector <Node> NodesContribute(NNodesC); //Vector class Node contribute
    NodesContribute[0] = Node(1,0.0,0.0,0.0);
    NodesContribute[1] = Node(2,1.0,0.0,0.0);
    NodesContribute[2] = Node(3,1.0,1.0,0.0);
    NodesContribute[3] = Node(4,0.0,1.0,0.0);

    int MaxIter=100; //maximum number of iteration
    
    double ctol= 1E-14; //convergence tolerance
    
    double gamma=2.0;; // parameter that controls the support size of the basic function

    double hnode = 1; //characteristic nodal spacing
 
    string prior = "cubic"; // prior in the weight function either "uniform", "cubic", "quartic", "gaussian"
    
    bool maxentprint = true; // print information
    
    Maxent maxent(dim, prior, MaxIter, ctol, maxentprint, gamma, hnode);
    
    maxent.BasicFunc(NNodesC,NodesContribute,point);

    VectorXd Phi = maxent.GetPhi();
    
    MatrixXd Phider =  maxent.GetPhider();

    return 0;
}

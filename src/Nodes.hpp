#ifndef NODES_HPP_INCLUDED
#define NODES_HPP_INCLUDED

class Node
{
protected:
  int TagNode; //Tag Node
  double x,y,z; //Coordenadas del nodo

public:
  Node(); 
  Node(int, double, double, double); 
  Node(const Node&); 
  void operator=(const Node&);
  ~Node(); 
  int GetTagNode();
  double GetCoordx();
  double GetCoordy();
  double GetCoordz();
};

//implementación

Node::Node()
{
  TagNode=0;
  x=0;
  y=0;
  z=0;
}

Node::Node(int auxTagNode, double Nx, double Ny, double Nz)
{
  TagNode=auxTagNode;
  x=Nx;
  y=Ny;
  z=Nz;
}

Node::~Node()
{
    
}

Node::Node(const Node& Nodeold)
{
  TagNode=Nodeold.TagNode;
  x=Nodeold.x;
  y=Nodeold.y;
  z=Nodeold.z;
}

void Node::operator=(const Node& Nodenew)
{
  TagNode=Nodenew.TagNode;
  x=Nodenew.x;
  y=Nodenew.y;
  z=Nodenew.z;
}

int Node::GetTagNode()
{
  return TagNode;
}

double Node::GetCoordx()
{
  return x;
}

double Node::GetCoordy()
{
  return y;
}

double Node::GetCoordz()
{
  return z;
}

#endif // NODES_HPP_INCLUDED

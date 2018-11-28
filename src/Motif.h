#ifndef MOTIF_H_
#define MOTIF_H_
#include <armadillo>
using namespace arma;

class Motif{
public:
     int id;
	   uvec edges;
     int degree;

public:
	   Motif(int count): id(count){ degree = 0;}

	   ~Motif(){}
};

#endif /* MOTIF_H_ */

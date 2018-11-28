#ifndef GRAPH_H_
#define GRAPH_H_
#include "Motif.h"
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
using namespace boost;
using namespace std;

class Graph{
public:
	int NodeSize;
	int EdgeSize;
  vector<dynamic_bitset<> > adjacencyList;
	int **adjacencyMatrix;
  vec C; //capability vector of edges
	vec used;

public:
	Graph(int nodeNum, char* fileName): NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
		  adjacencyMatrix = new int*[nodeNum];
		  for(int i=0;i<nodeNum;i++)
		     adjacencyMatrix[i] = new int[nodeNum]();
		  readNetworkFile(fileName);
	}

  ~Graph(){
		for(int i=0;i<NodeSize;i++)
		delete[] adjacencyMatrix[i];
		delete[] adjacencyMatrix;
		adjacencyList.clear();
	}

  void findMotifs(vector<Motif*>& pattern, int caseN);
private:
  void readNetworkFile(char* fileName);
};

#endif

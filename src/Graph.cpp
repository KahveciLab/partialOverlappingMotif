#include "Graph.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>

void Graph::readNetworkFile(char* fileName){
	int sourceNode, destNode, mappingCounter = 0;
	ifstream in_stream;
	string sourceName, destName, str1, count;
	int edgeIndex=0;
	in_stream.open(fileName);
	vector<string> strs;
	map<string,int> nodeNameToNoMap;
	vector<int> capa;

	while(!in_stream.eof()) {
		getline(in_stream, str1);

		if(str1.size()>0){
			strs.clear();
			boost::split(strs,str1,boost::is_any_of("\t"));

			sourceName=strs[0];
			destName=strs[1];
			count = strs[2];
			stringstream ss(count);

			if(nodeNameToNoMap.find(sourceName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(sourceName, mappingCounter));
				++mappingCounter;
			}

			if(nodeNameToNoMap.find(destName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(destName, mappingCounter));
				++mappingCounter;
			}

			sourceNode = nodeNameToNoMap.find(sourceName)->second;
			destNode = nodeNameToNoMap.find(destName)->second;

			if(sourceNode != destNode && adjacencyList[sourceNode][destNode]==0){
				adjacencyMatrix[sourceNode][destNode] = edgeIndex;
				adjacencyList[sourceNode][destNode] = 1;
        int num = 0;
				ss >> num;
				capa.push_back(num);
				edgeIndex++;
			}
		}
	}

	in_stream.close();
	EdgeSize=edgeIndex;
	strs.clear();
	nodeNameToNoMap.clear();

	//initialize capability vector
	C.resize(EdgeSize);
	used=zeros(EdgeSize);

  for(int i=0;i<EdgeSize;i++)
	  C(i)=capa[i];

}


void Graph::findMotifs(vector<Motif*>& pattern, int caseN){
	dynamic_bitset<>::size_type it,itt,iit;

	int index1,index2,index3, index4;
	index1=index2=index3=index4=-1;

	int count = 0;
  if(caseN == 1){ //Feed forward loop
		for(int i=0;i<NodeSize;i++){
			it=adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				itt=adjacencyList[i].find_next(it);
				while(itt!=dynamic_bitset<>::npos){
					if(adjacencyList[it][itt]==1){
						Motif* t=new Motif(count);
						index1=adjacencyMatrix[i][it];
						index2=adjacencyMatrix[i][itt];
						index3=adjacencyMatrix[it][itt];
						t->edges << index1 << index2 << index3 << endr;
						used.elem(t->edges) += 1;
						pattern.push_back(t);
						count++;
					}
					if(adjacencyList[itt][it]==1){
						Motif* t=new Motif(count);
						index1=adjacencyMatrix[i][it];
						index2=adjacencyMatrix[i][itt];
						index3=adjacencyMatrix[itt][it];
						t->edges << index1 << index2 << index3 << endr;
						used.elem(t->edges) += 1;
						pattern.push_back(t);
						count++;
					}

					itt=adjacencyList[i].find_next(itt);
				}
				it=adjacencyList[i].find_next(it);
			}
		}
	}
	else if(caseN==2){ //bifan
		for(int i=0;i<NodeSize-1;i++)
		for(int j=i+1;j<NodeSize;j++){
			dynamic_bitset<> tdb = adjacencyList[i];
			tdb &= adjacencyList[j];
			if(tdb.count()>=2){
				it = tdb.find_first();
				while(it!=dynamic_bitset<>::npos){
					itt=tdb.find_next(it);
					while(itt!=dynamic_bitset<>::npos){
						Motif* t=new Motif(count);
						index1=adjacencyMatrix[i][it];
						index2=adjacencyMatrix[i][itt];
						index3=adjacencyMatrix[j][it];
						index4=adjacencyMatrix[j][itt];
						t->edges << index1 << index2 << index3 << index4 << endr;
						used.elem(t->edges) += 1;
						pattern.push_back(t);
						count++;
						itt = tdb.find_next(itt);
					}
					it = tdb.find_next(it);
				}
			}
		}
	}
  else if(caseN==3){ //bi-parallel
		for(int i=0;i<NodeSize;i++){
			it = adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				itt = adjacencyList[i].find_next(it);
				while(itt != dynamic_bitset<>::npos){
					dynamic_bitset<> tdb=adjacencyList[it];
					tdb &= adjacencyList[itt];
					tdb[i]=0;

					iit = tdb.find_first();
          while(iit!=dynamic_bitset<>::npos){
						Motif* t=new Motif(count);
						index1=adjacencyMatrix[i][it];
						index2=adjacencyMatrix[i][itt];
						index3=adjacencyMatrix[it][iit];
						index4=adjacencyMatrix[itt][iit];
						t->edges << index1 << index2 << index3 << index4 << endr;
						used.elem(t->edges) += 1;
						pattern.push_back(t);
						count++;

						iit = tdb.find_next(iit);
					}

					itt = adjacencyList[i].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}
	}
	else if(caseN==4){ // cascade and delay
		for(int i=0;i<NodeSize;i++){
			it = adjacencyList[i].find_next(i);
			while(it!=dynamic_bitset<>::npos){
				itt = adjacencyList[it].find_next(i);
				while(itt != dynamic_bitset<>::npos){
					if(adjacencyList[itt][i]==1){
						Motif* t=new Motif(count);
						index1=adjacencyMatrix[i][it];
						index2=adjacencyMatrix[it][itt];
						index3=adjacencyMatrix[itt][i];
						t->edges << index1 << index2 << index3 << endr;
						used.elem(t->edges) += 1;
						pattern.push_back(t);
						count++;
					}
					itt = adjacencyList[it].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}

	}






}

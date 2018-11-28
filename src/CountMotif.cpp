#include "Graph.h"
#include <iostream>
#include <random>
#include <ctime> // Needed for the true randomization
#include <cstdlib>
#include <cmath>
#include <unordered_set>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/algorithm/string.hpp>
Graph* g;

struct CompareNode : public std::binary_function<Motif*, Motif*, bool>
{
	bool operator()(const Motif* lhs, const Motif* rhs) const{
	    return lhs->degree > rhs->degree;
	}
};

typedef heap::fibonacci_heap<Motif*,heap::compare<CompareNode> >::handle_type cn;
//randomly generate a solution
double generateSolution(vector<Motif*>& pattern, dynamic_bitset<>& solution, vec& CTemp, vector<int>& nodeSet){
  solution.reset();
  std::random_device r;
  std::default_random_engine generator(r());
  std::bernoulli_distribution distribution(0.5);
  for(int i=0;i<pattern.size();i++){
    bool isEx=distribution(generator);
    if(isEx){
      uvec z = (CTemp.elem(pattern[i]->edges) < g->C.elem(pattern[i]->edges));
      bool states = all(z);
      if(states){
        solution[i]=1;
        CTemp.elem(pattern[i]->edges) +=1;
        nodeSet.push_back(i);
      }
    }
  }

}

void localSearch(vector<Motif*>& pattern, vector<dynamic_bitset<> >& madjacencyList, dynamic_bitset<>& solution,  vec& sol_C, vector<int>& nodeSet ){
  dynamic_bitset<>::size_type it, itt;
  bool updated = true;
  std::random_device r;

  while(updated){
    updated = false;

    vector<int> nodes;
    dynamic_bitset<> noNeighborNodes(pattern.size());
    while(!nodeSet.empty()){
      std::default_random_engine generator(r());
      std::uniform_int_distribution<int> distribution(0,nodeSet.size()-1);
      int index = distribution(generator);
      int removedNode = nodeSet[index];

      vector<pair<int,int>> possSwap;
      dynamic_bitset<> solu_prim= solution;
      solu_prim[removedNode]=0;
      vec CTemp = sol_C;
      CTemp.elem(pattern[removedNode]->edges) -=1;
      vec CDiff = g->C - CTemp; //how much we can add

      //find neighbor
      dynamic_bitset<> tdb(pattern.size());
      for(int i=0;i<pattern[removedNode]->edges.n_elem;i++){
        tdb |= madjacencyList[pattern[removedNode]->edges(i)];
      }
      dynamic_bitset<> neigh = tdb;
      neigh ^= solu_prim;
      tdb &= neigh;

      it= tdb.find_first();
      while(it!=dynamic_bitset<>::npos){
        itt=tdb.find_next(it);
        while(itt!=dynamic_bitset<>::npos){
          vec incre=zeros(g->EdgeSize);
          incre.elem(pattern[it]->edges) += 1;
          incre.elem(pattern[itt]->edges) += 1;
          uvec z = (incre <= CDiff);
          bool states = all(z);
          if(states){
            possSwap.push_back(pair<int,int>(it,itt));
          }
          itt=tdb.find_next(itt);
        }
        it=tdb.find_next(it);
      }

      if(!possSwap.empty()){
        std::default_random_engine generator0(r());
        std::uniform_int_distribution<int> distribution0(0,possSwap.size()-1);
        int num = distribution0(generator0);
        int x = possSwap[num].first;
        int y = possSwap[num].second;
//        cout<<"remove "<<removedNode<<"\t insert "<<x<<" and "<<y<<endl;
        sol_C.elem(pattern[x]->edges) += 1;
        sol_C.elem(pattern[y]->edges) += 1;
        sol_C.elem(pattern[removedNode]->edges) -= 1;
        solution[removedNode]=0;
        solution[x]=1;
        solution[y]=1;
        nodeSet.push_back(x);
        nodeSet.push_back(y);
        updated = true;
      }
      else{
        nodes.push_back(removedNode);
        noNeighborNodes |= tdb;
      }

      //remove removeNode from nodes set
      nodeSet[index] = nodeSet[nodeSet.size()-1];
      nodeSet.pop_back();
    }
    noNeighborNodes.flip();

    if(noNeighborNodes.any()){
      it = noNeighborNodes.find_first();
      while(it!=dynamic_bitset<>::npos){
        nodeSet.push_back(it);
        it = noNeighborNodes.find_next(it);
      }

      while(!nodeSet.empty()){
        std::default_random_engine generator(r());
        std::uniform_int_distribution<int> distribution(0,nodeSet.size()-1);
        int index = distribution(generator);
        int selectNode = nodeSet[index];

        uvec z = (sol_C.elem(pattern[selectNode]->edges) < g->C.elem(pattern[selectNode]->edges));
        bool states = all(z);
        if(states){
          solution[selectNode]=1;
          sol_C.elem(pattern[selectNode]->edges) +=1;
          updated = true;
          nodes.push_back(selectNode);
//          cout<<"select node "<<selectNode<<" from free nodes"<<endl;
        }

        //remove removeNode from nodes set
        nodeSet[index] = nodeSet[nodeSet.size()-1];
        nodeSet.pop_back();
      }
    }

    nodeSet = nodes;

  }

}

void pertub(vector<Motif*>& pattern, vector<dynamic_bitset<> >& madjacencyList, dynamic_bitset<>& solution,  vec& sol_C, vector<int>& nodeSet){
  dynamic_bitset<>::size_type it;
  std::random_device r;

  std::default_random_engine generator1(r());
  std::bernoulli_distribution distribution1(0.3);
  vector<int> nodes;

  while(!nodeSet.empty()){
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> distribution(0,nodeSet.size()-1);
    int index = distribution(generator);
    int removedNode = nodeSet[index];

    bool isEx=distribution1(generator1);
    if(!isEx){
      dynamic_bitset<> solu_prim= solution;
      vec CTemp = sol_C;
      CTemp.elem(pattern[removedNode]->edges) -=1;

      //find neighbor
      dynamic_bitset<> tdb(pattern.size());
      for(int i=0;i<pattern[removedNode]->edges.n_elem;i++){
        tdb |= madjacencyList[pattern[removedNode]->edges(i)];
      }
      dynamic_bitset<> neigh = tdb;
      neigh ^= solu_prim;
      tdb &= neigh;


      vector<int> possSwap;
      it= tdb.find_first();
      while(it!=dynamic_bitset<>::npos){
        uvec z = (CTemp.elem(pattern[it]->edges) < g->C.elem(pattern[it]->edges));
        bool states = all(z);
        if(states){
          possSwap.push_back(it);
        }
        it=tdb.find_next(it);
      }

      if(!possSwap.empty()){
        std::default_random_engine generator2(r());
        std::uniform_int_distribution<int> distribution2(0,possSwap.size()-1);
        int swapNode = possSwap[distribution2(generator2)];
        sol_C.elem(pattern[swapNode]->edges) += 1;
        sol_C.elem(pattern[removedNode]->edges) -= 1;
        solution[removedNode]=0;
        solution[swapNode]=1;
        nodes.push_back(swapNode);
//        cout<<"insert "<<swapNode<<"\t remove "<<removedNode<<endl;
      }
      else nodes.push_back(removedNode);

    }
    else{
      nodes.push_back(removedNode);
    }
    nodeSet[index] = nodeSet[nodeSet.size()-1];
    nodeSet.pop_back();
  }

  nodeSet = nodes;
}

int partialCount(vector<Motif*>& patterns){
  int size=patterns.size();
  if(size==0) {
		return 0;
	}

  //some motifs must exist
  vector<Motif*> pattern;
  int index = 0;
  int count = 0;
  for(int i=0;i<size;i++){
    Motif* a = patterns[i];
    uvec z = (g->used.elem(a->edges) <= g->C.elem(a->edges));
    bool state = all(z);
    if(state){
      count++;
      g->C.elem(a->edges) -= 1;
      delete a;
    }
    else{
      a->id = index;
      index++;
      pattern.push_back(a);
    }
  }
  patterns.clear();
  if(index == 0){
		return count;
	}
  size = index;

  //build overlap Graph
  vector<dynamic_bitset<> > madjacencyList(g->EdgeSize, dynamic_bitset<>(size));
  dynamic_bitset<>::size_type it;
  for(int i=0;i<size;i++){
    Motif* a=pattern[i];
    for(int j=0;j<a->edges.n_elem;j++){
      madjacencyList[a->edges(j)][i]=1;
    }
  }

  dynamic_bitset<> best(size); //best solution
  vec CTemp = zeros(g->EdgeSize);
  vector<int> nodeSet;
  generateSolution(pattern, best, CTemp, nodeSet);
  localSearch(pattern, madjacencyList, best, CTemp, nodeSet);
  int con = 1;

  while(con < 20){
    int cur = best.count();
    pertub(pattern, madjacencyList, best, CTemp, nodeSet);
    localSearch(pattern, madjacencyList, best, CTemp, nodeSet);
    if(best.count()==cur){
      con++;
    }
    else{
      con = 1;
    }
  }

  for(int i=0;i<pattern.size();i++) {
		delete pattern[i];
	}

  pattern.clear();
  return count+best.count();
}

int F2count(vector<Motif*>& pattern){
  if (pattern.size() < 2) {
    return pattern.size();
  }
  int countF2 = 0;
  int size = pattern.size();
  //build overlap graph
  vector<dynamic_bitset<> > madjacencyList(g->EdgeSize, dynamic_bitset<>(size));
  dynamic_bitset<>::size_type it;
  for(int i=0;i<size;i++){
    Motif* a=pattern[i];
    for(int j=0;j<a->edges.n_elem;j++){
      madjacencyList[a->edges(j)][i]=1;
    }
  }

  dynamic_bitset<> neighbors(size);
  heap::fibonacci_heap<Motif*,heap::compare<CompareNode> > pq;//priority queue to store the priority of motifs
  map<int, cn > tool;

  for(int i=0;i<size;i++){
     Motif* a=pattern[i];
     neighbors.reset();
     for(int j=0;j<a->edges.n_elem;j++){
        int index = a->edges(j);
        neighbors |= madjacencyList[index];
     }
     neighbors[i]=0;

     if(neighbors.none()){
        countF2++;
     }
     else{
        a->degree = neighbors.count();
        tool.insert(pair<int,cn>(i,pq.push(a)));
     }
  }

  dynamic_bitset<> updateEdge(g->EdgeSize);
  dynamic_bitset<> nn(size);//record which embeddings are to update
  dynamic_bitset<> tdb(size);
  while (!pq.empty()) {
    Motif* a = pq.top();
    countF2++;

    //delete a and its neighbors
    neighbors.reset();
    for(int j=0;j<a->edges.n_elem;j++){
       int index = a->edges(j);
       neighbors |= madjacencyList[index];
    }

    //get updateEdge
    updateEdge.reset();
    it = neighbors.find_first();
    while(it!=dynamic_bitset<>::npos){
       Motif* b = pattern[it];
       for(int i=0;i<b->edges.n_elem;i++)
          updateEdge[b->edges(i)]=1;
       it=neighbors.find_next(it);
    }
    for(int i=0;i<a->edges.n_elem;i++)
       updateEdge[a->edges(i)]=0;

    it=neighbors.find_first();
   	 while(it!=dynamic_bitset<>::npos){
          pq.erase(tool[it]); //remove motif from the priority queue
   	      it=neighbors.find_next(it);
   	 }

     nn.reset();
     it=updateEdge.find_first();
     while(it!=dynamic_bitset<>::npos){
        tdb = madjacencyList[it];
        tdb &= neighbors;
        madjacencyList[it] ^= tdb;
        nn |= madjacencyList[it];
        it = updateEdge.find_next(it);
    }

    it=nn.find_first();
    while(it!=dynamic_bitset<>::npos){
        Motif* b=pattern[it];

        neighbors.reset();
        for(int j=0;j<b->edges.n_elem;j++){
          int index = b->edges(j);
          neighbors |= madjacencyList[index];
        }
        neighbors[it]=0;

        if(neighbors.none()){
            countF2++;
            pq.erase(tool[it]);
        }

       else{
          b->degree = neighbors.count();
          pq.update(tool[it]);
       }

       it=nn.find_next(it);
    }
  }

  updateEdge.clear();
  nn.clear();
  tdb.clear();
  neighbors.clear();
	madjacencyList.clear();
	tool.clear();
	return countF2;
}

void count(int i){
  const clock_t begin_f1 = clock();
  vector<Motif*> pattern;
  g->findMotifs(pattern,i);
  const clock_t end_f1 = clock();
  float f1_time = 1000 * float( end_f1 - begin_f1 )/CLOCKS_PER_SEC;
  cout << pattern.size()<<"\t"<<f1_time<<"\t";


  const clock_t begin_f2 = clock();
  int countF2 = F2count(pattern);
  const clock_t end_f2 = clock();
  float f2_time = 1000 * float(end_f2 - begin_f2) / CLOCKS_PER_SEC + f1_time;
  cout<<countF2<<"\t"<<f2_time<<"\t";

  const clock_t begin_partial = clock();
  int count = partialCount(pattern);
  const clock_t end_partial = clock();
  float partial_time = 1000 * float(end_partial - begin_partial) / CLOCKS_PER_SEC + f1_time;
  cout<<count<<"\t"<<partial_time<<"\t";

}

int main(int argc, char *argv[]){
  int numberNodes = atoi(argv[1]);
  cout<<"F1\ttime\tF2\ttime\tPartial\ttime"<<endl;
  string typeName[] = {"Feedforward","Bifan","Biparallel","Cascade"};

  for(int j=1;j<=4;j++){
		cout<<"motif pattern: "<<typeName[j-1]<<endl;
    g = new Graph(numberNodes, argv[2]);
    count(j);
    delete g;
		cout<<endl;
  }


}

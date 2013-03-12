#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <map>
#include <vector>

template<class T>
class Graph{
public:
  Graph() {}
	void AddNode (T node) {
		nodes.push_back(node);
		neighbors[node] = vector<T>();
	}
	
	void AddEdge (T node1, T node2) {
		if (neighbors.count(node1)*neighbors.count(node2) == 0) {
			cerr << "AddEdge should be called after AddNode" << endl;
			return;
		}
		neighbors[node1].push_back(node2);
		neighbors[node2].push_back(node1);
	}
	
public:
	vector<T> nodes;
	map<T, vector<T> > neighbors;
};

#endif

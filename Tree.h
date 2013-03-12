#ifndef __TREE_H__
#define __TREE_H__

#include <map>
#include <vector>

template<class T>
class Tree{
public:
  Tree() {}

	Tree (T _root): root(_root) {
		children[_root] = vector<T>();
	}
	
	// This could only be called when root is not set yet.
	void SetRoot (T _root) {
		if (children.count(_root) == 0) {
			root = _root;
		}
		else {
			cerr << "SetRoot is only allowed if root has not been set before!" << endl << endl;
			return;
		}
	}
	
	void AddChild (T parent, T child) {
		if (children.count(child) > 0) {
			cerr << "Tree should not contain loops!" << endl << endl;
			return;
		}
		else if (children.count(parent) == 0){
			cerr << "AddChild could only add a child to a existing node in the tree." << endl << endl;
			return;
		}
		children[parent].push_back(child);
		children[child] = vector<T>();
	}
	
	void DepthFirstTraverse (T cur, vector<T> & res){
		res.push_back(cur);
		for (int i = 0 ; i < children[cur].size(); i++) {
			DepthFirstTraverse (children[cur][i], res);
		}
	}
	
	vector<T> DepthFirstTraverse (){
		vector<T> res;
		DepthFirstTraverse (root, res);
		return res;
	}
	
public:
	T root;
	map<T, vector<T> > children;
};

#endif

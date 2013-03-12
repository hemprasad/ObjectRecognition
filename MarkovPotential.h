#include <Coordinate.h>
#include <SMessage.h>
#include <Tree.h>
#include <math.h>
#include <algorithm>
#include <utility>

// This class takes care of the markov potentials (phi_ij's) for you! 
// It learns parameters from the training data, and then can apply the
//  learned models to compute the unary potentials on new images.
//

class MarkovPotential
{
public:
  MarkovPotential(const vector<GTImage> &training_data, int _step = 3, int _template_size=15) : template_size(_template_size), step(_step){
	    // learn a little 31x31 template for each part. The template consists of a
	    //  estimate of mean gradient magnitue and variance in gradient magnitude 
	    //  at each position.
	    //
	    
	    Coordinate sqr_means[6][6]; // A helper that stores mean of (Li-Lj)^2
	    
		// train
		for(int img=0; img<training_data.size(); img++) {
			const GTImage &gt = training_data[img];	
			// now update the model for this image, for each of 6 parts
			for(int i=0; i<6; i++) {
				for(int j = i+1; j < 6; j ++) {
					Coordinate diff, delta, sqr_delta;
					// diff is the difference between Li and Lj in this image
					diff = gt.coords[i] - gt.coords[j];
					delta = diff - means[i][j];
					sqr_delta = diff.sqr() - sqr_means[i][j];
					means[i][j] = means[i][j] + delta/(img+1);
					sqr_means[i][j] = sqr_means[i][j] + sqr_delta/(img+1);
				}
			}
		}
		
		for(int i=0; i<6; i++) {
			for(int j=i+1; j<6; j++) {
				vars[i][j] = sqr_means[i][j] - means[i][j].sqr();
			}
			for(int j=0; j<i; j++) {
				means[i][j] = - means[j][i];
				vars[i][j] = vars[j][i];
			}
		}
		cerr << "Means of rows" << endl;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++)
				cerr << means[i][j].row << "\t";
			cerr << endl;
		}
		
		cerr << "Means of cols" << endl;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++)
				cerr << means[i][j].col << "\t";
			cerr << endl;
		}

		cerr << "Vars of rows" << endl;				
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++)
				cerr << vars[i][j].row << "\t";
			cerr << endl;
		}
		
		cerr << "Vars of cols" << endl;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++)
				cerr << vars[i][j].col << "\t";
			cerr << endl;
		}		
	}
		
	void compute_markov_logpotential (vector<SDoublePlane> &unary_potentials, vector<SMessage> & message, Tree<int> model, int node){
		message[node] = SMessage(unary_potentials[node].rows(), unary_potentials[node].cols(), model.children[node].size());
		SMessage &result = message[node];
		double psi_ij;
		if (model.children[node].empty()) {
		// Leaf Nodes
			for (int i = 0; i < result.rows(); i++) {
    			for (int j = 0; j < result.cols(); j++) {
    				result.vals[i][j] = unary_potentials[node][i][j];
    			}
    		}
    		cerr << "Finished computing markov logpotential for part "<< node<< endl;
    	}
		else {
			// Wait until it has received message from all its children
			for (int c_index = 0; c_index < model.children[node].size(); c_index++) {
				int child = model.children[node][c_index];
				compute_markov_logpotential(unary_potentials, message, model, child);
			}
			
			// Outer Loop
			for (int ix = 0; ix < result.rows(); ix+=step) {
				for (int iy = 0; iy < result.cols(); iy+=step) {
					double max_val = -1.0e100;
					Coordinate max_index;
					result.vals[ix][iy] = unary_potentials[node][ix][iy];
    				// Inner Loop
    				for (int c_index = 0; c_index < model.children[node].size(); c_index++) {
    					int child = model.children[node][c_index];
						// find the max and min of message from child.
						// then calculate the upper and lower bound of Li that might lead to a maximum value of Psi_i(Lj) + Psi_ij(Li,Lj)
						for (int jx = 0; jx < result.rows(); jx += step) {
							for (int jy = 0; jy < result.cols(); jy += step) {
    							psi_ij = - (ix - jx - means[node][child].row) * (ix - jx - means[node][child].row) / (2 * vars[node][child].row)
    									 - (iy - jy - means[node][child].col) * (iy - jy - means[node][child].col) / (2 * vars[node][child].col);
    							if (psi_ij + message[child].vals[jx][jy] > max_val) {
    								max_val = psi_ij + message[child].vals[jx][jy];
    								max_index = Coordinate(jx, jy);
    							}
    						}
    					}
	    					result.vals[ix][iy] += max_val;
	    					result.coords[c_index][ix][iy] = max_index;
    				}
    				// End of Inner Loop
//		    		cerr << "Finished computing argmax for ("<< ix <<", "<< iy << ")" << endl;
				}
//				cerr << "Finished computing argmax for row " << ix << endl;
			}
			// End of Outer Loop
    		cerr << "Finished computing markov logpotential for part "<< node<< endl;
		}
	}
		
	void compute_markov_logpotential (vector<SDoublePlane> & unary_potentials, vector<SMessage> & message, Tree<int> model){
		return compute_markov_logpotential (unary_potentials, message, model, model.root);
	}
	
	// computes log potentials given a new image, will be using unary_potentials as a parameter.
	void compute_markov_argmax(vector<SDoublePlane> &unary_potentials, Tree<int> model, vector<Coordinate> &result) 
	{
		cerr << "Computing Markov potentials..." << endl;
		
		// message[i] = message sent from i to its parent
		vector<SMessage> message(6);

		compute_markov_logpotential(unary_potentials, message, model);
		
		double max_val = -1.0e100;
		Coordinate max_index;
		for (int ix = 0; ix < message[model.root].rows(); ix+=step) {
			for (int iy = 0; iy < message[model.root].cols(); iy+=step) {
				if (message[model.root].vals[ix][iy] > max_val) {
					max_val = message[model.root].vals[ix][iy];
					max_index = Coordinate(ix, iy);
				}
			}
		}
		
//		cerr << "Depth First Traverse" << endl;
		vector<int> order = model.DepthFirstTraverse();
		result[model.root] = max_index;
//		cerr << "Starting from " << model.root << endl;
		for (int p_index = 0; p_index < order.size(); p_index++) {
			int parent = order[p_index];
			Coordinate p_max_index = result[parent];
//			cerr << parent << " expanding" << endl;
			for (int c_index = 0; c_index < model.children[parent].size(); c_index++) {
				int child = model.children[parent][c_index];
				result[child] = message[parent].coords[c_index][int(p_max_index.row)][int(p_max_index.col)];
			}
		}
		cerr << "Finished computing Markov potentials..." << endl;
	}

	void compute_loopy_message (vector<SDoublePlane> &unary_potentials, vector<vector<SDoublePlane> > & message, Graph<int> model, int node){
		vector<int> & neighbors = model.neighbors[node];
		for (int i = 0; i < neighbors.size(); i++) {
			int nb_to = neighbors[i];
			SDoublePlane &marginal = message[node][node];
			SDoublePlane &result = message[node][nb_to];
			double psi_ij;
			
			// Outer Loop
			for (int ix = 0; ix < result.rows(); ix+=step) {
				for (int iy = 0; iy < result.cols(); iy+=step) {
					double max_val = -1.0e100;
					Coordinate max_index;
					result[ix][iy] = unary_potentials[node][ix][iy];
					marginal[ix][iy] = unary_potentials[node][ix][iy];
    				// Inner Loop
    				for (int n_index = 0; n_index < neighbors.size(); n_index++) {

    					int nb_from = neighbors[n_index];

						for (int jx = 0; jx < result.rows(); jx += step) {
							for (int jy = 0; jy < result.cols(); jy += step) {
								psi_ij = - (ix - jx - means[node][nb_from].row) * (ix - jx - means[node][nb_from].row) / (2 * vars[node][nb_from].row)
										 - (iy - jy - means[node][nb_from].col) * (iy - jy - means[node][nb_from].col) / (2 * vars[node][nb_from].col);
								if (psi_ij + message[nb_from][node][jx][jy] > max_val) {
									max_val = psi_ij + message[nb_from][node][jx][jy];
								}
							}
						}
						if (nb_from != nb_to) {
		    				result[ix][iy] += max_val;
		    			}
	    				marginal[ix][iy] += max_val;
    				}
    				// End of Inner Loop
				}
			}
			// End of Outer Loop
    		cerr << "Finished computing Loopy BP message from part "<< node << " to part " << nb_to << endl;
		}
	}
	
	void compute_loopy_argmax(vector<SDoublePlane> &unary_potentials, Graph<int> model, vector<Coordinate> &result, int iteration) 
	{
		cerr << "Computing Loopy BP messages..." << endl;
		
		// message[i] = message sent from i to its parent
		vector<vector<SDoublePlane> > message(6);
		
		// initialize
		for (int node = 0; node < 6; node++) {
			message[node] = vector<SDoublePlane> (6);
			// Use message[node][node] to store marginal of node
			message[node][node] = SDoublePlane(unary_potentials[node].rows(),unary_potentials[node].cols());
			for (int ix = 0; ix < unary_potentials[node].rows(); ix++) {
				for (int iy = 0; iy < unary_potentials[node].cols(); iy++) {
					message[node][node][ix][iy] = unary_potentials[node][ix][iy];
				}
			}
			for (int n_index = 0; n_index < model.neighbors[node].size(); n_index++) {
				int nb_to = model.neighbors[node][n_index];
				message[node][nb_to] = SDoublePlane(unary_potentials[node].rows(),unary_potentials[node].cols());
				for (int ix = 0; ix < unary_potentials[node].rows(); ix++) {
					for (int iy = 0; iy < unary_potentials[node].cols(); iy++) {
						message[node][nb_to][ix][iy] = unary_potentials[node][ix][iy];
					}
				}
			}
		}
		// iteratively compute messages
		for (int i = 0; i < iteration; i++) {
			for (int node = 0; node < 6; node++) {
				compute_loopy_message (unary_potentials, message, model, node);
				cerr << "Iteration " << i << ": Computing Loopy BP messages for part " << node << endl;
			}
		}
		
		
		for(int ii=0; ii<6; ii++) {
			double max_val=-1e100;
			Coordinate max_index;
			SDoublePlane &marginal = message[ii][ii];
			for(int i=0; i<marginal.rows(); i++) {
				for(int j=0; j<marginal.cols(); j++) {
					if(marginal[i][j] > max_val) {
						max_val = marginal[i][j];
						max_index = Coordinate(i, j);
					}
				}
				result[ii] = max_index;
			}
		}
		cerr << "Finished computing Loopy BP argmax..." << endl;
	}
	public:
		Coordinate means[6][6], vars[6][6];
		int template_size; 	// default value is 15
		int step; 			// default value is 3
};

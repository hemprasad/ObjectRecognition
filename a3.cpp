#include <SImage.h>
#include <SImageIO.h>
#include <Coordinate.h>
#include <Tree.h>
#include <Graph.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

using namespace std;

// Ground truth information for training or testing image
class GTImage
{
public:
  string image_filename;
  Coordinate coords[6];
};

// Code for computing unary potential functions
//  Feel free to modify this if you want, but 
//  it should work as-is.
#include <UnaryPotential.h>
#include <MarkovPotential.h>

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlay rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
    input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}



// Load ground truth coordinates and image names from file
//
vector<GTImage> load_groundtruth(const string &gt_filename)
{
  // read in training images, one at a time, and put the ground truth coordinates
  //  into a vector.
  ifstream ifs(gt_filename.c_str());
  
  vector<GTImage> training_images;
  while(ifs.good())
    {
      string image_number;
      ifs >> image_number;
      if(!ifs.good())	break;
      
      GTImage gt;
      gt.image_filename = string("imgs/image_") + image_number + ".png";
      for(int ii=0; ii<6; ii++)
	ifs >> gt.coords[ii].col >> gt.coords[ii].row;

      training_images.push_back(gt);
    }

  return training_images;
}


// Draw a given part configuration on an image, using simple
//  colored rectangles.
// Color code is: 
// left eye: bright green
// right eye: dark green
// nose: bright red
// left mouth: bright blue
// right mouth: dark blue
// chin: dark red
void draw_configuration(const char *output_fname, const char *sample_filename, const vector<Coordinate> &configuration)
{
  int g_code[] = {255, 128,   0,   0,   0, 0};
  int r_code[] = {0, 0,     255,   0,   0, 128};
  int b_code[] = {0, 0,       0, 255, 128, 0};

  SDoublePlane img = SImageIO::read_png_file(sample_filename);
  SDoublePlane r = img, g = img, b = img;
  for(int ii=0; ii<configuration.size(); ii++)
    {
      overlay_rectangle(g, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, g_code[ii], 3);
      overlay_rectangle(r, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, r_code[ii], 3);
      overlay_rectangle(b, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, b_code[ii], 3);
    }

  SImageIO::write_png_file(output_fname, r, g, b);
}


int main(int argc, char *argv[])
{
  try {
    if(!(argc == 4)) {
      cerr << "usage:   " << argv[0] << " training_gt_file testing_gt_file model_id" << endl;
      cerr << "example: " << argv[0] << " gt.train gt.test c" << endl;
      return 1;
    }

    string train_gt_filename(argv[1]), test_gt_filename(argv[2]);
    
    const char model_id = argv[3][0];
    if (model_id < 'a' || model_id > 'e') {
    	cerr << "model_id should be a, b, c, d or e" << endl;
    	return 1;
    }
    
    // This is a parameter used to compute the unary potentials. It specifies
    //  how big the part appearance models should be.
    const int template_size=15;

    // Load in the training data
    vector<GTImage> training_data = load_groundtruth(train_gt_filename);
    vector<GTImage> testing_data = load_groundtruth(test_gt_filename);

    // Now set up the unary potential functions
    UnaryPotential unary_potentials(training_data, template_size);
    cerr << "Done learning unary potentials" << endl << endl;

    // Now set up the markov potential functions
    MarkovPotential markov_potentials(training_data, 5);
    cerr << "Done learning markov potentials" << endl << endl;

	// Define FacialModels 
	map<char,Tree<int> > Models;

	// Model b
	Models['b'] = Tree<int> (0);
	Models['b'].AddChild(0, 1);
	Models['b'].AddChild(0, 2);
	Models['b'].AddChild(0, 3);
	Models['b'].AddChild(0, 4);
	Models['b'].AddChild(0, 5);

	// Model c
	Models['c'] = Tree<int> (2);
	Models['c'].AddChild(2, 1);
	Models['c'].AddChild(2, 4);
	Models['c'].AddChild(2, 5);
	Models['c'].AddChild(1, 0);
	Models['c'].AddChild(4, 3);
	
	// Model d
	Models['d'] = Tree<int> (2);
	Models['d'].AddChild(2, 0);
	Models['d'].AddChild(2, 1);
	Models['d'].AddChild(2, 3);
	Models['d'].AddChild(2, 4);
	Models['d'].AddChild(2, 5);

	map<char,Graph<int> > GModels;
	GModels['e'] = Graph<int>();
	GModels['e'].AddNode(0);
	GModels['e'].AddNode(1);
	GModels['e'].AddNode(2);
	GModels['e'].AddNode(3);
	GModels['e'].AddNode(4);
	GModels['e'].AddNode(5);
	GModels['e'].AddEdge(0,1);
	GModels['e'].AddEdge(0,3);
	GModels['e'].AddEdge(1,2);
	GModels['e'].AddEdge(1,4);
	GModels['e'].AddEdge(2,5);
	GModels['e'].AddEdge(3,5);
	GModels['e'].AddEdge(4,5);
		
	
	//     
    // Now run naive inference on all the sample images
    // 
	ofstream outfile;
	outfile.open ("result.txt", ios::out | ios::app);
	
	outfile << "--------------------------------"
			<< " Model " << model_id << " "
			<< "--------------------------------"
			<< endl;
	int part_count = 6;
	int total_correct_count = 0;
    for(int img=0; img<testing_data.size(); img++) {
    	const GTImage &gt = testing_data[img];
        vector<Coordinate> best_configuration(part_count);

		// Compute using UnaryPotentials
        vector<SDoublePlane> my_potentials = unary_potentials.compute_unary_logpotential(gt.image_filename);
        
		if (model_id == 'a') {
			for(int ii=0; ii<part_count; ii++) {
				double max_val=-1e100;
				for(int i=0; i<my_potentials[ii].rows(); i++) {
					for(int j=0; j<my_potentials[ii].cols(); j++) {
						if(my_potentials[ii][i][j] > max_val) {
							best_configuration[ii].row = i;
							best_configuration[ii].col = j;
							max_val = my_potentials[ii][i][j];
						}
					}
				}
			}
		}
    	// Compute using MarkovPotentials
		else if (model_id == 'b') {
    		markov_potentials.compute_markov_argmax(my_potentials, Models['b'], best_configuration);
		}
		else if (model_id == 'c') {
    		markov_potentials.compute_markov_argmax(my_potentials, Models['c'], best_configuration);
		}    	    	
		else if (model_id == 'd') {
    		markov_potentials.compute_markov_argmax(my_potentials, Models['d'], best_configuration);
		}
		else if (model_id == 'e') {
			markov_potentials.compute_loopy_argmax(my_potentials, GModels['e'], best_configuration, 3);
		}

    	// Evaluate
		int image_correct_count = 0;    	
	    for (int i = 0; i < part_count; i ++){
    		Coordinate diff = best_configuration[i] - gt.coords[i];
	    	if (fabs(diff.row) <= 20 && fabs(diff.col) <= 20){
	    		image_correct_count ++;
	    	}
    	}
    	total_correct_count += image_correct_count;
    	
    	// Store result
    	string image_number = gt.image_filename.substr(gt.image_filename.find('_') + 1, 4);
    	outfile << image_number << " ";
	    for (int i = 0; i < part_count; i ++){
    		outfile << int(best_configuration[i].col) << " " << int(best_configuration[i].row) << " ";
    	}
    	outfile << image_correct_count << endl;

    	// Visualize result by drawing detected parts on the image, and 
    	//  outputing to a new image file.
    	//
    	string result_imagename = string("result/result_") + model_id + "_" + image_number + ".png";
    	draw_configuration(result_imagename.c_str(), gt.image_filename.c_str(), best_configuration);
    }
    outfile << "Accuracy: " << double(total_correct_count) / (testing_data.size() * part_count) << endl;
    outfile.close();

  } catch(std::string &str) {
    cerr << str << endl;
  }

  return 0;
}

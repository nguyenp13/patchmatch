
/*

This is our main executable that wraps our PatchMatch implementation. See patchmatch.h for code relating to the PatchMatch implementation. 

TODO: 
    check the types to make sure they all make sense
    Check to make sure everything works
    add timing mechanism

*/

#include <chrono>
#include <cmath>
#include "patchmatch.h"
#include "array.h"
#include "pfm.h"

using std::cout;
using std::endl;
using std::string;

void usage() {
    fprintf(stderr, "\n"
                    "usage: main <input_image_a>.png <input_image_b>.png <output_file>.pfm <options>\n"
                    "\n"
                    "This is an implementation of PatchMatch. Please see http://gfx.cs.princeton.edu/pubs/Barnes_2009_PAR/index.php for details about the PatchMatch algorithm and about how to tune and use these parameters. This implementation uses the SSD metric for patch comparison. \n"
                    "\n"
                    "\n"
                    "Required Parameters: \n"
                    "\n"
                    "    <input_image_a>.png: This is the file name of input image A. Must be a .png file. \n"
                    "\n"
                    "    <input_image_b>.png: This is the file name of input image B. Must be a .png file. \n"
                    "\n"
                    "    <output_file>.pfm: This string is the name of the desired output file. This will be a .pfm file containing the nearest neighbor field. The output file will be 3 dimensional. It will have the same X and Y dimensions as input image A (both minus the size of patch_dim). It's Z dimension will have length 3 (to describe the X and Y coordinates of input image b that represents the nearest neighbor patch as well as the patch distance). Patches of size patch_dim by patch_dim will be referred to by their upper left coordinate. In order to get the coordinates (x_b,y_b) in input image b that correspond to the coordinates (x_a,y_a) in input image A using this output file, we will have to use a pfm reader to extract the values at output_file[y_a,x_a,0] to get x_b and output_file[y_a,x_a,1] to get y_b. \n"
                    "\n"
                    "\n"
                    "There are also several other parameters with default values that can change if specified, e.g. \"./main input_a.png input_b.png output.pfm -patch_dim 5 -random_search_size_exponent 0\". \n"
                    "\n"
                    "\n"
                    "Options: \n"
                    "\n"
                    "    -patch_dim <patch_dim>: This int value determines the size of the patches. The patches will be size <patch_dim> by <patch_dim>. \n"
                    "\n"
                    "    -num_iterations <num_iterations>: This int value is the fixed number of times the PatchMatch algorithm will run before terminating. "
                    "\n"
                    "    -random_search_size_exponent <random_search_size_exponent>: This int value determines the largest neighborhood size around a patch in the nearest neighbor field to conduct the random search. When we are conducting the random search around patch p in our nearest neighbor field, the largest neighborhood around p we will search will be of size 2^<random_search_size_exponent>. \n"
                    "\n"
                    "    -random_search_attempts <random_search_attempts>: This int value defines how many neighbors we will compare to in each iteration of the random search. \n"
                    "\n"
           );
    exit(1);
}

int main(int argc, char* argv[]) {
    
    if (argc<4) {
        usage();
    }

    std::srand( std::time(NULL) );

    static const char* A_name = argv[1];
    static const char* B_name = argv[2];
    static const char* output_name = argv[3];
    static const int patch_dim = atoi(get_command_line_param_val_default_val(argc, argv, "-patch_dim", "5"));
    static const int num_iterations = atoi(get_command_line_param_val_default_val(argc, argv, "-num_iterations", "4"));
    static const int random_search_size_exponent = atoi(get_command_line_param_val_default_val(argc, argv, "-random_search_size_exponent", "3"));
    static const int num_random_search_attempts = atoi(get_command_line_param_val_default_val(argc, argv, "-random_search_attempts", "8"));
    
    static const png::image< png::rgba_pixel > A_image(A_name); 
    static const png::image< png::rgba_pixel > B_image(B_name); 
    static const int A_height = A_image.get_height();
    static const int A_width = A_image.get_width();
    static const int B_height = B_image.get_height();
    static const int B_width = B_image.get_width();
    static const int Ann_height = A_height-patch_dim+1;
    static const int Ann_width = A_width-patch_dim+1;
    static const Array<byte> A(A_image);
    static const Array<byte> B(B_image);
    Array<int> Ann(vector<int>{Ann_height, Ann_width, 3});
    
    NEWLINE;
    PRINT("Parameter Values");
    TEST(A_height);
    TEST(A_width);
    TEST(B_height);
    TEST(B_width);
    TEST(Ann_height);
    TEST(Ann_width);
    TEST(patch_dim);
    TEST(num_iterations);
    TEST(random_search_size_exponent);
    TEST(num_random_search_attempts);
    
    long total_patch_distance;
    double mean_patch_distance;
    
    long initial_total_patch_distance = 0;
    #pragma omp parallel for reduction(+:initial_total_patch_distance)
    for (int yy = 0; yy < Ann_height; ++yy) {
        for (int xx = 0; xx < Ann_width; ++xx) {
            initial_total_patch_distance += Ann(yy,xx,D_COORD);
        }
    }
    NEWLINE;
    cout << "Initial Total Patch Distance: " << initial_total_patch_distance << endl;
    cout << "Initial Mean Patch Distance:  " << DOUBLE(initial_total_patch_distance)/DOUBLE(Ann_height*Ann_width) << endl;
    fflush(stdout);
    
    PRINT(1);
    patchmatch(A, B, Ann, A_height, A_width, B_height, B_width, Ann_height, Ann_width, patch_dim, num_iterations, random_search_size_exponent, num_random_search_attempts, total_patch_distance, mean_patch_distance);
    
    // Write output to .pfm file
    
    float *depth = new float[Ann_height*Ann_width*3];
    for (int y = 0; y < Ann_height; y++) {
        for (int x = 0; x < Ann_width; x++) {
            int i = (Ann_height-1-y)*Ann_width*3+x*3;
            
            depth[i] = FLOAT(Ann(y,x,X_COORD));
            depth[i+1] = FLOAT(Ann(y,x,Y_COORD));
            depth[i+2] = FLOAT(Ann(y,x,D_COORD));
        }
    }
    write_pfm_file3(output_name, depth, Ann_width, Ann_height);
    
    return 0;
}


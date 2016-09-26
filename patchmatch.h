
/*

This is an implementation of PatchMatch. Please see http://gfx.cs.princeton.edu/pubs/Barnes_2009_PAR/index.php for details about the PatchMatch algorithm. 

*/

#pragma once

#ifndef PATCHMATCH_H
#define PATCHMATCH_H

#include "util.h"
#include "array.h"

#define Y_COORD 0
#define X_COORD 1
#define D_COORD 2

int patch_SSD(
            const Array<byte> &A, 
            const Array<byte> &B, 
            const int &ax, 
            const int &ay, 
            const int &bx, 
            const int &by, 
            const int &patch_dim
        ) {
    double score = 0;
    int channels = A.channels();
    for(int dy=0;dy<patch_dim;dy++) {
        for(int dx=0;dx<patch_dim;dx++) {
            for(int z=0;z<channels;z++) {
                score += SQUARE(A(ay+dy,ax+dx,z)-B(by+dy,bx+dx,z));
            }
        }
    }
    return score;
}

void patchmatch(
            const Array<byte> &A, 
            const Array<byte> &B, 
            Array<int> &Ann, 
            const int &A_height, 
            const int &A_width, 
            const int &B_height, 
            const int &B_width, 
            const int &Ann_height, 
            const int &Ann_width, 
            const int &patch_dim, 
            const int &num_iterations, 
            const int &random_search_size_exponent, 
            const int &num_random_search_attempts, 
            long &total_patch_distance, 
            double &mean_patch_distance
        ) {
    
    // Randomize the nearest neighbor field 
    PRINT("RANDOMIZE");
    for(int y=0;y<Ann_height;y++) { 
        for(int x=0;x<Ann_width;x++) { 
                Ann(y,x,Y_COORD)=RAND_INT(0, B_height-patch_dim+1); 
                Ann(y,x,X_COORD)=RAND_INT(0, B_width-patch_dim+1); 
                Ann(y,x,D_COORD)=patch_SSD(A, B, x, y, Ann(y,x,X_COORD), Ann(y,x,Y_COORD), patch_dim); 
        } 
        if (y%100==0) {
            TEST(y);
        }
    } 
    PRINT("RANDOMIZE END");
    
    bool going_down_and_right = true; 
    
    PRINT("ITERATION START");
    for(int iteration_index=0; iteration_index<num_iterations; iteration_index++) { 
        TEST(iteration_index);
        int delta, start_x, start_y; 
        if (going_down_and_right) { 
            delta = 1; 
            start_x = 0; 
            start_y = 0; 
        } else { 
            delta = -1; 
            start_x = Ann_height-1; 
            start_y = Ann_height-1; 
        }
        
        // Belief Propogation 
        PRINT("Belief Propogation START");
        for(int y=start_y; 0<=y && y<Ann_height; y+=delta) { 
            for(int x=start_x; 0<=x && x<Ann_width; x+=delta) { 
                int by;
                int bx;
                // Vertical offset
                by = Ann(y-delta,x,Y_COORD)+delta; 
                bx = Ann(y-delta,x,X_COORD); 
                if (0<=by && by<B_height-patch_dim+1) {
                    int ay = y; 
                    int ax = x; 
                    int old_patch_distance = Ann(y,x,D_COORD);
                    int new_patch_distance = patch_SSD(A, B, ax, ay, bx, by, patch_dim);
                    if (new_patch_distance<old_patch_distance) {
                        Ann(y,x,Y_COORD) = by;
                        Ann(y,x,X_COORD) = bx;
                        Ann(y,x,D_COORD) = new_patch_distance;
                    }
                } 
                // Horizontal offset
                by = Ann(y,x,Y_COORD); 
                bx = Ann(y,x-delta,X_COORD)+delta; 
                if (0<=bx && bx<B_width-patch_dim+1) {
                    int ay = y; 
                    int ax = x; 
                    int old_patch_distance = Ann(y,x,D_COORD);
                    int new_patch_distance = patch_SSD(A, B, ax, ay, bx, by, patch_dim);
                    if (new_patch_distance<old_patch_distance) {
                        Ann(y,x,Y_COORD) = by;
                        Ann(y,x,X_COORD) = bx;
                        Ann(y,x,D_COORD) = new_patch_distance;
                    }
                } 
            } 
        } 
        ///////////////////////////////////////
        total_patch_distance = 0;
        // Calculate total and mean patch distances
        for(int y=0; y<Ann_height; y++) { 
            for(int x=0; x<Ann_width; x++) { 
                total_patch_distance += LONG(Ann(y,x,D_COORD));
            }
        }
        mean_patch_distance = DOUBLE(total_patch_distance)/DOUBLE(Ann_height*Ann_width);
        NEWLINE;
        NEWLINE;
        PRINT("Pre-random search");
        TEST(total_patch_distance);
        TEST(mean_patch_distance);
        NEWLINE;
        ///////////////////////////////////////
        PRINT("Belief Propogation end");
        going_down_and_right = !going_down_and_right; 
        
        // Random Search
        for(int y=start_y;y<Ann_height;y+=delta) { 
            for(int x=start_x;x<Ann_width;x+=delta) { 
                int bx = Ann(y,x,X_COORD);
                int by = Ann(y,x,Y_COORD);
                for(int radius_index=random_search_size_exponent; 0<radius_index; radius_index--) {
                    int radius = INT(pow(2,radius_index));
                    int search_box_min_x = MAX(0,bx-radius);
                    int search_box_min_y = MAX(0,by-radius);
                    int search_box_max_x = MIN(B_width-patch_dim+1,bx+radius);
                    int search_box_max_y = MIN(B_height-patch_dim+1,by+radius);
                    for(int search_attempt_count=0; search_attempt_count<num_random_search_attempts; search_attempt_count++) {
                        int bx_new = RAND_INT(search_box_min_x,search_box_max_x+1);
                        int by_new = RAND_INT(search_box_min_y,search_box_max_y+1);
                        int old_patch_distance = Ann(y,x,D_COORD);
                        int new_patch_distance = patch_SSD(A, B, x, y, bx_new, by_new, patch_dim);
                        if (new_patch_distance<old_patch_distance) {
                            Ann(y,x,Y_COORD) = by_new;
                            Ann(y,x,X_COORD) = bx_new;
                            Ann(y,x,D_COORD) = new_patch_distance;
                        }
                        
                    }
                }
            }
        }
        ///////////////////////////////////////
        total_patch_distance = 0;
        // Calculate total and mean patch distances
        for(int y=0; y<Ann_height; y++) { 
            for(int x=0; x<Ann_width; x++) { 
                total_patch_distance += LONG(Ann(y,x,D_COORD));
            }
        }
        mean_patch_distance = DOUBLE(total_patch_distance)/DOUBLE(Ann_height*Ann_width);
        PRINT("Post-random search");
        TEST(total_patch_distance);
        TEST(mean_patch_distance);
        NEWLINE;
        NEWLINE;
        ///////////////////////////////////////
    }
    PRINT("ITERATION end");
    
    total_patch_distance = 0;
    // Calculate total and mean patch distances
    for(int y=0; y<Ann_height; y++) { 
        for(int x=0; x<Ann_width; x++) { 
            total_patch_distance += LONG(Ann(y,x,D_COORD));
        }
    }
    mean_patch_distance = DOUBLE(total_patch_distance)/DOUBLE(Ann_height*Ann_width);
}

#endif // PATCHMATCH_H


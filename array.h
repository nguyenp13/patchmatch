
/*

This is a general purpose multidimensional array class. 

It was written with image processing purposes in mind. 

The class has a png++ dependency.

*/

#pragma once

#ifndef ARRAY_H
#define ARRAY_H

#include <png++/png.hpp>
#include <cmath>
#include "util.h"

#define CLAMP_FLOAT(p) ( MIN(1.0,MAX(0.0,p)) )
#define CLAMP_INT(p) ( MIN(255,MAX(0,p)) )
#define CLAMP(val,min,max) ( MIN(max,MAX(min,val)) )
#define PIXEL(r,g,b) ( png::rgba_pixel(INT(CLAMP_INT(r)),INT(CLAMP_INT(g)),INT(CLAMP_INT(b)),255) )

#define R_COORD 0
#define G_COORD 1
#define B_COORD 2

using std::string;
using std::to_string;
using std::ostream;
using std::vector;
using std::cout;
using std::endl;
using std::round;

template<class real>
class Array {
    public:
        real *data;
        vector<int> sizes;
        vector<int> stride;
        int nelems;
        
        void resize(const vector<int> &sizes_) {
            if (sizes == sizes_) { return; }
            
            delete[] data;
            sizes = sizes_;
            
            stride.resize(sizes.size());
            nelems = 1;
            for (int i = stride.size()-1; i >= 0; i--) {
                stride[i] = nelems;
                nelems *= sizes[i];
            }
            data = new real[nelems];
        }
        
        void assign(const Array &other) {
            resize(other.sizes);
            #pragma omp parallel for
            for(unsigned i=0; i < nelems; i++){
                data[i] = other.data[i];
            }
        }
        
        Array() :data(NULL) {
            resize(vector<int>{1});
        }
        
        Array(const Array &other) :data(NULL) {
            assign(other);
        }
        
        Array(const vector<int> &sizes_) :data(NULL) {
            resize(sizes_);
        }
        
        Array(const png::image< png::rgba_pixel > &png_image) :data(NULL) {
            int h = png_image.get_height();
            int w = png_image.get_width();
            resize(vector<int>{h, w, 3});
            for (int y=0; y<h; y++) {
                for (int x=0; x<w; x++) {
                    (*this)(y,x,R_COORD) = (real)(png_image[y][x].red);
                    (*this)(y,x,G_COORD) = (real)(png_image[y][x].green);
                    (*this)(y,x,B_COORD) = (real)(png_image[y][x].blue);
                }
            }
        }
        
        ~Array() {
            delete[] data;
        }
        
        void save_to_png(char* output_name) {
            const int output_name_length = strlen(output_name);
            ASSERT(sizes.size() == 2 || (sizes.size() == 3 && channels()==3), "Need array either to be 3D with 3 channels in the third dimension or to be 2D to save to png");
            ASSERT(!strcmp(&output_name[output_name_length-4],".png"), "output_name must have the .png extension");
            png::image< png::rgba_pixel > output_image(this->width(), this->height());
            if ( IS_INT((real)0) || IS_BYTE((real)0) ) {
                for (int y=0; y<height(); y++) {
                    for (int x=0; x<width(); x++) {
                        if (sizes.size() == 2) {
                            output_image[y][x] = PIXEL((*this)(y,x), (*this)(y,x), (*this)(y,x));
                        } else {
                            output_image[y][x] = PIXEL((*this)(y,x,R_COORD), (*this)(y,x,G_COORD), (*this)(y,x,B_COORD));
                        }
                    }
                }
            } else if (IS_FLOAT((real)0) || IS_DOUBLE((real)0)) {
                for (int y=0; y<height(); y++) {
                    for (int x=0; x<width(); x++) {
                        if (sizes.size() == 2) {
                            output_image[y][x] = PIXEL(INT(round((*this)(y,x)*255)), INT(round((*this)(y,x)*255)), INT(round((*this)(y,x)*255)));
                        } else {
                            output_image[y][x] = PIXEL(INT(round((*this)(y,x,R_COORD)*255)), INT(round((*this)(y,x,G_COORD)*255)), INT(round((*this)(y,x,B_COORD)*255)));
                        }
                    }
                }
            } else {
                ASSERT(false, "save_to_png() is not supported for types other than int, unsigned char, float, or double");
            } 
            output_image.write(output_name);
        }
        
        void clear(const real &val=0) {
            #pragma omp parallel for
            for(unsigned i=0; i<nelems; i++){
                data[i]=val;
            }
        }
        
        real sum() {
            real sum = 0;
            for(unsigned i=0; i<nelems; i++){
                sum += data[i];
            }
            return sum;
        }
        
        real product() {
            real product = 0;
            for(unsigned i=0; i<nelems; i++){
                product *= data[i];
            }
            return product;
        }
        
        Array& normalize() {
            real sum = this->sum();
            #pragma omp parallel for
            for(unsigned i=0; i<nelems; i++){
                data[i] /= sum;
            }
            return (*this);
        }
        
        int height() const {
            ASSERT(sizes.size() >= 2, "Need at least 2 dimensions to get height");
            return sizes[0];
        }

        int width() const {
            ASSERT(sizes.size() >= 2, "Need at least 2 dimensions to get width");
            return sizes[1];
        }

        int channels() const {
            return (sizes.size()<3) ? 1 : sizes[2];
        }
        
        int dimensions() const {
            return sizes.size();
        }
        
        Array<real>& rgb2gray() {
        	ASSERT(sizes.size() >= 3, "Need at least 3 dimensions to convert to grayscale");
        	
        	real* data_grayscale = new real[height() * width()];
	        for (int y = 0; y < height() ; y++){
		        for (int x = 0; x < width() ; x++){
		            real cum_sum = 0;
		            for (int z = 0; z < channels(); z++){
                        double channel_weight = 1.0/channels();
                        if (channels() == 3) {
                            if (z == 0) { channel_weight = 0.299; }
                            else if (z == 1) { channel_weight = 0.5870; }
                            else if (z == 2) { channel_weight = 0.1140; }
                        }
				        cum_sum += (*this)(y,x,z)*channel_weight;
			        }
			        data_grayscale[y*width()+x] = cum_sum;
		        }
	        }
	        
            delete[] data;
            stride.resize(2);
            stride[0]=width();
            stride[1]=1;
            sizes.resize(2); 
            data = data_grayscale;
            return *this;
        }
        
        string str() const {
            string ans("[");
            if (dimensions() == 1) {
                for (int i = 0; i < nelems; i++) {
                    ans += to_string(data[i]);
                    if (i < nelems - 1) { ans += ", "; }
                }
            } else if (dimensions() == 2) {
                for (int y = 0; y < height(); y++) {
                    ans += "[";
                    for (int x = 0; x < width(); x++) {
                        ans += to_string((*this)(y,x));
                        if (x < width()-1) { ans += ", "; }
                    }
                    ans += "]";
                    if (y < height()-1) { ans += ",\n "; }
                }
            }else if (dimensions() == 3) {
                for (int z = 0; z < channels(); z++) {
                    ans += "\nChannel "+to_string(z)+":\n ";
                    for (int y = 0; y < height(); y++) {
                        ans += "[";
                        for (int x = 0; x < width(); x++) {
                            ans += to_string((*this)(y,x,z));
                            if (x < width()-1) { ans += ", "; }
                        }
                        ans += "]";
                        if (y < height()-1) { ans += ",\n "; }
                    }
                    ans += "\n";
                }
            } else {
                ASSERT(false, "Array::str() not supported for dimension > 3");
            }
            ans += "]";
            return ans;
        }
        
        real& operator()(int v0) const {
            ASSERT(sizes.size() == 1, (string("1D lookup in array of dimensionality ")+to_string(sizes.size())).c_str()); 
            ASSERT(v0 >= 0 && v0 < sizes[0], (string("1D lookup out of bounds (index=")+to_string(v0)+", length="+to_string(sizes[0])+")").c_str()); 
            return data[v0];
        }
        
        real& operator()(int v0, int v1) const {
            ASSERT(sizes.size() == 2, (string("2D lookup in array of dimensionality ")+to_string(sizes.size())).c_str()); 
            ASSERT(v0 >= 0 && v0 < sizes[0] && v1 >= 0 && v1 < sizes[1], (string("2D lookup out of bounds (row=")+to_string(v0)+", column="+to_string(v1)+", height="+to_string(sizes[0])+", width="+to_string(sizes[1])+")").c_str()); 
            return data[v0*stride[0]+v1];
        }
        
        real& operator()(int v0, int v1, int v2) const {
            ASSERT(sizes.size() == 3, (string("3D lookup in array of dimensionality ")+to_string(sizes.size())).c_str());
            ASSERT(v0 >= 0 && v0 < sizes[0] && v1 >= 0 && v1 < sizes[1] && v2 >= 0 && v2 < sizes[2], (string("2D lookup out of bounds (row=")+to_string(v0)+", column="+to_string(v1)+", channel="+to_string(v1)+", height="+to_string(sizes[0])+", width="+to_string(sizes[1])+", num_channels="+to_string(sizes[2])+")").c_str()); 
            return data[v0*stride[0]+v1*stride[1]+v2];
        }
};

template<class real>
ostream& operator<<(ostream& os, const Array<real>& obj) {
    os << obj.str();
    return os;
}

template<class real, class pixel_type_in, class pixel_type_out>
void convolution_filter(const Array<real> &kernel, const Array<pixel_type_in> &I, Array<pixel_type_out> &J) {
    ASSERT(I.sizes.size()>=2, "I must be at least 2 dimensional");
    ASSERT(kernel.sizes.size()==2, "Convolutions with kernels of dimension other than 2 have not been implemented yet");
    ASSERT(kernel.height()==kernel.width(), "Kernel must be square, i.e. the height and width must be the same");
    ASSERT(MOD(kernel.height(),2)==1, "Kernel height and width must be odd");
    
    J.resize(I.sizes);
    
    int const &I_height = I.height();
    int const &I_width = I.width();
    int const &I_channels = (I.sizes.size()==2) ? 1 : I.channels();
    int const kernel_half_dim = kernel.height()/2;
    
    #pragma omp parallel for
    for(int y=0; y<I_height; y++) {
        for(int x=0; x<I_width; x++) {
            pixel_type_out final_value[I_channels];
            for(int channel=0; channel<I_channels; channel++) {
                final_value[channel] = 0;
            }
            for (int row=-kernel_half_dim; row<kernel_half_dim+1; row++) {
                for (int column=-kernel_half_dim; column<kernel_half_dim+1; column++) {
                    double spatial_weight = kernel(row+kernel_half_dim, column  +kernel_half_dim);
                    if (I.sizes.size()==2) {
                        final_value[0] += spatial_weight*I(CLAMP(y+row,0,I_height-1),CLAMP(x+column,0,I_width-1));
                    } else {
                        for(int channel=0; channel<I_channels; channel++) {
                            final_value[channel] += (pixel_type_out)(spatial_weight*(pixel_type_out)I(CLAMP(y+row,0,I_height-1),CLAMP(x+column,0,I_width-1),channel));
                        }
                    }
                }
            }
            
            if (I.sizes.size()==2) {
                J(y,x) = final_value[0];
            } else {
                for(int channel=0; channel<I_channels; channel++) {
                    J(y,x,channel) = final_value[channel];
                }
            }
        }
    }
}

void get_gaussian_kernel(const int &dim, const double &sigma, Array<double> &kernel) {
    ASSERT(MOD(dim,2)==1, "Gaussian kernel must be of odd dimension");
    kernel.resize(vector<int>{dim,dim});
    
    #pragma omp parallel for
    for(int y=0; y<dim; y++) {
        for(int x=0; x<y+1; x++) {
            kernel(y,x) = G(y-dim/2, x-dim/2, sigma);
            if (x != y) {
                kernel(x,y) = kernel(y,x);
            }
        }
    }
    kernel.normalize();
}

#endif // ARRAY_H


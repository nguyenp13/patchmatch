
/* 

This header is used for the importation of several miscellaneous macros, functions, typedefs, etc. that I commonly find helpful. 

*/

#pragma once

#ifndef UTIL_H
#define UTIL_H

#include <cassert>
#include <iostream>
#include <cstring>
#include <typeinfo>
#include <vector>
#include <fstream>

using std::string;
using std::to_string;
using std::vector;
using std::cout;
using std::endl;

#define DEBUG BUILD_DEBUG
#define ASSERT(x, msg) if (DEBUG && !(x)) { fprintf(stderr, "Assertion failed in %s(%d): %s\n", __FILE__, __LINE__, msg); assert(false); exit(1); }

#define INT(x) ((int)(x))
#define DOUBLE(x) ((double)(x))
#define FLOAT(x) ((float)(x))
#define BYTE(x) ((unsigned char)(x))
#define SHORT(x) ((short)(x))
#define LONG(x) ((long)(x))
#define LONG_LONG(x) ((long long)(x))
#define U_LONG_LONG(x) ((unsigned long long)(x))

#define PI (3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709)

#define SQUARE(x) ((x)*(x))
#define MOD(a,b) (((INT(a)%INT(b))+INT(b))%INT(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))
#define IS_ODD(x) (MOD(x,2))
#define IS_EVEN(x) (!(MOD(x,2)))

#define RAND_FLOAT (DOUBLE(rand())/DOUBLE(RAND_MAX)) // Returns a random float in the range [0,1]
#define RAND_INT(lower, higher) (MOD(rand(),higher-lower)+lower) // Returns a random int in the range [lower,higher)

#define QUIT exit(1);
#define TEST(x) cout << (#x) << ": " << x << endl; fflush(stdout);
#define PRINT(x) cout << x << endl; fflush(stdout);
#define DIVIDER cout << "======================================================================================================================================" << endl; fflush(stdout);
#define NEWLINE cout << endl; fflush(stdout);

#define IS_INT(a) same_type(a,INT(0))
#define IS_DOUBLE(a) same_type(a,DOUBLE(0))
#define IS_FLOAT(a) same_type(a,FLOAT(0))
#define IS_BYTE(a) same_type(a,BYTE(0))
#define IS_LONG(a) same_type(a,LONG(0))
#define IS_LONG_LONG(a) same_type(a,LONG_LONG(0))
#define IS_U_LONG_LONG(a) same_type(a,U_LONG_LONG(0))

#define CHAR_STAR_EQUAL(a,b) !strcmp(a,b)

typedef unsigned char byte;

template <class type1, class type2>
bool same_type(const type1 &a, const type2 &b) {
    return typeid(a) == typeid(b);
}

double G(const double &x, const double &sigma){
    return exp(-(x*x)/(2*SQUARE(sigma))) / sqrt(2*PI*SQUARE(sigma));
}

double G(const double &x, const double &y, const double &sigma){
    return exp(-(SQUARE(y)+SQUARE(x))/(2*SQUARE(sigma))) / (2*PI*SQUARE(sigma));
}

template <class real>
real lerp(const real &a, const real &b, const double &t) { // Linear Interpolation
    ASSERT(t>=0 && t<=1, (string("Interpolation factor (")+to_string(t)+") must be in interval [0.0,1.0]").c_str());
    return a*t+(1-t)*b;
}

/*  String Manipulation and File Reading */

const char* get_command_line_param_val_default_val(const int &argc, char** &argv, const char* param_option, const char* default_val) {
    for (int i=0;i<argc; i++) {
        if (CHAR_STAR_EQUAL(argv[i],param_option)) {
            if (i+1 != argc) {
                return argv[i+1];
            } else {
                return default_val;
            }
        }
    }
    return default_val;
}

void split(const string &line, const string &delimiter, vector<string> &vector_of_strings, bool skip_empty=false) {
    auto start = 0U;
    auto end = line.find(delimiter);
    while (end != std::string::npos){
        if (!skip_empty || start != end) {
            vector_of_strings.push_back( line.substr(start, end - start) );
        } 
        start = end + delimiter.length();
        end = line.find(delimiter, start);
    }
    if (delimiter != string("\n")) { // since files have an extra newline at the end, we need to make sure we don't count the extra empty line we may get if our delimiter is the newline
        vector_of_strings.push_back( line.substr(start, end) );
    }
}

void write_file(const char* const &filename, string &input) { 
    std::ofstream out(filename);
    out << input;
    out.close();
}

void write_file(const string &filename, string &input) { 
    write_file(filename.c_str(), input);
}

void read_file(const char* const &filename, string &output) { 
    FILE* f = fopen(filename, "r");
    char* buffer;
    if (!f) {
        fprintf(stderr, "Unable to open %s for reading\n", filename);
        return;
    }
    
    fseek(f, 0, SEEK_END);
    int length = ftell(f);
    fseek(f, 0, SEEK_SET);
    
    buffer = (char*)malloc(length+1);
    length = fread((void*)buffer, 1, length, f);
    fclose(f);
    buffer[length] = '\0';
    
    output = string(buffer);
    
    free(buffer);
}

void read_file(const string &filename, string &output) { 
    read_file(filename.c_str(), output);
}

#endif // UTIL_H


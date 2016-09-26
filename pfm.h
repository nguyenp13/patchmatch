
/* 

The code below for operating on .pfm files is modified from that of Connelly Barnes 

*/

#pragma once

#ifndef PFM_H
#define PFM_H

#include "util.h"

int is_little_endian() {
    if (sizeof(float) != 4) { 
        PRINT("Bad float size."); 
        QUIT;
    }
    byte b[4] = { 255, 0, 0, 0 };
    return *((float*)b) < 1.0;
}

void write_pfm_file3(const char *filename, float *depth, int w, int h) {
    FILE *f = fopen(filename, "wb");
    static const int channels = 3;
    double scale = is_little_endian() ? -1.0 : 1.0;
    fprintf(f, "PF\n%d %d\n%lf\n", w, h, scale);
    for (int i = 0; i < w*h*channels; i++) {
        float d = depth[i];
        fwrite((void *) &d, 1, 4, f);
    }
    fclose(f);
}

#endif // PFM_H


# patchmatch

Implementation of PatchMatch. 

Finding correspondences between patches (small neighborhoods of pixels in an image) is useful in many graphics and vision applications such as optical flow, hole filling, or super-resolution. Naive exhaustive searches can be slow as they can be O(n**2) where n is the number of pixels in any of our input images. PatchMatch is a linear time patch correspondence algorithm that approximates the optimal patch correspondences between two images. For further details and source code, see http://gfx.cs.princeton.edu/gfx/pubs/Barnes_2009_PAR/index.php


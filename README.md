This is the repository for Samuel Bergé's bachelor's thesis in 2025.


I thank Prof.Dr.Jan-Henrik Haunert and Dr.Jonas Sauer for supervising this work and providing helpful comments. I also 
want to thank Alexander Naumann, guiding me through the first steps of using CGAL and helped with getting used to work
with geospatial data.


The project deals with polygon aggregation. We want to aggregate polygons (f.e. buildings in a street map) such that
the orientation of the edges is preserved as much as possible. In this thesis, we established different methods:

- version 0: Baseline-EE, edge extension until hiting another polygon or box
- version 1: Baseline-EE with additional thresholds to edge extension lengths
            - th_variant 0: Constant threshold
            - th_variant 1: Linear threshold
            - th_variant 2: Square root threshold
- version 3: EdgeRelink (merging extended edges via EdgeRelink before identifying the new faces)
- version 3.1: EdgeRelink with threshold option (th_variant 0,1,2)
- version 5: ConnectOuter (merging extended edges via ConnectOuter before identifying the new faces)
- version 5.1: ConnectOuter with threshold option (th_variant 0,1,2)

Further versions (2 and 4) were implemented but not evaluated in the thesis:
- version 2: pre-process by taking shortest/middle/longest edge of relaxed similar groups
- version 4: using subdivision of input data and running first on each subset, then on the results


USAGE

--alpha DOUBLE
Alpha in objective function (Standard: 0.01)

--version STRING 
As introduced above (Standard: "0")

--threshold DOUBLE
Length of extension for constant threshold (Standard: 50)

--th_scale DOUBLE
Factor for linear and rooted threshold extension (Standard: 10)

--degree DOUBLE
Degree (tau) limitation for merging via EdgeRelink or ConnectOuter (Standard: 5)

--distance DOUBLE
Distance (delta) limitation factor, that is multiplied by the avg. edge length of the input data. 
Also for EdgeRelink and ConnectOuter (Standard: 0.75)


--thresholdvariant INT
Threshold variant (constant, linear, square root) for edge extension (Standard: 0)

--power DOUBLE
For subdivision grid size. Size of subsets on which instance should be run first, should be pow(n, power) with n being
number of input edges. (Standard: 0.5)


Exectution example:
./aggregator /path/to/input.gpkg /path/to/output/folder --alpha=0.01 --version=3.1 --th_variant=2 --th_scale=27

This would run version 3.1 (EdgeRelink with threshold) with square root threshold variant and a scaling factor of 27 on
alpha 0.01. 


While this repository contains the code that was written for this project, for comparison with the delaunay triangulation
in the bachelorsthesis, we used code provided internal by Rottmann et al. (as referenced in the thesis).





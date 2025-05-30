This is the repository for Samuel Bergé's bachelor's thesis in 2025.

The project deals with polygon aggregation. We want to aggregate polygons (f.e. buildings in a street map) such that
the orientation of the edges is preserved as much as possible. 

You will find extern code pieces used in this project, of which most parts come from colleagues of the University of 
Bonn. They already built a framework in which polygons get aggregated by using delauney triangulation. For that they
established a graph representation of the polygons such that a min cut algorithm will find an solution that is optimal 
for their objective function. 

In this project we use that min cut algorithm as a black box and try to implement and optimize a variation to the delauney
triangulation in which, as said, the edge orientation of the polygons should be preserved. 
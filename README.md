# 3-approximation-for-Minimum-Convex-Partition-Problem
A C++ implementation for 3-approximation for Minimum Convex Partition Problem for INF562 course at École Polytechnique

### Described in the article:
Christian Knauer and Andreas Spillner: Approximation Algorithms for the Minimum Convex Partition Problem, SWAT 2006, pp. 232-241

It uses this library (https://github.com/nlohmann/json) to parse json files.

It uses this wrapper (https://github.com/lava/matplotlib-cpp) of python's matplotlib library to print the phases of the recursion.

### To compile
```g++ -o main main.cpp -I/usr/include/python2.7 -lpython2.7```

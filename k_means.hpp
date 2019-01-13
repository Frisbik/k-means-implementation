#ifndef K_MEANS_H
#define K_MEANS_H

#include<tuple>
#include "center_initialization.hpp"

/* number of data points / observations */
int n;
/* dimension of space */
int d;

/* infinite value */
typedef std::vector<std::vector<int>> Clusters;
typedef std::tuple<std::vector<double>, double, int> Results;

/* initializes data points by reading the data file */
void init_datapoints(std::vector<DataPoint*>& dtp, std::string& datafile);

/* void solves the instance with the specified algorithm */
Results solve_k_means(const std::vector<DataPoint*>& points, const int& num_center, const int& num_iter,
   const char& method, const int& num_candidates);

/* creates the assignment */
Clusters cluster_clients(const std::vector<int>& a);

/* computes the cost of the given assignement */
double compute_cost(const std::vector<double>& dist);

/* applies Lloyds algorithm with datapoints, centers and assignment */
void Lloyds(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assignment, Clusters& clusters);

/* computes the centroid of a given assignement to a center */
DataPoint compute_centroid(const std::vector<DataPoint*>& points, const std::vector<int>& assigned);

#endif
#ifndef DATAPOINT_H
#define DATAPOINT_H

#include<limits>
#include "datapoint.hpp"

const double infty = std::numeric_limits<double>::infinity();

/* chooses center at random */
void init_centers_random(std::vector<DataPoint*>&centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers);

/* chooses centers using D2 seeding */
void init_centers_D2(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers);

/* chooses centers using GREEDY fashion way */
void init_centers_T_greedy(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers, const int num_cand);

/* applies GREEDY algorithm to choose centers */
void init_centers_greedy(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers);

/* computes cost of solution with additional center */
double candidate_contribution(const std::vector<DataPoint*>& points, const std::vector<double>& dist,
    const DataPoint& candidate);

/* stores new distances and assignement, given a new center */
void update_distance(const std::vector<DataPoint*>& points, std::vector<double>& dist, std::vector<int>& assign,
    const DataPoint& c, const int& idx);

#endif
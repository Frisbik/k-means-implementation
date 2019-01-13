#include<random>
#include<assert.h>
#include "center_initialization.hpp"

extern int n;
extern int d;

/* chooses center at random */
void init_centers_random(std::vector<DataPoint*>&centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers)
{
    /* chooses k centers uniformly at random, we cannot take the same center multiple times */
     /* random library */
    std::random_device rd;
    std::mt19937 gen(rd());

    /* Start by choosing a first center at random from the set of points */
    std::uniform_int_distribution<> dis(0,n-1);
    std::vector<int> center_idx;

    for( unsigned int i = 0; i < num_centers; i++ )
    {
        bool isIn = true;
        while( isIn )
        {
            int idx = dis(gen);
            isIn = false;

            for( unsigned int j = 0; j < center_idx.size(); j++ )
                if( center_idx[j] == idx )
                    isIn = true;

            if( !isIn )
                center_idx.emplace_back(idx);
        }

        /* add center to vector, updates distances and assignment */
        centers.emplace_back(new DataPoint(points[center_idx.back()]->get_coordinates()));
        update_distance(points, dist, assign, *centers.back(), i);
    }

    assert(centers.size() == num_centers);
}

/* chooses centers using D2 seeding */
void init_centers_D2(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers)
{
    /* random library */
    std::random_device rd;
    std::mt19937 gen(rd());

    /* Start by choosing a first center at random from the set of points */
    std::uniform_int_distribution<> dis(0,n-1);

    /* Add the chosen point to the set of centers */
    DataPoint* dp = new DataPoint(points[dis(gen)]->get_coordinates());
    centers.emplace_back(dp);
    //centers.back()->print();

    /* Update the distance of points to the chosen center */
    for( unsigned int i = 0; i < n; i++ )
    {
        /* update distance function: minimun{current distance, distance to added center} */
        dist[i] = distance(*points[i], *centers.back());
        /* since there is only one center, there is also only one cluster. Each point is assigned to it */
        assign[i] = 0;
    }

    /* Then for 1,...,k-1, choose a center with D2-seeding */
    for( unsigned int i = 1; i < num_centers; i++ )
    {
        /* select a center at random depending on weights and add it to centers set */
        std::discrete_distribution<int> rand_center(dist.begin(), dist.end());
        DataPoint* dp2 = new DataPoint(points[rand_center(gen)]->get_coordinates());
        centers.emplace_back(dp2);
        //centers.back()->print();

        /* update distance = min{old_distance, distance_to_added_center}, cluster index and weights */
        for( unsigned int j = 0; j < n; j++ )
        {
            if( dist[j] > distance(*points[j], *centers.back()) )
            {
                dist[j] = distance(*points[j], *centers.back());
                assign[j] = i;
            }
        }
    }
}

/* chooses centers using GREEDY fashion way */
void init_centers_T_greedy(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers, const int num_cand)
{
    /* random library */
    std::random_device rd;
    std::mt19937 gen(rd());

    /* Start by choosing a first center at random from the set of points */
    std::uniform_int_distribution<> dis(0,n-1);

    /* Add the chosen point to the set of centers */
    centers.emplace_back(new DataPoint(points[dis(gen)]->get_coordinates()));
    //centers.back()->print();

    /* Update the distance of points to the chosen center */
    for( unsigned int i = 0; i < points.size(); i++ )
    {
        /* since it is the first center added, all distances are infinite, update them */
        dist[i] = distance(*points[i], *centers.back());
        /* since there is only one center, there is also only one cluster. Each point is assigned to it */
        assign[i] = 0;
    }

    /* greedily choose the center out of #num_cand candidates, which contributes the most to the decrease */
    for( unsigned int i = 1; i < num_centers; i++ )
    {
        /* vectors to store: index of candidates, objective cost of each candidate */
        std::vector<int> cand_idx;
        std::vector<double> cand_costs;

        /* sets the random distribution given distances */
        std::discrete_distribution<int> rand_center(dist.begin(), dist.end());

        /* select candidates */
        for( unsigned int T = 0; T < num_cand; T++ )
        {
            /* pick a candidate at random and ensure that it was not picked before */
            bool isIn = true;
            //while(isIn)
            {
                int idx = rand_center(gen);
                isIn = false;

                for( unsigned int l = 0; l < cand_idx.size(); l++ )
                    if( cand_idx[l] == idx )
                       isIn = true;

                if( !isIn )
                {
                    cand_idx.emplace_back(idx);
                    cand_costs.emplace_back(candidate_contribution(points, dist, *points[idx]));
                }
            }

            /* computes the objective cost by adding the selected candidate */
            //cand_costs.emplace_back(candidate_contribution(points, dist, *points[cand_idx[T]]));
        }
        assert(cand_idx.size() <= num_cand);
        assert(cand_costs.size() <= num_cand);

        /* Find minimum of all costs */
        double min = cand_costs[0];
        int chosed = cand_idx[0];
        for( unsigned int j = 0; j < cand_costs.size(); j++ )
            if( min > cand_costs[j] )
            {
                chosed = cand_idx[j];
                min = cand_costs[j];
            }

        /* add center */
        centers.emplace_back(new DataPoint(points[chosed]->get_coordinates()));

        /* since a new center has been added, the distance vector changes as well as the clustering vector */
        update_distance(points, dist, assign, *centers.back(), i);
    }

    assert(centers.size() == num_centers);
}

/* applies GREEDY algorithm to choose centers */
void init_centers_greedy(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assign, const int num_centers)
{
    /* add k centers to center list */
    for( unsigned int j = 0; j < num_centers ; j++ )
    {
        int best_center_idx;
        double best_center_cost = infty;
        for( unsigned int i = 0; i < points.size(); i++ )
        {
            double d = candidate_contribution(points, dist, *points[i]);
            if( best_center_cost > d )
            {
                best_center_idx = i;
                best_center_cost = d;
            }
        }

        /* adds best candidate to set of centers */
        centers.emplace_back(new DataPoint(points[best_center_idx]->get_coordinates()));

        /* update distance since number we added a center */
        update_distance(points, dist, assign, *centers.back(), j);
    }
    assert(centers.size() == num_centers);
}

/* computes cost of solution with additional center */
double candidate_contribution(const std::vector<DataPoint*>& points, const std::vector<double>& dist,
    const DataPoint& candidate)
{
    double value = 0.0;
    for( unsigned int i = 0; i < points.size(); i++ )
        value += std::min(dist[i], distance(*points[i], candidate));

    return value;
}

/* stores new distances and assignement, given a new center */
void update_distance(const std::vector<DataPoint*>& points, std::vector<double>& dist, std::vector<int>& assign,
    const DataPoint& c, const int& idx)
{
    /* the distance is equal to min{previous_distance, distance_to_new_center} and assignment changes only if
     * the distance changes.
     */

    // assert points.size() = dist.size() = n;
    for( unsigned int i = 0; i < points.size(); i++ )
    {
        double d = distance(*points[i], c);
        if( dist[i] > d )
        {
            dist[i] = d;
            assign[i] = idx;
        }
    }
}
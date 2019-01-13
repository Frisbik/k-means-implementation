#include<iostream>
#include<vector>
#include<sstream>
#include<fstream>
#include<string>
#include<random>            // uniform_int_distribution
#include<limits>            // infinity
#include<algorithm>         //std::min
#include<chrono>            //time elapsed
#include<tuple>


#include "k_means.hpp"

int main(int argc, char const *argv[])
{

    /* in this function, we need to call k-means solver */
    /* as an input, we need to give multiples values of k, which seeding algorithm should be used */
    /* when using greedy, we should also give the number of candidates that should be used */
    /* Should be returned: time, cost, #iterations, cost at each iterations */
    /* An input of the following function is:
     * INPUT:   (datafile / #number_executions / k / max_iterations / use Random_seeding(Y/N) /
     *          use D2-seeding(Y/N) / use T_Greedy-seeding(Y/N)) / use GREEDY (Y/N).
     *
     * OUTPUT:  file containing desired results, located in "data/datafile_k.dat.
     */

    if( argc < 8 )
    {
        std::cout << "Invalid number of arguments in " << argv[0] << " file." << std::endl;
        return 1;
    }

    std::string file = argv[1];
    int num_executions = atoi(argv[2]);
    int k = atoi(argv[3]);
    int max_iter = atoi(argv[4]);
    bool use_random = (*argv[5] == 'Y') ? true : false;
    bool use_D2 = (*argv[6] == 'Y') ? true : false;
    bool use_T_greedy = (*argv[7] == 'Y') ? true : false;
    bool use_greedy = (*argv[8] == 'Y') ? true : false;
    std::vector<int> num_candidates;

    /* at least one seeding method has to used. */
    if( !use_random && !use_D2 && !use_T_greedy && !use_greedy )
        return 0;

    if( use_T_greedy )
    {
        num_candidates.emplace_back(2);
        num_candidates.emplace_back(int(ceil(log(k))));
        num_candidates.emplace_back(int(k));
        num_candidates.emplace_back(int(ceil(k*log(k))));
        num_candidates.emplace_back(int(k*k));
    }

    std::vector<DataPoint*> datapoints;

    /* initialize datapoints */
    init_datapoints(datapoints, file);
    std::cout << "dtp_size: " << datapoints.size() << std::endl;

    std::vector<Results> random_result;
    std::vector<Results> D2_result;
    std::vector<std::vector<Results>> T_greedy_result(num_candidates.size());
    Results greedy_result;

    for( unsigned int i = 0; i < num_executions; i++ )
    {
        std::cout << "Iteration: " << i << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        if( use_random )
            random_result.emplace_back(solve_k_means(datapoints, k, max_iter, 'N', 0));
        if( use_D2 )
            D2_result.emplace_back(solve_k_means(datapoints, k, max_iter, 'P', 0));
        if( use_T_greedy )
            for( unsigned int j = 0; j < num_candidates.size(); j++ )
                T_greedy_result[j].emplace_back(solve_k_means(datapoints, k, max_iter, 'T', num_candidates[j]));
    }
    /* greedy is not random */
    if( use_greedy )
            greedy_result = solve_k_means(datapoints, k, max_iter, 'G', 0);

    for( unsigned int i = 0; i < datapoints.size(); i++ )
        delete datapoints[i];

    /* create a file containining the cost at each iteration */
    int max = std::get<0>(random_result[0]).size();
    if( max < std::get<0>(D2_result[0]).size() )
        max = std::get<0>(D2_result[0]).size();

    for( unsigned int j = 0; j < num_candidates.size(); j++ )
        if( max < std::get<0>((T_greedy_result[j])[0]).size() )
            max = std::get<0>((T_greedy_result[j])[0]).size();

    if( use_greedy )
        if( max < std::get<0>(greedy_result).size() )
            max = std::get<0>(greedy_result).size();

    std::ofstream cost_evolution;
    std::string name;
    name = "experiments/" + file.substr(5, file.length()-14 ) + "_" + std::to_string(k) + "_means.dat";
    cost_evolution.open(name);


    for( unsigned int i = 0; i < max; i++ )
    {
        std::string str;
        std::ostringstream streamObj;
        str = std::to_string(i) + "                                                                                                                   ";
        if( i >= std::get<0>(random_result[0]).size() )
            streamObj << std::get<0>(random_result[0]).back();
        else
            streamObj <<  std::get<0>(random_result[0])[i];

        str.insert(4, streamObj.str());
        streamObj.str("");
        streamObj.clear();

        if( i >= std::get<0>(D2_result[0]).size() )
            streamObj << std::get<0>(D2_result[0]).back();
        else
            streamObj << std::get<0>(D2_result[0])[i];

        str.insert(18, streamObj.str());
        streamObj.str("");
        streamObj.clear();


        for( unsigned int j = 0; j < num_candidates.size(); j++ )
        {
            if( i >= std::get<0>((T_greedy_result[j])[0]).size() )
                streamObj << std::get<0>((T_greedy_result[j])[0]).back();
            else
                streamObj << std::get<0>((T_greedy_result[j])[0])[i];

            str.insert( 33 + 15*j, streamObj.str());
            streamObj.str("");
            streamObj.clear();
        }

        if( use_greedy )
        {
            if( i >= std::get<0>(greedy_result).size() )
                streamObj << std::get<0>(greedy_result).back();
            else
                streamObj << std::get<0>(greedy_result)[i];

            str.insert(33 + 15*(num_candidates.size()), streamObj.str());
            streamObj.str("");
            streamObj.clear();
        }

        cost_evolution  << str << "\n";
    }
    cost_evolution.close();

    /* prints average cost of each solution */
    double average_random = 0.0;
    double time_random = 0.0;
    double av_random_first = 0.0;
    double min_random = infty;
    double average_D2 = 0.0;
    double time_D2 = 0.0;
    double av_D2_first = 0.0;
    double min_D2 = infty;
    std::vector<double> average_T_greedy(num_candidates.size(), 0.0);
    std::vector<double> time_T_greedy(num_candidates.size(), 0.0);
    std::vector<double> av_T_greedy_first(num_candidates.size(), 0.0);
    std::vector<double> min_T_greedy(num_candidates.size(), infty);

    for( unsigned int i = 0; i < num_executions; i++ )
    {
        average_random += std::get<1>(random_result[i]);
        time_random += std::get<2>(random_result[i]);
        av_random_first += std::get<0>(random_result[i])[0];
        average_D2 += std::get<1>(D2_result[i]);
        time_D2 += std::get<2>(D2_result[i]);
        av_D2_first += std::get<0>(D2_result[i])[0];
        for( unsigned int j = 0; j < num_candidates.size(); j++ )
        {
            average_T_greedy[j] += std::get<1>((T_greedy_result[j])[i]);
            time_T_greedy[j] += std::get<2>((T_greedy_result[j])[i]);
            av_T_greedy_first[j] += std::get<0>((T_greedy_result[j])[i])[0];
        }

        if( min_random > std::get<1>(random_result[i]) )
            min_random = std::get<1>(random_result[i]);
        if( min_D2 > std::get<1>(D2_result[i]) )
            min_D2 = std::get<1>(D2_result[i]);
        for( unsigned int j = 0; j < num_candidates.size(); j++ )
            if( min_T_greedy[j] > std::get<1>((T_greedy_result[j])[i]) )
                min_T_greedy[j] = std::get<1>((T_greedy_result[j])[i]);
    }

    average_random = average_random / num_executions;
    time_random = time_random / num_executions;
    av_random_first /= num_executions;
    average_D2 = average_D2 / num_executions;
    time_D2 = time_D2 / num_executions;
    av_D2_first /= num_executions;
    for( unsigned int i = 0; i < average_T_greedy.size(); i++)
    {
        average_T_greedy[i] = average_T_greedy[i] / num_executions;
        time_T_greedy[i] = time_T_greedy[i] / num_executions;
        av_T_greedy_first[i] /= num_executions;
    }

    std::ofstream average_result;
    std::string name_result;
    name_result = "experiments/" + file.substr(5, file.length()-14 ) + "_" + std::to_string(k) + "_means_results.dat";
    average_result.open(name_result);

    average_result << "# average initial cost after " << num_executions << " trials. \n";
    average_result << av_random_first << "    " << av_D2_first << "   ";
    for( unsigned int i = 0; i < av_T_greedy_first.size(); i++ )
        average_result << av_T_greedy_first[i] << "   ";
    if( use_greedy )
        average_result << std::get<0>(greedy_result)[0];

    average_result << "\n" << "\n";
    average_result << "# average final cost after " << num_executions << " trials. \n";
    average_result << average_random << "    " << average_D2 << "   ";
    for( unsigned int i = 0; i < average_T_greedy.size(); i++ )
        average_result << average_T_greedy[i] << "   ";
    if( use_greedy )
        average_result << std::get<1>(greedy_result);

    average_result << "\n" << "\n";
    average_result << "# best objective value found after " << num_executions << " trials \n";
    average_result << min_random << "   " << min_D2 << "   ";
    for( unsigned int i = 0; i < time_T_greedy.size(); i++ )
        average_result << min_T_greedy[i] << "   ";
    if( use_greedy )
        average_result << std::get<1>(greedy_result);

    average_result << "\n" << "\n";
    average_result << "# average time after " << num_executions << " trials \n";
    average_result << time_random << "   " << time_D2 << "   ";
    for( unsigned int i = 0; i < time_T_greedy.size(); i++ )
        average_result << time_T_greedy[i] << "   ";
    if( use_greedy )
        average_result << std::get<2>(greedy_result);

    average_result.close();

    return 0;
}


/* void solves the instance with the specified algorithm */
Results solve_k_means(const std::vector<DataPoint*>& points,
    const int& num_center, const int& num_iter, const char& method, const int& num_candidates)
{
    std::vector<DataPoint*> centers;

    /* vector containing distance to closest center */
    std::vector<double> distances(points.size(), infty);
    /* vector containing the cluster a datapoint belongs to */
    std::vector<int> assignment(points.size());
    /* vector containing cost at each iteration */
    std::vector<double> potential;

    /* start measuring time */
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if( method == 'N' )
    {
        std::cout << "Random seeding for center selection" << std::endl;
        init_centers_random(centers, points, distances, assignment, num_center);
    }
    else if( method == 'P' )
    {
        /* chooses centers with D2 weighthing */
        std::cout << "D2-seeding for center selection" << std::endl;
        init_centers_D2(centers, points, distances, assignment, num_center);
    }
    else if( method == 'T' )
    {
        //if( num_candidates >= points.size() )
        //    return std::make_tuple(potential, 0.0, 0);

        /* chooses centers in a Greedy way with T candidates */
        std::cout << "Greedy selection with " << num_candidates << " candidates for centers" << std::endl;
        init_centers_T_greedy(centers, points, distances, assignment, num_center, num_candidates);
    }
    else if( method == 'G' )
    {
        /* chooses centers with greedy method */
        std::cout << "GREEDY-seeding for center selection" << std::endl;
        init_centers_greedy(centers, points, distances, assignment, num_center);
    }

    /* constructs clustering */
    Clusters clusters = cluster_clients(assignment);

    /* computes cost of the current assignement */
    potential.emplace_back(compute_cost(distances));

    /* apply Lloyds algorithm */
    for( unsigned int i = 1; i < num_iter; i++ )
    {
        Lloyds(centers, points, distances, assignment, clusters);

        // for( unsigned int i = 0; i < centers.size(); i++ )
        //     centers[i]->print();

        potential.emplace_back(compute_cost(distances));
        // assert new cost lower or equal
        clusters = cluster_clients(assignment);

        /* stop applying Lloyds algorithm if the cost does not change */
        if( potential[i-1] == potential[i] )
            break;

    }

    /* ends time */
    end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    for( unsigned int i = 0; i < centers.size(); i++ )
        delete centers[i];

    // for( unsigned int i= 0; i < potential.size(); i++ )
    //     std::cout << "Cost[" << i << "] is equal to: " << potential[i] << std::endl;

    std::cout << "Elapsed Time: " << elapsed_seconds << "ms\n" << std::endl;

    // return std::make_tuple(potential, potential.back(), elapsed_seconds);

    return Results(potential, potential.back(), elapsed_seconds);
}

/* initializes all data points */
void init_datapoints(std::vector<DataPoint*>& dtp, std::string& datafile)
{
    std::ifstream file(datafile);
    std::string line;
    std::stringstream linechunk;
    std::string param;

    /* get number of data points */
    std::getline(file, param);
    n = stoi(param);
    std::cout << "n: " << n << std::endl;

    /* get dimension size */
    std::getline(file, param);
    d = stoi(param);
    std::cout << "d: " << d << std::endl;

    /* read coordinates of point line by line. Each line corresponds to a new data point */
    for( int i = 0; i < n; i++ )
    {
        std::getline(file, line);
        std::stringstream linechunk(line);
        std::vector<double> coords (d);
        /* reads first column that is not a number */
        std::getline(linechunk, param, ',');

        for( unsigned int j = 0; j < d; j++ )
        {
            std::getline(linechunk, param, ',');
            coords[j] = stod(param);
        }

        DataPoint* dp = new DataPoint(coords);
        dtp.emplace_back(dp);
    }

    std::cout << "Number elements is : " << n << ", dimension size is : " << d << std::endl;
    file.close();
}

/* creates the assignment */
Clusters cluster_clients(const std::vector<int>& a)
{
   Clusters c(n);
    //assert size of a = n;
    for( unsigned int i = 0; i < a.size(); i++ )
        (c)[a[i]].emplace_back(i);

    return c;
}

/* applies Lloyds algorithm with datapoints, centers and assignment */
void Lloyds(std::vector<DataPoint*>& centers, const std::vector<DataPoint*>& points,
    std::vector<double>& dist, std::vector<int>& assignment, Clusters& clusters)
{
    /* compute mass center of all clusters */
    for( unsigned int i = 0; i < centers.size(); i++ )
    {
        if( (clusters)[i].size() == 0 )
            continue;
        else
            (*centers[i]) = compute_centroid(points, (clusters)[i]);
    }

    /* update distance, and closest center */
    for( unsigned int j = 0; j < points.size(); j++ )
    {
        double tmp_dist = infty;
        int tmp_index = centers.size();

        for( unsigned int l = 0; l < centers.size(); l++ )
        {
            if( tmp_dist > distance(*points[j], *centers[l]) )
            {
                tmp_dist = distance(*points[j], *centers[l]);
                tmp_index = l;
            }
        }
        //assert distance[i] != infty;
        //assert tmp_index!> centers.size();
        dist[j] = tmp_dist;
        assignment[j] = tmp_index;
    }
}

/* computes the centroid of a given assignement to a center */
DataPoint compute_centroid(const std::vector<DataPoint*>& points, const std::vector<int>& assigned)
{
    std::vector<double> coords(d, 0.0);
    double size = assigned.size();
    //assert size != 0
    for( const int& i : assigned )
        for( unsigned int j = 0; j < d; j++ )
            coords[j] += (points[i]->get_coordinates()[j]);

    //assert coords has size d

    for( unsigned int r = 0; r < d; r++ )
        coords[r] /= size;

    return DataPoint(coords);
}

/* computes the cost of the given assignement */
double compute_cost(const std::vector<double>& dist)
{
    double value = 0.0;
    //assert dis.size() = n;
    for( unsigned int i = 0; i < dist.size(); i++ )
        value += dist[i];

    return value;
}
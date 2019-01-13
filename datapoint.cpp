#include<vector>
#include<assert.h>
#include<numeric>               //std::inner_product
#include<assert.h>              // assert
#include "datapoint.hpp"


/* default constructor */
DataPoint::DataPoint()
{
}

/* inititializes data point with coordinates, weight and distance */
DataPoint::DataPoint(const std::vector<double>& coords):
    coordinates(coords)
{}

/* copy constructor */
DataPoint::DataPoint(const DataPoint& point):
    coordinates(point.coordinates)
{}

/* destroys the datapoint */
DataPoint::~DataPoint()
{}

/* sets coordinates */
void DataPoint::set_coordinates(const std::vector<double>& buffer)
{
    coordinates = buffer;
}

/* gets coordinates and return const reference to avoid copy */
const std::vector<double>& DataPoint::get_coordinates() const
{
    return this->coordinates;
}

/* prints informations or data point (position/D2-weighting/distance) */
void DataPoint::print() const
{
    std::cout << "coordinates of current data point are:" << std::endl;
    for( int i = 0; i < coordinates.size(); i++ )
    {
        std::cout << coordinates[i] << ", ";
    }
    std::cout << std::endl;
}

/* computes distance between two data points */
double distance(DataPoint const& p, DataPoint const& q)
{
    double dist = 0.0;
    double dist2 = 0.0;
    dist += std::inner_product(p.get_coordinates().begin(), p.get_coordinates().end(),
        p.get_coordinates().begin(), 0.0);
    dist += std::inner_product(q.get_coordinates().begin(), q.get_coordinates().end(),
        q.get_coordinates().begin(), 0.0);
    dist2 += std::inner_product(p.get_coordinates().begin(), p.get_coordinates().end(),
        q.get_coordinates().begin(), dist2);

    assert(dist - (2*dist2) >= 0);
    return dist - (2*dist2);
}
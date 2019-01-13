#include<iostream>
#include<vector>

class DataPoint
{
    public:
        DataPoint();
        DataPoint(const std::vector<double>& coords);
        ~DataPoint();

        /* copy constructor */
        DataPoint(const DataPoint& point);

        /* setters */
        void set_coordinates(const std::vector<double>& buffer);

        /* getters */
        const std::vector<double>& get_coordinates() const;

        /* miscellaneous */
        void print() const;

    private:
        /* vector storing coordinates of point */
        std::vector<double> coordinates;
};

/* computes distance between two data points */
double distance(DataPoint const& p, DataPoint const& q);

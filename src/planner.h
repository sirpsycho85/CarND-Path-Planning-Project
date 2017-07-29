#ifndef PLANNER_H
#define PLANNER_H

#include <vector>
#include <map>

using namespace std;

class Planner {
  public:
    Planner();
    virtual ~Planner();
    void update_telemetry(map<std::string,double> car_data_doubles, map<std::string,auto> car_data_vectors);
  private:
    vector<double> next_x_vals;
    vector<double> next_y_vals;
    
    map<std::string,double> car_data_doubles_;
    map<std::string,auto> car_data_vectors_;
    
    // Previous path data given to the Planner
    vector<double> previous_path_x;
    vector<double> previous_path_y;
    
    // Previous path's end s and d values 
    double end_path_s;
    double end_path_d;
    
    // Sensor Fusion Data, a list of all other cars on the same side of the road.
    vector<double> sensor_fusion;
};

#endif
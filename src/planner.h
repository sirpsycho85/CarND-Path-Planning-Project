#ifndef PLANNER_H
#define PLANNER_H

#include <vector>

using namespace std;

class Planner {
  public:
    Planner();
    virtual ~Planner();
    void update_telemetry(auto j);
  private:
    vector<double> next_x_vals;
    vector<double> next_y_vals;
    
    // Main car's localization Data
    double car_x;
    double car_y;
    double car_s;
    double car_d;
    double car_yaw;
    double car_speed;
    
    // Previous path data given to the Planner
    auto previous_path_x;
    auto previous_path_y;
    
    // Previous path's end s and d values 
    double end_path_s;
    double end_path_d;
    
    // Sensor Fusion Data, a list of all other cars on the same side of the road.
    auto sensor_fusion;
};

#endif
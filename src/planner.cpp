#include "planner.h"
#include "json.hpp"

using namespace std;

Planner::Planner() {}
Planner::~Planner() {}

void Planner::update_telemetry(map<std::string,double> car_data_doubles, map<std::string,auto> car_data_vectors) {

  car_data_doubles_ = car_data_doubles;
  car_data_vectors_ = car_data_vectors;
  
  // // Main car's localization Data
  // car_x = j[1]["x"];
  // car_y = j[1]["y"];
  // car_s = j[1]["s"];
  // car_d = j[1]["d"];
  // car_yaw = j[1]["yaw"];
  // car_speed = j[1]["speed"];

  // // Previous path data given to the Planner
  // previous_path_x = j[1]["previous_path_x"];
  // previous_path_y = j[1]["previous_path_y"];
  
  // // Previous path's end s and d values 
  // end_path_s = j[1]["end_path_s"];
  // end_path_d = j[1]["end_path_d"];

  // // Sensor Fusion Data, a list of all other cars on the same side of the road.
  // sensor_fusion = j[1]["sensor_fusion"];
}

/*
// find next waypoint
int next_waypoint = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);

// make sure the next waypoint is far enough away that the trajectory would approximate the road
next_waypoint++;

// start in current position
vector<double> start = {car_s, car_speed, 0};

// in a particular plan, target 25 m/s or current speed + 5.
double speed_limit = 25;
double target_speed = min(car_speed+5, speed_limit);

// end at the next waypoint, accelerating to the target speed, with no acceleration
vector<double> end = {map_waypoints_s[next_waypoint], target_speed, 0};

// the time to get to the location the solver should use is is distance to waypoint divided by target velocity
double time_target = (map_waypoints_s[next_waypoint]-car_s)/target_speed;

// get the polynomial of the path
vector<double> poly = JMT(start, end, time_target);

// find points along the path in 20ms increments
double t_inc = 0.02;

for(int i = 0; i < 50; i++)
{
  //each next point is built up from the polynomials and the time increment
  double next_s = 0;
  for (int j = 0; j < poly.size(); ++j)
  {
    next_s += poly[j] * pow(t_inc * i, j);
  }

  vector<double> next_xy = getXY(next_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

  next_x_vals.push_back(next_xy[0]);
  next_y_vals.push_back(next_xy[1]);
}

*/
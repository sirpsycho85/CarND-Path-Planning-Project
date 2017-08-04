#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"
#include "spline.h"

/*
TODO: if other car is in same lane and (other_s - car_s)%max_s < safe_distance, set max_v lower

          TODO: Kostas Oreopoulos [1:29 PM] 
In case some else was struggling, The naive solution is always the best. To deal with the car speeding because s is linear and the car is in curves, the following seems to work pretty well.
* once you have your splines and you can convert from s,d to xy in continuous fashion, then do the following
* take for example from the point you are a distance of 100 meters (the next) and divide it in 4-5-6 points you want.
Then for those all those points calculate the XY coodinates and add the distances  (1-2, 2-3.. etc) . You will get a distance Dxy.
If its a straight line you are driving, then Ds (which is just the difference of the first and last points) should be equal to Dxy.
Dxy > Ds, means you are turning, thus travelling a bigger distance. So the speed you said you liked... Dv, which was Ds/dt will in fact to Dxy/dt.
So to compensate you should tell the vehicle to scale down the velocity by that factor (Ds / Dxy) OR , you can apply that directy to the distance you tell it to travel (more accurate )

// TODO: things to try:
  // Spline assumes evenly spaced waypoints but our points are not.......
  // Marcus Erbar suggested: Appending deltas of new_path to previous_path[0]. Then, as a second step, smoothing appended_path with previous_path for a couple of timesteps (I think 20). That finally got rid of the seemingly random noise.
*/

#include <typeinfo>
#include <map>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<double> JMT(vector< double> start, vector <double> end, double T) {
  
  MatrixXd A = MatrixXd(3, 3);
  
  A << T*T*T, T*T*T*T, T*T*T*T*T,
          3*T*T, 4*T*T*T,5*T*T*T*T,
          6*T, 12*T*T, 20*T*T*T;
    
  MatrixXd B = MatrixXd(3,1);     
  B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
          end[1]-(start[1]+start[2]*T),
          end[2]-start[2];
          
  MatrixXd Ai = A.inverse();
  
  MatrixXd C = Ai*B;
  
  vector <double> result = {start[0], start[1], .5*start[2]};
  
  for(int i = 0; i < C.size(); i++)
  {
      result.push_back(C.data()[i]);
  }
  
  return result;
}

int nearest_lane(double d) {
  int lane = (int)((d - 2)/4 + 0.5);
  return lane;
}

int main() {

  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;

  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;
  double offset = 2 + 4 * lane;

  vector<double> map_waypoints_x2;
  vector<double> map_waypoints_y2;

  for (int i = 0; i < map_waypoints_x.size(); ++i)
  {
    map_waypoints_x2.push_back(map_waypoints_x[i] += map_waypoints_dx[i] * offset);
    map_waypoints_y2.push_back(map_waypoints_y[i] += map_waypoints_dy[i] * offset);
  }

  tk::spline sx, sy, sx2, sy2;
  sx.set_points(map_waypoints_s, map_waypoints_x);
  sy.set_points(map_waypoints_s, map_waypoints_y);
  sx2.set_points(map_waypoints_s, map_waypoints_x2);
  sy2.set_points(map_waypoints_s, map_waypoints_y2);

  int num_messages = 0;
  vector<double> prev_s_vals = {};
  double start_v = 0;

  h.onMessage([&max_s,&start_v,&prev_s_vals,&sx,&sy,&sx2,&sy2,&num_messages,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          
        	// Main car's localization Data
        	double car_x = j[1]["x"];
        	double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];

        	// Previous path's end s and d values 
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"]; // [ id, x, y, vx, vy, s, d]

        	json msgJson;


          //----------- SOLUTION -------------


          vector<double> next_x_vals;
          vector<double> next_y_vals;
          double max_a;
          double max_v;
          double traj_t;
          double end_v;
          double end_s;
          vector<double> start;
          vector<double> end;
          vector<double> poly;
          vector<double> next_s_vals;
          vector<double> new_s_vals;
          int num_pts;
          int max_pts_to_reuse;
          int num_pts_used;
          int num_pts_unused;
          int num_pts_to_reuse;
          double t_inc;
          double safe_distance;


        	
          num_pts = 100;

          max_pts_to_reuse = 20;

          t_inc = 0.02;

          num_pts_used = prev_s_vals.size() - previous_path_x.size();

          num_pts_unused = previous_path_x.size();

          num_pts_to_reuse = min(num_pts_unused, max_pts_to_reuse);

          if (num_pts_used >= 2) {
            start_v = (prev_s_vals[num_pts_used - 1] - prev_s_vals[num_pts_used - 2]) / t_inc;
          }

          max_a = 15;
          
          max_v = 30;

          safe_distance = 50;

          for (int car_i = 0; car_i < sensor_fusion.size(); ++car_i)
          {
            int other_lane = nearest_lane(sensor_fusion[car_i][6]);
            
            int current_lane = nearest_lane(car_d);

            int other_s = sensor_fusion[car_i][5];

            double current_distance = fmod((car_s - other_s), max_s);

            if (other_lane == current_lane && current_distance < safe_distance)
            {
              cout << "TOO CLOSE!" << endl;
            }
          }
          
          traj_t = t_inc * num_pts;

          end_v = min(start_v + max_a * traj_t, max_v);
          
          end_s = car_s + (start_v + end_v)/2 * traj_t;
          
          start = {car_s, start_v, 0};
          
          end = {end_s, end_v, 0};

          poly = JMT(start, end, traj_t);

          next_s_vals = {};

          new_s_vals = {};

          for(int i = 0; i < num_pts; i++) {

            double new_s;

            new_s = 0;

            for (int j = 0; j < poly.size(); ++j) {
              new_s += poly[j] * pow(t_inc * i, j);
            }

            new_s_vals.push_back(new_s);

            if (i == 0 && num_pts_to_reuse > 0) {
              next_s_vals.push_back(prev_s_vals[i + num_pts_used]);
            }

            else if (i == 0 && num_pts_to_reuse <= 0) {
              next_s_vals.push_back(new_s_vals[i]);
            }

            else {
              next_s_vals.push_back(next_s_vals[i-1] + (new_s_vals[i] - new_s_vals[i-1]));
            }
          }

          for(int i = 0; i < num_pts; i++) {
            double next_x;
            double next_y;

            next_x = sx2(next_s_vals[i]);
            next_y = sy2(next_s_vals[i]);

            next_x_vals.push_back(next_x);
            next_y_vals.push_back(next_y);
          }

          prev_s_vals = next_s_vals;

        	msgJson["next_x"] = next_x_vals;
        	msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

        	//this_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      }
      else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
    num_messages++;
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

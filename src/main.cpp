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

// #include "planner.h"

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

/*
vector<double> get_pos_after_reuse(vector<double> car_current, auto previous_path_x, auto previous_path_y, int num_pts_to_reuse) {
  
  double car_x = car_current[0];
  double car_y = car_current[1];
  double car_yaw = car_current[2];

  double pos_x, pos_y, angle;

  if(num_pts_to_reuse == 0)
  {
    pos_x = car_x;
    pos_y = car_y;
    angle = deg2rad(car_yaw);
  }
  else if (num_pts_to_reuse == 1)
  {
    pos_x = previous_path_x[num_pts_to_reuse-1];
    pos_y = previous_path_y[num_pts_to_reuse-1];
    angle = atan2(pos_y-car_y,pos_x-car_x);
  }
  else
  {
    pos_x = previous_path_x[num_pts_to_reuse-1];
    pos_y = previous_path_y[num_pts_to_reuse-1];
    double pos_x2 = previous_path_x[num_pts_to_reuse-2];
    double pos_y2 = previous_path_y[num_pts_to_reuse-2];
    angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
  }

  vector<double> pos_after_reuse;

  pos_after_reuse.push_back(pos_x);
  pos_after_reuse.push_back(pos_y);
  pos_after_reuse.push_back(angle);

  return pos_after_reuse;
}
*/

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


  // double map_x = 784.6001;
  // double map_y = 1135.571;
  // double map_s = 90.4504146575928;
  // double d = 0;

  // vector<double> gotten_XY = getXY(map_s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  // cout << gotten_XY[0]<< " " << gotten_XY[1] << endl;
  // cout << map_x<< " " << map_y << endl;

  // vector<double> gotten_SD = getFrenet(map_x, map_y, 0, map_waypoints_x, map_waypoints_y);
  // cout << gotten_SD[0] << endl;
  // cout << map_s << endl;

  // exit(EXIT_FAILURE);



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

  vector<double> cached_s = {};
  double start_v = 0;

  h.onMessage([&start_v,&cached_s,&sx,&sy,&sx2,&sy2,&num_messages,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
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
        	auto sensor_fusion = j[1]["sensor_fusion"];

        	json msgJson;

          //----------- SOLUTION -------------

          vector<double> next_x_vals;
          vector<double> next_y_vals;
        	
          int num_pts = 100;
          int max_pts_to_reuse = 0;

          double t_inc = 0.02;
          // int num_pts_used = previous_path_x.size();
          int num_pts_used = cached_s.size() - previous_path_x.size(); // 0 - 0 = 0
          int num_pts_unused = cached_s.size() - num_pts_used;

          int num_pts_to_reuse = min(num_pts_unused, max_pts_to_reuse); // min(0-0,3) = 0
          
          // for(int i = 0; i < num_pts_to_reuse; i++)
          // {
          //   next_x_vals.push_back(previous_path_x[i]);
          //   next_y_vals.push_back(previous_path_y[i]);
          // }

          int start_index = num_pts_used - 1; // -1
          int end_index = start_index + num_pts_to_reuse; // -1
          vector<double>temp_cached_s = {}; // {}

          for (int i = start_index; i < end_index; ++i) // -1 !< -1, loop does not run
          {
            temp_cached_s.push_back(cached_s[i]);

            double next_x = sx2(cached_s[i]);
            double next_y = sy2(cached_s[i]);
            next_x_vals.push_back(next_x);
            next_y_vals.push_back(next_y);
          }

          if (num_pts_to_reuse >= 2)
          {
            start_v = (cached_s[num_pts_to_reuse - 1] - cached_s[num_pts_to_reuse - 2]) / t_inc;
          }

          cached_s = temp_cached_s; // {}

          double pos_x;
          double pos_y;
          double angle;
          double pos_s;

          if(num_pts_to_reuse == 0) // this
          {
            pos_x = car_x;
            pos_y = car_y;
            angle = deg2rad(car_yaw);

            pos_s = car_s;
          }
          else if (num_pts_to_reuse == 1)
          {
            pos_x = previous_path_x[num_pts_to_reuse-1];
            pos_y = previous_path_y[num_pts_to_reuse-1];
            angle = atan2(pos_y-car_y,pos_x-car_x);

            pos_s = cached_s[cached_s.size() - 1];
          }
          else
          {
            pos_x = previous_path_x[num_pts_to_reuse-1];
            pos_y = previous_path_y[num_pts_to_reuse-1];
            double pos_x2 = previous_path_x[num_pts_to_reuse-2];
            double pos_y2 = previous_path_y[num_pts_to_reuse-2];
            angle = atan2(pos_y-pos_y2,pos_x-pos_x2);

            pos_s = cached_s[cached_s.size() - 1];
          }

          // TODO: things to try:
          // cache s values, because car_s is not equal to pos_s derived from pos_x, pos_y - even with no reuse
          // spline to pick end, JMT to generate path, spline to convert s -> x,y
          // ----need to cache speed and calculate it myself for the target
          // generate new path from car current path every step and average
          // average old values with the new ones
          // Marcus Erbar suggested: Appending deltas of new_path to previous_path[0]. Then, as a second step, smoothing appended_path with previous_path for a couple of timesteps (I think 20). That finally got rid of the seemingly random noise.

          int next_waypoint = NextWaypoint(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
            
          next_waypoint++;

          // vector<double> pos_sd = getFrenet(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);

          // double pos_s = pos_sd[0];
          // double pos_s = car_s;

          double max_a = 2;
          double max_v = 10;
          double traj_t = t_inc * num_pts;

          
          // double end_v = min(car_speed + max_a * traj_t, max_v);
          // double end_v = 5;
          double end_v = min(start_v + max_a * traj_t, max_v);

          // TODO:
          /*
          Kostas Oreopoulos [1:29 PM] 
In case some else was struggling, The naive solution is always the best. To deal with the car speeding because s is linear and the car is in curves, the following seems to work pretty well.
* once you have your splines and you can convert from s,d to xy in continuous fashion, then do the following
* take for example from the point you are a distance of 100 meters (the next) and divide it in 4-5-6 points you want.
Then for those all those points calculate the XY coodinates and add the distances  (1-2, 2-3.. etc) . You will get a distance Dxy.
If its a straight line you are driving, then Ds (which is just the difference of the first and last points) should be equal to Dxy.
Dxy > Ds, means you are turning, thus travelling a bigger distance. So the speed you said you liked... Dv, which was Ds/dt will in fact to Dxy/dt.
So to compensate you should tell the vehicle to scale down the velocity by that factor (Ds / Dxy) OR , you can apply that directy to the distance you tell it to travel (more accurate )
          */
          

          // double traj_t = (map_waypoints_s[next_waypoint]-car_s)/end_v;
          // double end_s = pos_s + 0.5 * max_a * traj_t * traj_t;
          double end_s = pos_s + (start_v + end_v)/2 * traj_t;
          

          // vector<double> start = {pos_s, car_speed, 0};
          vector<double> start = {pos_s, start_v, 0};
          vector<double> end = {end_s, end_v, 0};
          // vector<double> end = {map_waypoints_s[next_waypoint], end_v, 0};

          double last_v = min(start_v + max_a * t_inc * num_pts_used, max_v);
          start_v = last_v;

          vector<double> poly = JMT(start, end, traj_t);

          // double last_s = 0;
          // for (int j = 0; j < poly.size(); ++j)
          // {
          //   last_s += poly[j] * pow(t_inc * num_pts, j);
          // }

          // cout << car_s << "\t" << pos_s << "\t" << end_s << "\t" << last_s << endl;

          // exit(EXIT_FAILURE);

          // double dist_inc = 0.1;

          for(int i = 0; i < num_pts - num_pts_to_reuse; i++)
          {
            // ----- JMT APPROACH
            double next_s = 0;
            for (int j = 0; j < poly.size(); ++j)
            {
              next_s += poly[j] * pow(t_inc * i, j);
            }
            cached_s.push_back(next_s);

            // vector<double> next_xy = getXY(next_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            // double next_x = next_xy[0];
            // double next_y = next_xy[1];


            // -----SPLINE ONLY APPROACH
            
            // double next_s = pos_s + dist_inc * (i); //TODO: why not +1?
            
            double next_x = sx2(next_s);
            double next_y = sy2(next_s);
            

            next_x_vals.push_back(next_x);
            next_y_vals.push_back(next_y);
          }

          // for (int i = 0; i < next_x_vals.size(); ++i)
          // {
          //   cout << next_x_vals[i] << ", ";
          // }
          // cout<<endl;

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

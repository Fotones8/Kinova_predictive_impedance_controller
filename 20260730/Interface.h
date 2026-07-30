//-------------------------------------------------
// Header file for Interface.cpp
//-------------------------------------------------
#include <Trajectory.h>

#ifndef EXPERIMENTALCODES_CHANNEL_H
#define EXPERIMENTALCODES_CHANNEL_H


// struct rail_noise with noise of rotation center, angle and the lengths to the channel
struct rail_noise {
    int translation;    // rotation center noise
    double ori_N;   // rotation angle noise
    double y_N;     // length noise to side 1 of the channel
    double z_N;     // length noise to side 2 of the channel

    rail_noise(double translation, double ori_N, double y_N, double z_N) : translation(translation),
                 ori_N(ori_N), y_N(y_N), z_N(z_N) {}
};


//----------------------------------------------
// Functions used in the main codes
//----------------------------------------------
// draw_interface: draw the information in the interface
// ini_noise: initialize the visual noise
// ini_draw: initialize the drawing
// start_display: displayed information at start
// end_display: displayed information at end
// draw_results: draw the final results
//----------------------------------------------
namespace Interface {
    void draw_interface(double block_y, double block_z, double block_angle, vector<Point>& traj, double noise_mean,
                   double noise_dev, pair<int, int>& visual_wd, int display, int count, int c_f, int inverse,
                   int& task, VectorXd& interface_center);
    void ini_noise(double& mean, double& stddev);
    void ini_draw();
    void start_display(int type, int second);
    void end_display();
    void draw_results (vector<Point>& traj, vector<Point>& results, std::vector<Point>& Curve1, std::vector<Point>& Curve2, int inverse);
};


#endif //EXPERIMENTALCODES_CHANNEL_H

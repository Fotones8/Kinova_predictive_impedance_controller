//----------------------------------------------------------
// Function of visual interface
//----------------------------------------------------------
// Description: function to display information in the visual interface
// Copyright: Yihan Liu 2024
//----------------------------------------------------------

#define _USE_MATH_DEFINES
#include "Interface.h"
#include <graphics.h>
#include <math.h>
#include <vector>
#include <Trajectory.h>
#include <noise.h>
#include <random>

using namespace noise;
using namespace std;

//-------------display parameters------------------
int WIDTH = 1900;                  // width of the display window (unit: pixel)
int HEIGHT = 1000;                 // height of the display window (unit: pixel)
double gain = 1.1;                 // gain of amplification factor of the real width/height
double real_width = gain * 0.6;    // Real width of moving area (unit: m)
double real_height = gain * 0.32;  // Real height of moving area (unit: m)
double central_y = 0;
double central_z = 0;
int change_center = 0;
double normal_centeral_y = 0.5 * real_width;      // normalized y axis of displayed starting point (unit: m)
double normal_centeral_z = 0.5 * real_height;  // normalized z axis of displayed starting point (unit: m)
double rail_wd = 0.02;             // width of the displayed rail (unit: m)
int display_rail = 0;              // switch to display the green rail (1:ON 0:OFF)

double real_block_wd = 0.055;      // real width of the displayed target block (unit: m)
double real_block_hg = 0.02;      // real height of the displayed target block (unit: m)

int num_block = 8;                  // number of displayed rail (including the noised rails)

double refresh_duration = 0.02;    // refresh duration of the visual noise
int refresh_rail = 0;              // index of the refresh rail

vector<rail_noise> track_noise;    // buffer to save the rail noise


//----------------------------------------------------------
// Function to initialize the visual noise
//----------------------------------------------------------
// 2 inputs:
// mean: mean of the visual noise
// stddev: standard deviation of the visual noise
//----------------------------------------------------------
void Interface :: ini_noise(double& mean, double& stddev){
    track_noise.clear();
    int translation = 0;
    double ori_noise = 0;
    double y_noise = 0;
    double z_noise = 0;

    // Initialize the visual noise for the first step
    for (int c = 0; c < num_block; c++) {
        translation = get_rand_index(0, 10);
        ori_noise = 10 * Gaussian_noise(mean, stddev);
        y_noise = Gaussian_noise(mean, stddev);
        z_noise = Gaussian_noise(mean, stddev);
        track_noise.emplace_back(translation, ori_noise, y_noise, z_noise);
    }
}


//----------------------------------------------------------
// Function to generate the rotated curve
//----------------------------------------------------------
// 7 inputs:
// input_points: desired input trajectory
// distance: distance between the trajectory and channel
// noise_mean: mean of visual noise
// noise_dev: standard deviation of visual noise
// visual_wd: visual width of trajectory with <start index, end index>
// refresh: index of the refresh rail
// drawn_rail: index of current drawn rail
//
// 2 outputs:
// <parallel_curve_1, parallel_curve_2>: 2 generated parallel curves
//----------------------------------------------------------
Point gen_noised_block(
        vector<Point>& input_points, double distance,
        double noise_mean, double noise_dev,
        pair<int, int>& visual_wd, int refresh, int drawn_rail) {

    // If refresh is needed & current rail is the desired refreshing rail
    if ((refresh == 1) && (drawn_rail == refresh_rail)) {
        track_noise[refresh_rail].translation = get_rand_index(visual_wd.first - 50, visual_wd.first + 50);
        track_noise[refresh_rail].ori_N = 10 * Gaussian_noise(noise_mean, noise_dev);
        track_noise[refresh_rail].y_N = Gaussian_noise(noise_mean, noise_dev);
        track_noise[refresh_rail].z_N = Gaussian_noise(noise_mean, noise_dev);
    }

    Point noised_block;

    noised_block.y = input_points[visual_wd.first].y + track_noise[drawn_rail].y_N;
    noised_block.z = input_points[visual_wd.first].z + track_noise[drawn_rail].z_N;
    noised_block.ori = input_points[visual_wd.first].ori + track_noise[drawn_rail].ori_N;

    return noised_block;
}


//----------------------------------------------------------
// Function to map the real position into the displayed position in the interface
//----------------------------------------------------------
// 4 inputs:
// y: current y-axis position (unit: m)
// z: current z-axis position (unit: m)
// start_y: y-axis position of start (unit: m)
// start_z: z-axis position of start (unit: m)
//
// 2 outputs:
// <updated_x, updated_y>: x and y axis position in the interface (unit: pixel)
//----------------------------------------------------------
pair<int, int> mapping(double y, double z, double& start_y, double& start_z){
    int updated_x = (y - central_y + normal_centeral_y)*(WIDTH/real_width);
    int updated_y = (HEIGHT - ((z - central_z + normal_centeral_z)*(HEIGHT/real_height)));

    return make_pair(updated_x, updated_y);
}

//----------------------------------------------------------
// Function to convert the sequence of the points in curve (unit: m)
// into the sequence of displayed curve POINTs in the interface
//----------------------------------------------------------
// 3 inputs:
// curve: sequence of the Points in curve
// start_y: y-axis position of start (unit: m)
// start_z: z-axis position of start (unit: m)
//
// 1 outputs:
// POINT_array: sequence of displayed curve POINTs
//----------------------------------------------------------
vector<POINT> convert_2_POINT(const vector<Point>& curve, double& start_y, double& start_z) {
    // Define the array to save the converted POINT
    vector<POINT> POINT_array(curve.size());

    for (size_t i = 0; i < curve.size(); ++i) {
        // Scaling and shifting
        pair<int, int> display_xy = mapping(curve[i].y, curve[i].z, start_y, start_z);
        POINT_array[i].x = display_xy.first;
        POINT_array[i].y = display_xy.second;
    }
    return POINT_array;
}


//----------------------------------------------------------
// Function to draw the curve in the interface
//----------------------------------------------------------
// 4 inputs:
// curve: Sequence of the Points in curve
// color: the color of drawing curve
// start_y: y-axis position of start (unit: m)
// start_z: z-axis position of start (unit: m)
//----------------------------------------------------------
void draw_curve(const vector<Point>& curve, COLORREF color, double& start_y, double& start_z) {
    // Set color and line style
    setlinecolor(color);
    setlinestyle(PS_SOLID, 3);

    // Convert real point (in m) into POINT (in pixels) & draw in the interface
    vector<POINT> point_array = convert_2_POINT(curve, start_y, start_z);
    polyline(point_array.data(), static_cast<int>(point_array.size()));
}


//----------------------------------------------------------
// Function to draw the current block
//----------------------------------------------------------
// 3 inputs:
// x: x-axis position of the current block (in the interface & unit: pixel)
// y: y-axis position of the current block (in the interface & unit: pixel)
// angle: current orientation of block (in the interface & unit: rad)
//----------------------------------------------------------
void draw_current_block(long x, long y, double angle) {
    // Set the color and line style
    setfillcolor(RED);
    setlinecolor(RED);
    setlinestyle(PS_SOLID, 1);

    // Define the full block dimensions
    const long full_width = real_block_wd * (WIDTH / real_width);   // Total width
    const long full_height = real_block_hg * (HEIGHT / real_height); // Total height

    // Define the four corners of the diamond (square rotated 45 degrees)
    POINT diamond[4] = {
            {x, y - full_height / 2},  // Top
            {x + full_width / 2, y},   // Right
            {x, y + full_height / 2},  // Bottom
            {x - full_width / 2, y}    // Left
    };

    // Rotate the diamond around (x, y)
    POINT rotatedDiamond[4];
    for (int i = 0; i < 4; i++) {
        long dx = diamond[i].x - x;
        long dy = diamond[i].y - y;
        rotatedDiamond[i].x = x + (long)(dx * cos(angle) - dy * sin(angle));
        rotatedDiamond[i].y = y + (long)(dx * sin(angle) + dy * cos(angle));
    }

    // Draw the rotated diamond
    fillpolygon(rotatedDiamond, 4);
}



//----------------------------------------------------------
// Function to draw the target block (desired state)
//----------------------------------------------------------
// 3 inputs:
// x: x-axis position of the current block (in the interface & unit: pixel)
// y: y-axis position of the current block (in the interface & unit: pixel)
// angle: current orientation of block (in the interface & unit: rad)
//----------------------------------------------------------
void draw_target_block(long x, long y, double angle) {
    // Set the color and line style
    setfillcolor(RGB(144, 170, 220));
    setlinecolor(RGB(144, 170, 220));
    setlinestyle(PS_SOLID, 1);

    // Define the full block dimensions
    const long full_width = real_block_wd * (WIDTH / real_width);   // Total width
    const long full_height = real_block_hg * (HEIGHT / real_height); // Total height

    // Define the four corners of the diamond (square rotated 45 degrees)
    POINT diamond[4] = {
            {x, y - full_height / 2},  // Top
            {x + full_width / 2, y},   // Right
            {x, y + full_height / 2},  // Bottom
            {x - full_width / 2, y}    // Left
    };

    // Rotate the diamond around (x, y)
    POINT rotatedDiamond[4];
    for (int i = 0; i < 4; i++) {
        long dx = diamond[i].x - x;
        long dy = diamond[i].y - y;
        rotatedDiamond[i].x = x + (long)(dx * cos(angle) - dy * sin(angle));
        rotatedDiamond[i].y = y + (long)(dx * sin(angle) + dy * cos(angle));
    }

    // Draw the rotated diamond
    fillpolygon(rotatedDiamond, 4);
}


//----------------------------------------------------------
// Function to display the start information
//----------------------------------------------------------
// 2 inputs:
// type: type of displayed information
//       1: display "Ready!"
//       2: display the count backward number
//       3: display "Start"
// second: the clock time
//----------------------------------------------------------
void Interface::start_display(int type, int second){
    // Define the text parameters
    int text_height = 240;
    int text_width = 0;
    int x = 0;
    int y = 0;
    settextcolor(RED);
    settextstyle(text_height, 0, _T("Arial"));

    // Display different types of information with the input
    char time[10];
    if (type == 1){
        text_width = textwidth("Ready!");
        x = (WIDTH - text_width) / 2;
        y = (HEIGHT - text_height) / 2;
        outtextxy(x, y, _T("Ready!"));
    }
    else if (type == 2){
        _itoa_s(second, time, 10);  // Convert 'num' to string and store in 'str'
        text_width = textwidth(time);
        x = (WIDTH - text_width) / 2;
        y = (HEIGHT - text_height) / 2;
        outtextxy(x, y,time);
    }
    else{
        text_width = textwidth("Start!");
        x = (WIDTH - text_width) / 2;
        y = (HEIGHT - text_height) / 2;
        outtextxy(x, y, _T("Start!"));
    }
}


//----------------------------------------------------------
// Function to display the end information (Display "Finish!" in the interface)
//----------------------------------------------------------
void Interface::end_display(){
    // Set the displayed words parameters
    int text_height = 240;
    settextcolor(RED);
    settextstyle(text_height, 0, _T("Arial"));
    int text_width = textwidth("Finish!");
    int x = (WIDTH - text_width) / 2;
    int y = (HEIGHT - text_height) / 2;

    // Display the information
    outtextxy(x, y, _T("Finish!"));
}


//----------------------------------------------------------
// Function to initialize the display setup
//----------------------------------------------------------
void Interface::ini_draw(){
    // Clear the screen and redraw
    // Initialize graphics window
    initgraph(WIDTH, HEIGHT);
    setbkcolor(WHITE);
    cleardevice();
    // Start batch drawing
    BeginBatchDraw();
}


//----------------------------------------------------------
// Function to draw the interface on the screen
//----------------------------------------------------------
// 12 inputs:
// block_y, block_z & block_angle: block y and z-axis position and orientation angle
// traj: desired designed trajectory
// noise_mean & noise_dev: mean and standard deviation of visual noise
// visual_wd: index width of the visual display
// display: type of display information
// count: count number shown in the interface
// c_f: control frequency
// inverse: if the movement of end effector in thr inverse direction
// task: task type (trajectory tracking or fix in a point)
//----------------------------------------------------------
void Interface::draw_interface(double block_y, double block_z, double block_angle, vector<Point>& traj, double noise_mean,
                        double noise_dev, pair<int, int>& visual_wd, int display, int count, int c_f, int inverse,
                        int& task, VectorXd& interface_center) {

//    // Check if the task is trajectory tracking / fixing in a point
//    if(task == 1){
//        auto result = minmax_element(traj.begin(), traj.end(),
//                                          [](const Point& a, const Point& b) {
//                                              return a.z < b.z; // Compare based on the z value
//                                          });
//        normal_start_z = ((traj[0].z/(result.first->z + result.second->z))) * real_height;
//        normal_start_y = 0.05;
//    }
//    else{
//        normal_start_y = 0.5*(real_width);
//        normal_start_z = 0.5*(real_height);
//    }

    // Clear the display
    cleardevice();
    // Generate the left and right parallel curves
    vector<Point> target_blocks(num_block);

    // Set refresh rail to be 0
    int refresh = 0;

    // Check if the visual noise needs to be refreshed
    if((visual_wd.first%(int)(refresh_duration*c_f) == 0) && (visual_wd.first != 0)){
        if(refresh_rail >= (num_block - 1)){
            refresh_rail = 0;
        }
        else{
            refresh_rail += 1;
        }
        refresh = 1;
    }

    // Generate the 5 rotated channel with boundary curves
    for (int n = 0; n < num_block; n++) {
        Point noised_target = gen_noised_block(traj, rail_wd/2, noise_mean, noise_dev,
                                                   visual_wd, refresh, n);

        target_blocks[n] = noised_target;
    }

    // Get the y and z-axis position of start from trajectory
    double start_y = 0;
    double start_z = 0;

    if (inverse == 1){
        start_y = traj[traj.size() - 1].y;
        start_z = traj[traj.size() - 1].z;
    }
    else {
        start_y = traj[0].y;
        start_z = traj[0].z;
    }

//    if (change_center == 0){
//        central_y = traj[0].y;
//        central_z = traj[0].z;
//        change_center = 1;
//    }

        central_y = interface_center[1];
        central_z = interface_center[2];

    // Draw the generated channel and target block
    if (visual_wd.second == 0){ // if the block (end effector) is at the start
        pair<int, int> target_xy = mapping(traj[0].y, traj[0].z, start_y, start_z);
        draw_target_block(target_xy.first, target_xy.second, -traj[0].ori);
    }
    else {
        for (int k = 0; k < target_blocks.size(); k++) {
                pair<int, int> target_xy = mapping(target_blocks[k].y, target_blocks[k].z, start_y, start_z);
                draw_target_block(target_xy.first, target_xy.second, -target_blocks[k].ori);
        }
    }

    // Draw the current block states in the interface
    pair<int,int> block_xy = mapping(block_y, block_z, start_y, start_z);
    draw_current_block(block_xy.first, block_xy.second, block_angle);

    // display the information at the start and end
    if ((display == 1) || (display == 2) || (display == 3)){
        start_display(display, count);
    }
    else if(display == 4){
        end_display();
    }

    FlushBatchDraw();
}

// Function to compute the average errors for y, z, and ori
void computeAverageErrors(const vector<Point>& traj, const vector<Point>& results) {
    // Use the smaller size to prevent out-of-bounds access
    int n = min(traj.size(), results.size());

    if (n == 0) {
        cerr << "Error: One or both input vectors are empty!" << std::endl;
        return;
    }

    double totalErrorY = 0.0;
    double totalErrorZ = 0.0;
    double totalErrorOri = 0.0;

//    cout << n << endl;

    for (int i = 0; i < n; ++i) {
        totalErrorY += abs(traj[i].y - results[i].y);
        totalErrorZ += abs(traj[i].z - results[i].z);
        totalErrorOri += abs(traj[i].ori + results[i].ori);
    }

    double avgErrorY = totalErrorY / n;
    double avgErrorZ = totalErrorZ / n;
    double avgErrorOri = totalErrorOri / n;

    cout << "Average Error in y: " << avgErrorY << endl;
    cout << "Average Error in z: " << avgErrorZ << endl;
    cout << "Average Error in ori: " << avgErrorOri << endl;
}

//----------------------------------------------------------
// Function to display the human-created and desired trajectories
//----------------------------------------------------------
// 5 inputs:
// traj: desired designed trajectory
// results: tracking results created by human subjects
// Curve1 & Curve2: two curves to form the channel for the desired trajectory
// inverse: if the moving direction is forward or inverse
//----------------------------------------------------------
void Interface::draw_results(vector<Point>& traj, vector<Point>& results, vector<Point>& Curve1, vector<Point>& Curve2, int inverse) {
    string input;

    // Display the desired and subject-created results until 'y' is pressed
    while (true) {
//        cleardevice();
//
//        double start_y = 0;
//        double start_z = 0;
//
//        if (inverse == 1){
//            start_y = traj[traj.size() - 1].y;
//            start_z = traj[traj.size() - 1].z;
//        }
//        else {
//            start_y = traj[0].y;
//            start_z = traj[0].z;
//        }
//
//        // Draw the generated left and right curves
//        draw_curve(traj, GREEN, start_y, start_z);
//        draw_curve(results, RED, start_y, start_z);

        computeAverageErrors(traj, results);

//        FlushBatchDraw();

        // Display the results until 'y' is pressed
        cout << "Press 'y' to exit the results display: ";
        cin >> input;

        if (input == "y") {  // Check if the input is 'y'
            cout << "Exiting..." << endl;
            break;
        }
    }
}
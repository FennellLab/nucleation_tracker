/*
   Program to calculate the hbond rings in aqueous environments 
   */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

// some global parameters
double distance_tol = 3.5; // hbond distance tolerance (in Ångstrom)
double distance_tol2 = distance_tol*distance_tol; // hbond square distance tolerance
double angle_tol = 30; // hbond angle tolerance (in degrees)
double angle_tol_rad = angle_tol * 3.1415926536 / 180.0; // hbond angle tolerance in radians
double OH_length = 0.926; // OH bond length (in Ångstrom)
double OH_length_i = 1.0 / OH_length; // OH bond length (in Ångstrom)
double rad_to_degree = 180.0 * 3.1415926536; // convert radians to degrees
double ring3_dot_size = 0.5; // nucleation dot size
double ring4_dot_size = 0.6; // nucleation dot size
double ring5_dot_size = 0.7; // nucleation dot size
double ring6_dot_size = 0.8; // nucleation dot size
double ring7_dot_size = 0.9; // nucleation dot size
double scale_factor = 20; // scale the pixels!
double buffer_size = 2; // the imaging buffer size
double slab_thickness = 40; // 0.5*slab thickness
double transparency = 0.6; // the transparency of ring dots
double outline = 0.04;
double axes_width = 0.08;

class Calcs {
    public:
        Calcs();
        ~Calcs();
        int IsWat3HBond(double (&pos1)[9], double (&pos2)[9], double (&box)[3]);
    private:
        bool printFlag;
};

Calcs::Calcs() {
    printFlag = true;
}

Calcs::~Calcs() {
}

int Calcs::IsWat3HBond(double (&pos1)[9], double (&pos2)[9], double (&box)[3]){
    // a function to identify if there is an HBond between 3 point water models
    int isHBond = 0;
    double r, r2, ri, xVal, yVal, zVal, tempAngle;
    double vec1[3], vec2[3];

    // the 9 positions are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z)
    // first, we do the O-O vector
    xVal = pos2[0]-pos1[0];
    yVal = pos2[1]-pos1[1];
    zVal = pos2[2]-pos1[2];

    // do minimum image filter...
    xVal -= box[0]*copysign(1.0,xVal)*floor(fabs(xVal/box[0])+0.5);
    yVal -= box[1]*copysign(1.0,yVal)*floor(fabs(yVal/box[1])+0.5);
    zVal -= box[2]*copysign(1.0,zVal)*floor(fabs(zVal/box[2])+0.5);

    // here, we take care of normalization 
    r2 = xVal*xVal + yVal*yVal + zVal*zVal;
    r = sqrt(r2);
    ri = 1.0 / r;

    vec1[0] = xVal * ri;
    vec1[1] = yVal * ri;
    vec1[2] = zVal * ri;

    // the O-H1
    xVal = pos1[3]-pos1[0];
    yVal = pos1[4]-pos1[1];
    zVal = pos1[5]-pos1[2];
    // do minimum image filter...
    xVal -= box[0]*copysign(1.0,xVal)*floor(fabs(xVal/box[0])+0.5);
    yVal -= box[1]*copysign(1.0,yVal)*floor(fabs(yVal/box[1])+0.5);
    zVal -= box[2]*copysign(1.0,zVal)*floor(fabs(zVal/box[2])+0.5);

    // here, we take care of normalization 
    r2 = xVal*xVal + yVal*yVal + zVal*zVal;
    r = sqrt(r2);
    ri = 1.0 / r;

    vec2[0] = xVal * ri;
    vec2[1] = yVal * ri;
    vec2[2] = zVal * ri;
    // vec2[0] = xVal * OH_length_i;
    // vec2[1] = yVal * OH_length_i;
    // vec2[2] = zVal * OH_length_i;

    tempAngle = acos((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));  

    if (tempAngle <= angle_tol_rad){
        isHBond = 1;
    } else {
        // continue down the rabbit hole

        // the O-H2
        xVal = pos1[6]-pos1[0];
        yVal = pos1[7]-pos1[1];
        zVal = pos1[8]-pos1[2];
        // do minimum image filter...
        xVal -= box[0]*copysign(1.0,xVal)*floor(fabs(xVal/box[0])+0.5);
        yVal -= box[1]*copysign(1.0,yVal)*floor(fabs(yVal/box[1])+0.5);
        zVal -= box[2]*copysign(1.0,zVal)*floor(fabs(zVal/box[2])+0.5);

        // here, we take care of normalization 
        r2 = xVal*xVal + yVal*yVal + zVal*zVal;
        r = sqrt(r2);
        ri = 1.0 / r;

        vec2[0] = xVal * ri;
        vec2[1] = yVal * ri;
        vec2[2] = zVal * ri;
        //    vec2[0] = xVal * OH_length_i;
        //    vec2[1] = yVal * OH_length_i;
        //    vec2[2] = zVal * OH_length_i;

        tempAngle = acos((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));  

        if (tempAngle <= angle_tol_rad){
            isHBond = 1;
        } else {
            // continue further down the rabbit hole

            // flip the O-O vector
            vec1[0] *= -1;
            vec1[1] *= -1;
            vec1[2] *= -1;

            // the O-H1
            xVal = pos2[3]-pos2[0];
            yVal = pos2[4]-pos2[1];
            zVal = pos2[5]-pos2[2];
            // do minimum image filter...
            xVal -= box[0]*copysign(1.0,xVal)*floor(fabs(xVal/box[0])+0.5);
            yVal -= box[1]*copysign(1.0,yVal)*floor(fabs(yVal/box[1])+0.5);
            zVal -= box[2]*copysign(1.0,zVal)*floor(fabs(zVal/box[2])+0.5);

            // here, we take care of normalization 
            r2 = xVal*xVal + yVal*yVal + zVal*zVal;
            r = sqrt(r2);
            ri = 1.0 / r;

            vec2[0] = xVal * ri;
            vec2[1] = yVal * ri;
            vec2[2] = zVal * ri;
            //        vec2[0] = xVal * OH_length_i;
            //        vec2[1] = yVal * OH_length_i;
            //        vec2[2] = zVal * OH_length_i;

            tempAngle = acos((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));  
            
            if (tempAngle <= angle_tol_rad){
                isHBond = 1;
            } else {
                // last level... I promise...

                // the O-H2
                xVal = pos2[6]-pos2[0];
                yVal = pos2[7]-pos2[1];
                zVal = pos2[8]-pos2[2];
                // do minimum image filter...
                xVal -= box[0]*copysign(1.0,xVal)*floor(fabs(xVal/box[0])+0.5);
                yVal -= box[1]*copysign(1.0,yVal)*floor(fabs(yVal/box[1])+0.5);
                zVal -= box[2]*copysign(1.0,zVal)*floor(fabs(zVal/box[2])+0.5);

                // here, we take care of normalization 
                r2 = xVal*xVal + yVal*yVal + zVal*zVal;
                r = sqrt(r2);
                ri = 1.0 / r;

                vec2[0] = xVal * ri;
                vec2[1] = yVal * ri;
                vec2[2] = zVal * ri;
                //            vec2[0] = xVal * OH_length_i;
                //            vec2[1] = yVal * OH_length_i;
                //            vec2[2] = zVal * OH_length_i;

                tempAngle = acos((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));  

                if (tempAngle <= angle_tol_rad){
                    isHBond = 1;
                }
            }
        }
    }

    return isHBond;
}


class PovObjects {
    public:
        PovObjects();
        ~PovObjects();
        void printHeader(ofstream& outputFile, double (&boxLength)[3]);
        void printHeader2(ofstream& outputFile);
        void printRing3(ofstream& outputFile, double (&boxLength)[3]);
        void printRing4(ofstream& outputFile, double (&boxLength)[3]);
        void printRing5(ofstream& outputFile, double (&boxLength)[3]);
        void printRing6(ofstream& outputFile, double (&boxLength)[3]);
        void printRing7(ofstream& outputFile, double (&boxLength)[3]);
        void printBar(ofstream& outputFile);
        void printAxes(ofstream& outputFile);
    private:
};

PovObjects::PovObjects(){
}

PovObjects::~PovObjects(){
}

void PovObjects::printHeader(ofstream& outputFile, double (&boxLength)[3]){
    outputFile <<
        "  //**************************************\n"
        "  // generated by nucleation_tracker\n"
        "  //**************************************\n"
        "\n"
        "  //**************************************\n"
        "  // Lights, camera, resolution!\n"
        "  //**************************************\n"
        "  global_settings{ max_trace_level 100 }\n"
        "\n"
        "#declare Ratio = " << boxLength[0]/boxLength[1] << ";\n"
        "#declare zoom = " << 5*boxLength[2] << ";\n"
        "#declare RAD = off;\n"
        "global_settings {\n"
        "#if(RAD)\n"
        "radiosity {\n"
        "pretrace_start 0.08\n"
        "pretrace_end   0.01\n"
        "count 500\n"
        "nearest_count 10\n"
        "error_bound 0.02\n"
        "recursion_limit 1\n"
        "low_error_factor 0.2\n"
        "gray_threshold 0.0\n"
        "minimum_reuse 0.015\n"
        "brightness 1.4\n"
        "adc_bailout 0.01/2\n"
        "}\n"
        "#end\n"
        "}\n"
        "\n"
        "camera{\n"
        "orthographic\n"
        "location < 0, 0, zoom >\n"
        "direction < 0, 0, 2 >\n"
        "up < 0, " << boxLength[1] << ", 0 >\n"
        "// Ratio is negative to switch povray to a right hand coordinate system.\n"
        "right < -" << boxLength[0] << ", 0, 0 >\n"
        "look_at < 0, 0, 0 >\n"
//        "    location < 0,0, 0.5333333*" << boxLength << "*zoom >\n"
//        //"    location < 0,0, 0.8*" << boxLength << "*zoom >\n"
//        "        direction < 0,0, zoom >\n"
//        "            // Ratio is negative to switch povray to a right hand coordinate system.\n"
//        "                right < -Ratio ,0 , 0 >\n"
//        "                    look_at < 0, -3.5, 0 >\n"
        "}\n"
        "\n"
        "background { color rgb 1 }\n"
        "\n"
        "global_settings { ambient_light rgb 1 }\n"
        "light_source { < 0, 0, 10*zoom > rgb < 1.0, 1.0, 1.0 > }\n"
        "light_source { < 0, 0, 10*zoom > rgb < 0.1, 0.1, 0.1 > }\n\n";

    return;
}

void PovObjects::printHeader2(ofstream& outputFile){
    outputFile <<
        "  //**************************************\n"
        "  // generated by nucleation_tracker\n"
        "  //**************************************\n"
        "\n"
        "  //**************************************\n"
        "  // Lights, camera, resolution!\n"
        "  //**************************************\n"
        "  global_settings{ max_trace_level 100 }\n"
        "\n"
        "#declare Ratio = 1;\n"
        "#declare zoom = 10;\n"
        "#declare RAD = off;\n"
        "global_settings {\n"
        "#if(RAD)\n"
        "radiosity {\n"
        "pretrace_start 0.08\n"
        "pretrace_end   0.01\n"
        "count 500\n"
        "nearest_count 10\n"
        "error_bound 0.02\n"
        "recursion_limit 1\n"
        "low_error_factor 0.2\n"
        "gray_threshold 0.0\n"
        "minimum_reuse 0.015\n"
        "brightness 1.4\n"
        "adc_bailout 0.01/2\n"
        "}\n"
        "#end\n"
        "}\n"
        "\n"
        "camera{\n"
        "orthographic\n"
        "location < 0, 0, zoom >\n"
        "direction < 0, 0, 2 >\n"
        "up < 0, 6, 0 >\n"
        "// Ratio is negative to switch povray to a right hand coordinate system.\n"
        "right < -6, 0, 0 >\n"
        "look_at < 0, 0, 0 >\n"
        "}\n"
        "\n"
        "background { color rgb 1 }\n"
        "\n"
        "global_settings { ambient_light rgb 1 }\n"
        "light_source { < 0, 0, 10*zoom > rgb < 1.0, 1.0, 1.0 > }\n"
        "light_source { < 0, 0, 10*zoom > rgb < 0.1, 0.1, 0.1 > }\n\n";

    return;
}

void PovObjects::printRing3(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring3 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring3_dot_size << "\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < " << -boxLength[0] << ", " << -boxLength[1] << ", center_z-0.01 >,< "
        << boxLength[0] << ", " << boxLength[0] << ", center_z+0.01 >\n"
        "        texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "  }\n"
        "#end\n";
    return;
}

void PovObjects::printRing4(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring4 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring4_dot_size << "\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < " << -boxLength[0] << ", " << -boxLength[1] << ", center_z-0.01 >,< "
        << boxLength[0] << ", " << boxLength[0] << ", center_z+0.01 >\n"
        "        texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "  }\n"
        "#end\n";
    return;
}

void PovObjects::printRing5(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring5 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring5_dot_size << "\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < " << -boxLength[0] << ", " << -boxLength[1] << ", center_z-0.01 >,< "
        << boxLength[0] << ", " << boxLength[0] << ", center_z+0.01 >\n"
        "        texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "  }\n"
        "#end\n";
    return;
}

void PovObjects::printRing6(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring6 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring6_dot_size << "\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < " << -boxLength[0] << ", " << -boxLength[1] << ", center_z-0.01 >,< "
        << boxLength[0] << ", " << boxLength[0] << ", center_z+0.01 >\n"
        "        texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "  }\n"
        "#end\n";
    return;
}

void PovObjects::printRing7(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring7 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring7_dot_size << "\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < " << -boxLength[0] << ", " << -boxLength[1] << ", center_z-0.01 >,< "
        << boxLength[0] << ", " << boxLength[0] << ", center_z+0.01 >\n"
        "        texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "  }\n"
        "#end\n";
    return;
}

void PovObjects::printBar(ofstream& outputFile){
    outputFile << "#macro bar (bottom_x, high, wide, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "    box{\n"
        "      < bottom_x+(1-wide)/2, -2, 0 >,< bottom_x+(1-wide)/2+wide, -2+high, 0.1>\n"
        "       texture{\n"
        "          pigment{ rgbt < outerRed, outerGreen, outerBlue, dot_transparency > }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < bottom_x+(1-wide)/2-" << outline << ", -2-" << outline << ", 0 >,< bottom_x+(1-wide)/2+wide+" << outline << ", -2+high+" << outline << ", 0 >\n"
        "        texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "#end\n";
    return;
}

void PovObjects::printAxes(ofstream& outputFile){
    outputFile << "#macro axes (wide)\n " 
        "    box{\n"
        "      < -2-wide, -2-wide, 0.2 >,< -2, 2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2-wide, -2-wide, 0.2 >,< 2, -2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2-2*wide, 2-wide, 0.2 >,< -2, 2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2-2*wide, 2-wide, -0.2 >,< 2, 2, -0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0.65 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2-2*wide, -wide, 0.2 >,< -2, 0, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2-1.5*wide, -wide, -0.2 >,< 2, 0, -0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0.65 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -2+0.5*(1-wide), -2-2*wide, 0.2 >,< -2+0.5*(1+wide), -2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < -1+0.5*(1-wide), -2-2*wide, 0.2 >,< -1+0.5*(1+wide), -2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < 0.5*(1-wide), -2-2*wide, 0.2 >,< 0.5*(1+wide), -2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "    box{\n"
        "      < 1+0.5*(1-wide), -2-2*wide, 0.2 >,< 1+0.5*(1+wide), -2, 0.2>\n"
        "       texture{\n"
        "          pigment{ rgbt 0 }\n"
        "          finish{\n"
        "            ambient .2\n"
        "              diffuse .6\n"
        "              specular .1\n"
        "              roughness .1\n"
        "          }\n"
        "        }\n"
        "    }\n"
        "#end\n";
    return;
}


int main(int argc, char *argv[]) {
    bool trigger;
    bool isPos;
    char povDistFileName[200], povFileName[200], trajFileName[200], fileName[200], inLine[1000], inValue[200], tag[100];
    char *token;
    const char *file;
    const char *delimit = " \t\n";
    const char *period = ".";
    const char *OWAtom = "   OW";
    const char *HWAtom1 = "  HW1";
    const char *HWAtom2 = "  HW2";
    const char *OAtom = "    O";
    const char *HAtom1 = "   H1";
    const char *HAtom2 = "   H2";
    int count = 0;
    int i, j, k, l, m, n, o, p, q, nAtoms, waterModel, waterNumber;
    int oCount, hCount, mCount, totalWaterCount, hbondCount, neighborCount, frameCount;
    int index2, index3, index4, index5, index6, index7, index8, index9;
    int lastRing3, lastRing4, lastRing5, lastRing6, lastRing7, lastRing8;
    int x_frame, y_frame;
    int totalRings;
    double xVal, yVal, zVal;
    double temp_dist2;
    double ionPosition[3], position[3];
    double water1[9], water2[9];
    double boxLength[3];
    double tempPosX[8], tempPosY[8], tempPosZ[8];
    double diffVal;
    double ringProbabilities[9];
    vector<int> atomCount;
    vector<int> neighborListIndex;
    vector<int> hbondListIndex;
    vector<int> tempVec;
    vector<double> oPosX, oPosY, oPosZ;
    vector<double> hPosX, hPosY, hPosZ;
    vector<vector<int> > neighborList;
    vector<vector<int> > hbondList;
    vector<vector<int> > ring3members;
    vector<vector<int> > ring4members;
    vector<vector<int> > ring5members;
    vector<vector<int> > ring6members;
    vector<vector<int> > ring7members;
    vector<vector<int> > ring8members;
    string strungName;
    string povName;
    string povName2;
    string povDistName;
    string red_val, green_val, blue_val;
    string lineString;
    string shortString;
    string frameCountString, xFrameString, yFrameString;
    stringstream frameInt, xFrameInt, yFrameInt;
    ofstream pov_out;
    ofstream povDistOut;
    ofstream ringTrajOut;

    if (argc != 2) {
        cout << "\nUsage: " << argv[0] << " [file name].gro\n\n";
        return 0;
    }

    // Now try opening the file
    ifstream prayer(argv[1]);

    // Make sure the file exists
    if (!prayer) { 
        cout << "Unable to open " << argv[1] << " for reading.\n";
        return 0;
    }
    prayer.close();

    // create directory for pov files
    const int dir_err = system("mkdir -p pov_files");
    if (-1 == dir_err)
    {
        printf("Error creating directory!n");
        exit(1);
    }


    // Build a filename string from the .status file name
    strungName = argv[1];
    file = strungName.c_str();
    strcpy(fileName, file);

    token = strtok(fileName, period);
    strcpy(fileName, token);
    strcpy(trajFileName, token);
    strcpy(povFileName, token);
    strcpy(povDistFileName, token);
    povName2 = povFileName;
    povDistName = povDistFileName;

    strcat(fileName,"_nuc_info.txt");
    strcat(trajFileName,"_ringtrj.xyz");
    strcat(povFileName,"_pov.txt");
    strungName = fileName;

    // Build readers and writers
    ifstream inputer(argv[1]);
    ofstream outputer(fileName);
    ofstream outputer_xyz(trajFileName);
    ofstream outputer_pov(povFileName);
    ofstream outputer_dist_pov(povDistName);

    // initialize our calculator
    Calcs *calculator = new Calcs(); 

    // output ring info. in pov files 
    PovObjects *povrayObjects = new PovObjects();

    // Read the .gro file and load the positions
    cout << "\nLoading and processing trajectory...\n";
    outputer << setw(8) << "# frame" <<  setw(8) << "1" << setw(8) << "3" << setw(8) << "4" << setw(8) << "5" << setw(8) << "6" << setw(8) << "7" << "\n"; //setw(8) << "8" << "\n";

    frameCount = 1;
    totalWaterCount = 0;
    hbondCount = 0;
    neighborCount = 0;
    lastRing3 = 0;
    lastRing4 = 0;
    lastRing5 = 0;
    lastRing6 = 0;
    lastRing7 = 0;
    lastRing8 = 0;
    totalRings = 0;

    inputer.getline(inLine,999,'\n');
    inputer.getline(inLine,999,'\n');
    token = strtok(inLine,delimit);
    strcpy(inValue,token);
    nAtoms = atoi(inValue);
    atomCount.push_back(nAtoms);
    while (!inputer.eof()) {
        // grab the water atom positions
        for(i=0; i<nAtoms; i++){
            inputer.getline(inLine,999,'\n');
            lineString = inLine;
            shortString = lineString.substr(10,5);
            //token = strtok(inLine,delimit);
            //token = strtok(NULL,delimit);
            strcpy(token, shortString.c_str());
            // test if it is an oxygen atom and load in the appropriate vector
            if (!strcmp(OAtom, token) || !strcmp(OWAtom, token)){
                //token = strtok(NULL,delimit);
                //token = strtok(NULL,delimit);
                shortString = lineString.substr(20,8);
                strcpy(inValue,shortString.c_str());
                //strcpy(inValue,token);
                xVal = atof(inValue) * 10.0;
                oPosX.push_back(xVal);
                //token = strtok(NULL,delimit);
                shortString = lineString.substr(28,8);
                strcpy(inValue,shortString.c_str());
                //strcpy(inValue,token);
                yVal = atof(inValue) * 10.0;
                oPosY.push_back(yVal);
                shortString = lineString.substr(36,8);
                strcpy(inValue,shortString.c_str());
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                zVal = atof(inValue) * 10.0;
                oPosZ.push_back(zVal);
                // test if it is a hydrogen atom and load in the appropriate vector
            } else if (!strcmp(HAtom1, token) || !strcmp(HWAtom1, token)){
                //token = strtok(NULL,delimit);
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                shortString = lineString.substr(20,8);
                strcpy(inValue,shortString.c_str());
                xVal = atof(inValue) * 10.0;
                hPosX.push_back(xVal);
                shortString = lineString.substr(28,8);
                strcpy(inValue,shortString.c_str());
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                yVal = atof(inValue) * 10.0;
                hPosY.push_back(yVal);
                shortString = lineString.substr(36,8);
                strcpy(inValue,shortString.c_str());
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                zVal = atof(inValue) * 10.0;
                hPosZ.push_back(zVal);
                // test if it is a hydrogen atom and load in the appropriate vector
            } else if (!strcmp(HAtom2, token) || !strcmp(HWAtom2, token)){
                //token = strtok(NULL,delimit);
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                shortString = lineString.substr(20,8);
                strcpy(inValue,shortString.c_str());
                xVal = atof(inValue) * 10.0;
                hPosX.push_back(xVal);
                shortString = lineString.substr(28,8);
                strcpy(inValue,shortString.c_str());
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                yVal = atof(inValue) * 10.0;
                hPosY.push_back(yVal);
                shortString = lineString.substr(36,8);
                strcpy(inValue,shortString.c_str());
                //token = strtok(NULL,delimit);
                //strcpy(inValue,token);
                zVal = atof(inValue) * 10.0;
                hPosZ.push_back(zVal);
                totalWaterCount++;
            }
        }

        // the box information - used in minimum image wrapping 
        inputer.getline(inLine,999,'\n');
        token = strtok(inLine,delimit);
        strcpy(inValue,token);
        boxLength[0] = atof(inValue) * 10.0;	
        token = strtok(NULL,delimit);
        strcpy(inValue,token);
        boxLength[1] = atof(inValue) * 10.0;	
        token = strtok(NULL,delimit);
        strcpy(inValue,token);
        boxLength[2] = atof(inValue) * 10.0;	

        // okay, let's do some calculations on this frame

        // okay, now we loop over the waters in the frame to determine if they
        // are neighbors
        for (i=0; i<oPosX.size(); i++){
            // can have up to 10 neighbors
            // vec2D.push_back(std::vector<int>(4, 11));
            neighborListIndex.push_back(0);
            neighborList.push_back(vector<int>(10,-1));
            // can have up to 10 HBonds
            hbondListIndex.push_back(0);
            hbondList.push_back(vector<int>(10,-1));
        }

        for (i=0; i<oPosX.size()-1; i++){
            // just use o and h vectors
            // load the reference water
            water1[0] = oPosX[i];
            water1[1] = oPosY[i];
            water1[2] = oPosZ[i];
            water1[3] = hPosX[2*i];
            water1[4] = hPosY[2*i];
            water1[5] = hPosZ[2*i];
            water1[6] = hPosX[2*i+1];
            water1[7] = hPosY[2*i+1];
            water1[8] = hPosZ[2*i+1];

            for (j=i+1; j<oPosX.size(); j++){
                // load the water
                water2[0] = oPosX[j];
                water2[1] = oPosY[j];
                water2[2] = oPosZ[j];
                water2[3] = hPosX[2*j];
                water2[4] = hPosY[2*j];
                water2[5] = hPosZ[2*j];
                water2[6] = hPosX[2*j+1];
                water2[7] = hPosY[2*j+1];
                water2[8] = hPosZ[2*j+1];

                // calculate if they are neighbors
                xVal = water2[0]-water1[0];
                yVal = water2[1]-water1[1];
                zVal = water2[2]-water1[2];
                // do minimum image filter...
                xVal -= boxLength[0]*copysign(1.0,xVal)*floor(fabs(xVal/boxLength[0])+0.5);
                yVal -= boxLength[1]*copysign(1.0,yVal)*floor(fabs(yVal/boxLength[1])+0.5);
                zVal -= boxLength[2]*copysign(1.0,zVal)*floor(fabs(zVal/boxLength[2])+0.5);

                temp_dist2 = xVal*xVal + yVal*yVal + zVal*zVal;
                if(temp_dist2 <= distance_tol2){
                    neighborList[i][neighborListIndex[i]] = j;
                    neighborListIndex[i]++;
                    neighborList[j][neighborListIndex[j]] = i;
                    neighborListIndex[j]++;
                    neighborCount++;
                }
            }
        }

        // now build our hbond lists
        for (i=0; i<neighborList.size(); i++){
            // load the reference water
            water1[0] = oPosX[i];
            water1[1] = oPosY[i];
            water1[2] = oPosZ[i];
            water1[3] = hPosX[2*i];
            water1[4] = hPosY[2*i];
            water1[5] = hPosZ[2*i];
            water1[6] = hPosX[2*i+1];
            water1[7] = hPosY[2*i+1];
            water1[8] = hPosZ[2*i+1];

            for (j=0; j<=neighborListIndex[i]; j++){
                if (neighborList[i][j] > i){
                    // load the water
                    water2[0] = oPosX[neighborList[i][j]];
                    water2[1] = oPosY[neighborList[i][j]];
                    water2[2] = oPosZ[neighborList[i][j]];
                    water2[3] = hPosX[2*neighborList[i][j]];
                    water2[4] = hPosY[2*neighborList[i][j]];
                    water2[5] = hPosZ[2*neighborList[i][j]];
                    water2[6] = hPosX[2*neighborList[i][j]+1];
                    water2[7] = hPosY[2*neighborList[i][j]+1];
                    water2[8] = hPosZ[2*neighborList[i][j]+1];

                    if(calculator->IsWat3HBond(water1, water2, boxLength)){
                        hbondList[i][hbondListIndex[i]] = neighborList[i][j];
                        hbondListIndex[i]++;
                        hbondList[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = i;
                        hbondListIndex[neighborList[i][j]]++;
                        hbondCount++;
                    }
                }
            }
        }

        // now we start identifying rings to build those lists
        for(i=0; i<hbondListIndex.size(); i++){
            for(j=0; j<hbondListIndex[i]; j++){
                index2 = hbondList[i][j]; 
                if (index2 != i){
                    for(k=0; k<hbondListIndex[index2]; k++){
                        index3 = hbondList[index2][k];
                        if (index3 != index2 && index3 != i){
                            for(l=0; l<hbondListIndex[index3]; l++){
                                index4 = hbondList[index3][l];
                                // test if a 3 member ring
                                if (index4 == i){
                                    tempVec.push_back(i); 
                                    tempVec.push_back(index2); 
                                    tempVec.push_back(index3);
                                    sort(tempVec.begin(),tempVec.end());
                                    ring3members.push_back(tempVec);
                                    tempVec.clear();
                                    break;
                                }
                                if (index4 != index3 && index4 != index2){
                                    for(m=0; m<hbondListIndex[index4]; m++){
                                        index5 = hbondList[index4][m];
                                        // test if a 4 member ring
                                        if (index5 == i){
                                            tempVec.push_back(i); 
                                            tempVec.push_back(index2); 
                                            tempVec.push_back(index3);
                                            tempVec.push_back(index4);
                                            sort(tempVec.begin(),tempVec.end());
                                            ring4members.push_back(tempVec);
                                            tempVec.clear();
                                            break;
                                        }
                                        if (index5 != index4 && index5 != index3 && index5 != index2){
                                            for(n=0; n<hbondListIndex[index5]; n++){
                                                index6 = hbondList[index5][n];
                                                // test if a 5 member ring
                                                if (index6 == i){
                                                    tempVec.push_back(i); 
                                                    tempVec.push_back(index2); 
                                                    tempVec.push_back(index3);
                                                    tempVec.push_back(index4);
                                                    tempVec.push_back(index5);
                                                    sort(tempVec.begin(),tempVec.end());
                                                    ring5members.push_back(tempVec);
                                                    tempVec.clear();
                                                    break;
                                                }
                                                if (index6 != index5 && index6 != index4 && index6 != index3 && index6 != index2){
                                                    for(o=0; o<hbondListIndex[index6]; o++){
                                                        index7 = hbondList[index6][o];
                                                        // test if a 5 member ring
                                                        if (index7 == i){
                                                            tempVec.push_back(i); 
                                                            tempVec.push_back(index2); 
                                                            tempVec.push_back(index3);
                                                            tempVec.push_back(index4);
                                                            tempVec.push_back(index5);
                                                            tempVec.push_back(index6);
                                                            sort(tempVec.begin(),tempVec.end());
                                                            ring6members.push_back(tempVec);
                                                            tempVec.clear();
                                                            break;
                                                        }
                                                        if (index7 != index6 && index7 != index5 && index7 != index4 && index7 != index3 && index7 != index2){
                                                            for(p=0; p<hbondListIndex[index7]; p++){
                                                                index8 = hbondList[index7][p];
                                                                // test if a 5 member ring
                                                                if (index8 == i){
                                                                    tempVec.push_back(i); 
                                                                    tempVec.push_back(index2); 
                                                                    tempVec.push_back(index3);
                                                                    tempVec.push_back(index4);
                                                                    tempVec.push_back(index5);
                                                                    tempVec.push_back(index6);
                                                                    tempVec.push_back(index7);
                                                                    sort(tempVec.begin(),tempVec.end());
                                                                    ring7members.push_back(tempVec);
                                                                    tempVec.clear();
                                                                    break;
                                                                }
                                                                /*
                                                                   if (index8 != index7 && index8 != index6 && index8 != index5 && index8 != index4 && index8 != index3 && index8 != index2){
                                                                   for(q=0; q<hbondListIndex[index8]; q++){
                                                                   index9 = hbondList[index8][q];
                                                                // test if a 5 member ring
                                                                if (index9 == i){
                                                                tempVec.push_back(i); 
                                                                tempVec.push_back(index2); 
                                                                tempVec.push_back(index3);
                                                                tempVec.push_back(index4);
                                                                tempVec.push_back(index5);
                                                                tempVec.push_back(index6);
                                                                tempVec.push_back(index7);
                                                                tempVec.push_back(index8);
                                                                sort(tempVec.begin(),tempVec.end());
                                                                ring8members.push_back(tempVec);
                                                                tempVec.clear();
                                                                break;
                                                                }
                                                                }
                                                                }
                                                                */
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // purge atom i from hbond lists 
            for(j=0; j<hbondListIndex[i]; j++){
                index2 = hbondList[i][j];
                for(k=0; k<hbondListIndex[index2]; k++){
                    hbondList[index2][k] = hbondList[index2][k+1];
                }
                hbondListIndex[index2]--;
            }
        }

        // prune the subset information
        // ...self ring lists first...
        for(i=ring3members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring3members[j].begin(),ring3members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring3members.erase(ring3members.begin() + i);
                    break;
                }
            }
        }
        for(i=ring4members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring4members[j].begin(),ring4members[j].end(), ring4members[i].begin(), ring4members[i].end())){
                    ring4members.erase(ring4members.begin() + i);
                    break;
                }
            }
        }
        for(i=ring5members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring5members[j].begin(),ring5members[j].end(), ring5members[i].begin(), ring5members[i].end())){
                    ring5members.erase(ring5members.begin() + i);
                    break;
                }
            }
        }
        for(i=ring6members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring6members[j].begin(),ring6members[j].end(), ring6members[i].begin(), ring6members[i].end())){
                    ring6members.erase(ring6members.begin() + i);
                    break;
                }
            }
        }
        for(i=ring7members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring7members[j].begin(),ring7members[j].end(), ring7members[i].begin(), ring7members[i].end())){
                    ring7members.erase(ring7members.begin() + i);
                    break;
                }
            }
        }
        /*
           for(i=ring8members.size()-1; i>0; i--){
           for (j=i-1; j>=0; j--){
           if(includes(ring8members[j].begin(),ring8members[j].end(), ring8members[i].begin(), ring8members[i].end())){
           ring8members.erase(ring8members.begin() + i);
           break;
           }
           }
           }
           */
        // now we do the pruning tree...
        // ...3s in 4s...
        for(i=ring3members.size()-1; i>=0; i--){
            for (j=ring4members.size()-1; j>=0; j--){
                if(includes(ring4members[j].begin(),ring4members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring4members.erase(ring4members.begin() + j);
                    break;
                }
            }
        }
        // ...3s in 5s...
        for(i=ring3members.size()-1; i>=0; i--){
            for (j=ring5members.size()-1; j>=0; j--){
                if(includes(ring5members[j].begin(),ring5members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring5members.erase(ring5members.begin() + j);
                    break;
                }
            }
        }
        // ...4s in 5s...
        for(i=ring4members.size()-1; i>=0; i--){
            for (j=ring5members.size()-1; j>=0; j--){
                if(includes(ring5members[j].begin(),ring5members[j].end(), ring4members[i].begin(), ring4members[i].end())){
                    ring5members.erase(ring5members.begin() + j);
                    break;
                }
            }
        }
        // ...3s in 6s...
        for(i=ring3members.size()-1; i>=0; i--){
            for (j=ring6members.size()-1; j>=0; j--){
                if(includes(ring6members[j].begin(),ring6members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring6members.erase(ring6members.begin() + j);
                    break;
                }
            }
        }
        // ...4s in 6s...
        for(i=ring4members.size()-1; i>=0; i--){
            for (j=ring6members.size()-1; j>=0; j--){
                if(includes(ring6members[j].begin(),ring6members[j].end(), ring4members[i].begin(), ring4members[i].end())){
                    ring6members.erase(ring6members.begin() + j);
                    break;
                }
            }
        }
        // ...5s in 6s...
        for(i=ring5members.size()-1; i>=0; i--){
            for (j=ring6members.size()-1; j>=0; j--){
                if(includes(ring6members[j].begin(),ring6members[j].end(), ring5members[i].begin(), ring5members[i].end())){
                    ring6members.erase(ring6members.begin() + j);
                    break;
                }
            }
        }
        // ...3s in 7s...
        for(i=ring3members.size()-1; i>=0; i--){
            for (j=ring7members.size()-1; j>=0; j--){
                if(includes(ring7members[j].begin(),ring7members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring7members.erase(ring7members.begin() + j);
                    break;
                }
            }
        }
        // ...4s in 7s...
        for(i=ring4members.size()-1; i>=0; i--){
            for (j=ring7members.size()-1; j>=0; j--){
                if(includes(ring7members[j].begin(),ring7members[j].end(), ring4members[i].begin(), ring4members[i].end())){
                    ring7members.erase(ring7members.begin() + j);
                    break;
                }
            }
        }
        // ...5s in 7s...
        for(i=ring5members.size()-1; i>=0; i--){
            for (j=ring7members.size()-1; j>=0; j--){
                if(includes(ring7members[j].begin(),ring7members[j].end(), ring5members[i].begin(), ring5members[i].end())){
                    ring7members.erase(ring7members.begin() + j);
                    break;
                }
            }
        }
        // ...6s in 7s...
        for(i=ring6members.size()-1; i>=0; i--){
            for (j=ring7members.size()-1; j>=0; j--){
                if(includes(ring7members[j].begin(),ring7members[j].end(), ring6members[i].begin(), ring6members[i].end())){
                    ring7members.erase(ring7members.begin() + j);
                    break;
                }
            }
        }
        /*
        // ...3s in 8s...
        for(i=ring3members.size()-1; i>=0; i--){
        for (j=ring8members.size()-1; j>=0; j--){
        if(includes(ring8members[j].begin(),ring8members[j].end(), ring3members[i].begin(), ring3members[i].end())){
        ring8members.erase(ring8members.begin() + j);
        break;
        }
        }
        }
        // ...4s in 8s...
        for(i=ring4members.size()-1; i>=0; i--){
        for (j=ring8members.size()-1; j>=0; j--){
        if(includes(ring8members[j].begin(),ring8members[j].end(), ring4members[i].begin(), ring4members[i].end())){
        ring8members.erase(ring8members.begin() + j);
        break;
        }
        }
        }
        // ...5s in 8s...
        for(i=ring5members.size()-1; i>=0; i--){
        for (j=ring8members.size()-1; j>=0; j--){
        if(includes(ring8members[j].begin(),ring8members[j].end(), ring5members[i].begin(), ring5members[i].end())){
        ring8members.erase(ring8members.begin() + j);
        break;
        }
        }
        }
        // ...6s in 8s...
        for(i=ring6members.size()-1; i>=0; i--){
        for (j=ring8members.size()-1; j>=0; j--){
        if(includes(ring8members[j].begin(),ring8members[j].end(), ring6members[i].begin(), ring6members[i].end())){
        ring8members.erase(ring8members.begin() + j);
        break;
        }
        }
        }
        // ...7s in 8s...
        for(i=ring7members.size()-1; i>=0; i--){
        for (j=ring8members.size()-1; j>=0; j--){
        if(includes(ring8members[j].begin(),ring8members[j].end(), ring7members[i].begin(), ring7members[i].end())){
        ring8members.erase(ring8members.begin() + j);
        break;
        }
        }
        }
        */

        // output ring data
        lastRing3 = ring3members.size()-lastRing3;
        lastRing4 = ring4members.size()-lastRing4;
        lastRing5 = ring5members.size()-lastRing5;
        lastRing6 = ring6members.size()-lastRing6;
        lastRing7 = ring7members.size()-lastRing7;
        lastRing8 = ring8members.size()-lastRing8;
        outputer << setw(8) << frameCount ;
        outputer << setw(8) << hbondCount ;
        outputer << setw(8) << lastRing3 ;
        outputer << setw(8) << lastRing4 ;
        outputer << setw(8) << lastRing5 ;
        outputer << setw(8) << lastRing6 ;
        outputer << setw(8) << lastRing7 << "\n";
        //totalRings = lastRing3+lastRing4+lastRing5+lastRing6+lastRing7;
        totalRings = lastRing4+lastRing5+lastRing6+lastRing7;
        ringProbabilities[0] = 0;
        ringProbabilities[1] = 0;
        ringProbabilities[2] = 0;
        ringProbabilities[3] = (float)lastRing3/(totalRings+lastRing3);
        ringProbabilities[4] = (float)lastRing4/totalRings;
        ringProbabilities[5] = (float)lastRing5/totalRings;
        ringProbabilities[6] = (float)lastRing6/totalRings;
        ringProbabilities[7] = (float)lastRing7/totalRings;
        ringProbabilities[8] = (float)lastRing8/totalRings;
        outputer_xyz << totalRings << "\n";
        outputer_xyz <<fixed<<setprecision(5)<<setw(8) << boxLength[0] <<" "<< boxLength[1] <<" "<< boxLength[2] << "\n";

        // outputer << setw(8) << lastRing7 ;
        // outputer << setw(8) << lastRing8 << "\n";

        // open a pov outputer
        frameInt.str("");
        frameInt << frameCount;
        frameCountString = frameInt.str();
        povName = "pov_files/" + povName2 + "_" + frameCountString + ".pov";
        povDistName = "pov_files/" + povName2 + "_dist_" + frameCountString + ".pov";
        pov_out.open(povName.c_str());
        povDistOut.open(povDistName.c_str());

        // print out the header and object info. to the .pov file
        povrayObjects->printHeader(pov_out, boxLength);
        povrayObjects->printRing3(pov_out, boxLength);
        povrayObjects->printRing4(pov_out, boxLength);
        povrayObjects->printRing5(pov_out, boxLength);
        povrayObjects->printRing6(pov_out, boxLength);
        povrayObjects->printRing7(pov_out, boxLength);

        // print out the header and object info. to the distribution .pov file
        povrayObjects->printHeader2(povDistOut);
        povrayObjects->printBar(povDistOut);
        povrayObjects->printAxes(povDistOut);
        povDistOut << "axes(" << axes_width << ")\n";
        // now print the distribution bars, since we have them!
        //povDistOut << "bar(-2.5, " <<setprecision(8)<<setw(10)<< 5.0*ringProbabilities[3] << ", 0.7, 0.2, 0.2, 1, 0.2)\n";
        povDistOut << "bar(-2.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[4] << ", 0.7, 0.5, 0.0, 1, 0.0)\n";
        povDistOut << "bar(-1.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[5] << ", 0.7, 0.5, 0.5, 1, 0.0)\n";
        povDistOut << "bar(0.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[6] << ", 0.7, 1.0, 0.5, 0.5, 0.0)\n";
        povDistOut << "bar(1.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[7] << ", 0.7, 1.0, 0.75, 0.5, 0.0)\n";

        for (i=0; i<ring3members.size(); i++){
            tempPosX[0]=oPosX[ring3members[i][0]];
            diffVal = oPosX[ring3members[i][1]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring3members[i][1]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring3members[i][1]]+boxLength[0];
            } else {
                tempPosX[1] = oPosX[ring3members[i][1]];
            }
            diffVal = oPosX[ring3members[i][2]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring3members[i][2]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring3members[i][2]]+boxLength[0];
            } else {
                tempPosX[2] = oPosX[ring3members[i][2]];
            }
            xVal = tempPosX[0]+tempPosX[1]+tempPosX[2];
            xVal /= 3;

            tempPosY[0]=oPosY[ring3members[i][0]];
            diffVal = oPosY[ring3members[i][1]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring3members[i][1]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring3members[i][1]]+boxLength[1];
            } else {
                tempPosY[1] = oPosY[ring3members[i][1]];
            }
            diffVal = oPosY[ring3members[i][2]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring3members[i][2]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring3members[i][2]]+boxLength[1];
            } else {
                tempPosY[2] = oPosY[ring3members[i][2]];
            }
            yVal = tempPosY[0]+tempPosY[1]+tempPosY[2];
            yVal /= 3;

            tempPosZ[0]=oPosZ[ring3members[i][0]];
            diffVal = oPosZ[ring3members[i][1]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring3members[i][1]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring3members[i][1]]+boxLength[2];
            } else {
                tempPosZ[1] = oPosZ[ring3members[i][1]];
            }
            diffVal = oPosZ[ring3members[i][2]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring3members[i][2]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring3members[i][2]]+boxLength[2];
            } else {
                tempPosZ[2] = oPosZ[ring3members[i][2]];
            }
            zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2];
            zVal /= 3;
        }
        //outputer_xyz << setw(5) << "Li " << setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
        outputer_xyz <<left<< setw(5) << "Li" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";

        for (i=0; i<ring4members.size(); i++){
            tempPosX[0]=oPosX[ring4members[i][0]];
            diffVal = oPosX[ring4members[i][1]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring4members[i][1]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring4members[i][1]]+boxLength[0];
            } else {
                tempPosX[1] = oPosX[ring4members[i][1]];
            }
            diffVal = oPosX[ring4members[i][2]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring4members[i][2]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring4members[i][2]]+boxLength[0];
            } else {
                tempPosX[2] = oPosX[ring4members[i][2]];
            }
            diffVal = oPosX[ring4members[i][3]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring4members[i][3]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring4members[i][3]]+boxLength[0];
            } else {
                tempPosX[3] = oPosX[ring4members[i][3]];
            }
            xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3];
            xVal /= 4;

            tempPosY[0]=oPosY[ring4members[i][0]];
            diffVal = oPosY[ring4members[i][1]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring4members[i][1]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring4members[i][1]]+boxLength[1];
            } else {
                tempPosY[1] = oPosY[ring4members[i][1]];
            }
            diffVal = oPosY[ring4members[i][2]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring4members[i][2]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring4members[i][2]]+boxLength[1];
            } else {
                tempPosY[2] = oPosY[ring4members[i][2]];
            }
            diffVal = oPosY[ring4members[i][3]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring4members[i][3]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring4members[i][3]]+boxLength[1];
            } else {
                tempPosY[3] = oPosY[ring4members[i][3]];
            }
            yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3];
            yVal /= 4;

            tempPosZ[0]=oPosZ[ring4members[i][0]];
            diffVal = oPosZ[ring4members[i][1]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring4members[i][1]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring4members[i][1]]+boxLength[2];
            } else {
                tempPosZ[1] = oPosZ[ring4members[i][1]];
            }
            diffVal = oPosZ[ring4members[i][2]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring4members[i][2]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring4members[i][2]]+boxLength[2];
            } else {
                tempPosZ[2] = oPosZ[ring4members[i][2]];
            }
            diffVal = oPosZ[ring4members[i][3]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring4members[i][3]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring4members[i][3]]+boxLength[2];
            } else {
                tempPosZ[3] = oPosZ[ring4members[i][3]];
            }
            zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3];
            zVal /= 4;
        }
        //outputer_xyz << setw(5) << "Be" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
        outputer_xyz <<left<< setw(5) << "Be" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";

        red_val = "0.0";
        green_val = "0.3";
        blue_val = "1.0";
        for (i=0; i<ring5members.size(); i++){
            tempPosX[0]=oPosX[ring5members[i][0]];
            diffVal = oPosX[ring5members[i][1]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring5members[i][1]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring5members[i][1]]+boxLength[0];
            } else {
                tempPosX[1] = oPosX[ring5members[i][1]];
            }
            diffVal = oPosX[ring5members[i][2]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring5members[i][2]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring5members[i][2]]+boxLength[0];
            } else {
                tempPosX[2] = oPosX[ring5members[i][2]];
            }
            diffVal = oPosX[ring5members[i][3]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring5members[i][3]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring5members[i][3]]+boxLength[0];
            } else {
                tempPosX[3] = oPosX[ring5members[i][3]];
            }
            diffVal = oPosX[ring5members[i][4]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring5members[i][4]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring5members[i][4]]+boxLength[0];
            } else {
                tempPosX[4] = oPosX[ring5members[i][4]];
            }
            xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4];
            xVal /= 5;

            tempPosY[0]=oPosY[ring5members[i][0]];
            diffVal = oPosY[ring5members[i][1]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring5members[i][1]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring5members[i][1]]+boxLength[1];
            } else {
                tempPosY[1] = oPosY[ring5members[i][1]];
            }
            diffVal = oPosY[ring5members[i][2]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring5members[i][2]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring5members[i][2]]+boxLength[1];
            } else {
                tempPosY[2] = oPosY[ring5members[i][2]];
            }
            diffVal = oPosY[ring5members[i][3]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring5members[i][3]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring5members[i][3]]+boxLength[1];
            } else {
                tempPosY[3] = oPosY[ring5members[i][3]];
            }
            diffVal = oPosY[ring5members[i][4]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring5members[i][4]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring5members[i][4]]+boxLength[1];
            } else {
                tempPosY[4] = oPosY[ring5members[i][4]];
            }
            yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4];
            yVal /= 5;

            tempPosZ[0]=oPosZ[ring5members[i][0]];
            diffVal = oPosZ[ring5members[i][1]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring5members[i][1]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring5members[i][1]]+boxLength[2];
            } else {
                tempPosZ[1] = oPosZ[ring5members[i][1]];
            }
            diffVal = oPosZ[ring5members[i][2]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring5members[i][2]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring5members[i][2]]+boxLength[2];
            } else {
                tempPosZ[2] = oPosZ[ring5members[i][2]];
            }
            diffVal = oPosZ[ring5members[i][3]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring5members[i][3]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring5members[i][3]]+boxLength[2];
            } else {
                tempPosZ[3] = oPosZ[ring5members[i][3]];
            }
            diffVal = oPosZ[ring5members[i][4]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring5members[i][4]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring5members[i][4]]+boxLength[2];
            } else {
                tempPosZ[4] = oPosZ[ring5members[i][4]];
            }
            zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4];
            zVal /= 5;
            //outputer_xyz << setw(5) << "B" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
            outputer_xyz <<left<< setw(5) << "B" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";

            if (zVal <= slab_thickness && zVal >= -slab_thickness){
                zVal = -0.1;

                pov_out << "ring5(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";

                if (xVal < (-0.5*boxLength[0]+buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring5(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring5(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring5(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring5(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring5(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring5(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring5(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                    pov_out << "ring5(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                    pov_out << "ring5(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                }
            } 
        }

        red_val = "1.0";
        green_val = "0.0";
        blue_val = "0.0";
        for (i=0; i<ring6members.size(); i++){
            tempPosX[0]=oPosX[ring6members[i][0]];
            diffVal = oPosX[ring6members[i][1]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring6members[i][1]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring6members[i][1]]+boxLength[0];
            } else {
                tempPosX[1] = oPosX[ring6members[i][1]];
            }
            diffVal = oPosX[ring6members[i][2]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring6members[i][2]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring6members[i][2]]+boxLength[0];
            } else {
                tempPosX[2] = oPosX[ring6members[i][2]];
            }
            diffVal = oPosX[ring6members[i][3]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring6members[i][3]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring6members[i][3]]+boxLength[0];
            } else {
                tempPosX[3] = oPosX[ring6members[i][3]];
            }
            diffVal = oPosX[ring6members[i][4]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring6members[i][4]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring6members[i][4]]+boxLength[0];
            } else {
                tempPosX[4] = oPosX[ring6members[i][4]];
            }
            diffVal = oPosX[ring6members[i][5]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[5] = oPosX[ring6members[i][5]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[5] = oPosX[ring6members[i][5]]+boxLength[0];
            } else {
                tempPosX[5] = oPosX[ring6members[i][5]];
            }
            xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4]+tempPosX[5];
            xVal /= 6;

            tempPosY[0]=oPosY[ring6members[i][0]];
            diffVal = oPosY[ring6members[i][1]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring6members[i][1]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring6members[i][1]]+boxLength[1];
            } else {
                tempPosY[1] = oPosY[ring6members[i][1]];
            }
            diffVal = oPosY[ring6members[i][2]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring6members[i][2]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring6members[i][2]]+boxLength[1];
            } else {
                tempPosY[2] = oPosY[ring6members[i][2]];
            }
            diffVal = oPosY[ring6members[i][3]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring6members[i][3]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring6members[i][3]]+boxLength[1];
            } else {
                tempPosY[3] = oPosY[ring6members[i][3]];
            }
            diffVal = oPosY[ring6members[i][4]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring6members[i][4]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring6members[i][4]]+boxLength[1];
            } else {
                tempPosY[4] = oPosY[ring6members[i][4]];
            }
            diffVal = oPosY[ring6members[i][5]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[5] = oPosY[ring6members[i][5]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[5] = oPosY[ring6members[i][5]]+boxLength[1];
            } else {
                tempPosY[5] = oPosY[ring6members[i][5]];
            }
            yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4]+tempPosY[5];
            yVal /= 6;

            tempPosZ[0]=oPosZ[ring6members[i][0]];
            diffVal = oPosZ[ring6members[i][1]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring6members[i][1]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring6members[i][1]]+boxLength[2];
            } else {
                tempPosZ[1] = oPosZ[ring6members[i][1]];
            }
            diffVal = oPosZ[ring6members[i][2]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring6members[i][2]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring6members[i][2]]+boxLength[2];
            } else {
                tempPosZ[2] = oPosZ[ring6members[i][2]];
            }
            diffVal = oPosZ[ring6members[i][3]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring6members[i][3]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring6members[i][3]]+boxLength[2];
            } else {
                tempPosZ[3] = oPosZ[ring6members[i][3]];
            }
            diffVal = oPosZ[ring6members[i][4]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring6members[i][4]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring6members[i][4]]+boxLength[2];
            } else {
                tempPosZ[4] = oPosZ[ring6members[i][4]];
            }
            diffVal = oPosZ[ring6members[i][5]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[5] = oPosZ[ring6members[i][5]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[5] = oPosZ[ring6members[i][5]]+boxLength[2];
            } else {
                tempPosZ[5] = oPosZ[ring6members[i][5]];
            }
            zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4]+tempPosZ[5];
            zVal /= 6;
            //outputer_xyz << setw(5) << "C" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
            outputer_xyz <<left<< setw(5) << "C" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";

            if (zVal <= slab_thickness && zVal >= -slab_thickness){
                zVal = 0.0;

                pov_out << "ring6(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";

                if (xVal < (-0.5*boxLength[0]+buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring6(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring6(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring6(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring6(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring6(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring6(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring6(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                    pov_out << "ring6(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                    pov_out << "ring6(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                }
            } 
        }

        red_val = "1.0";
        green_val = "0.5";
        blue_val = "0.0";
        for (i=0; i<ring7members.size(); i++){
            tempPosX[0]=oPosX[ring7members[i][0]];
            diffVal = oPosX[ring7members[i][1]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring7members[i][1]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[1] = oPosX[ring7members[i][1]]+boxLength[0];
            } else {
                tempPosX[1] = oPosX[ring7members[i][1]];
            }
            diffVal = oPosX[ring7members[i][2]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring7members[i][2]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[2] = oPosX[ring7members[i][2]]+boxLength[0];
            } else {
                tempPosX[2] = oPosX[ring7members[i][2]];
            }
            diffVal = oPosX[ring7members[i][3]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring7members[i][3]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[3] = oPosX[ring7members[i][3]]+boxLength[0];
            } else {
                tempPosX[3] = oPosX[ring7members[i][3]];
            }
            diffVal = oPosX[ring7members[i][4]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring7members[i][4]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[4] = oPosX[ring7members[i][4]]+boxLength[0];
            } else {
                tempPosX[4] = oPosX[ring7members[i][4]];
            }
            diffVal = oPosX[ring7members[i][5]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[5] = oPosX[ring7members[i][5]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[5] = oPosX[ring7members[i][5]]+boxLength[0];
            } else {
                tempPosX[5] = oPosX[ring7members[i][5]];
            }
            diffVal = oPosX[ring7members[i][6]]-tempPosX[0];
            if (diffVal > 0.5*boxLength[0]){
                tempPosX[6] = oPosX[ring7members[i][6]]-boxLength[0];
            } else if (diffVal < -0.5*boxLength[0]){
                tempPosX[6] = oPosX[ring7members[i][6]]+boxLength[0];
            } else {
                tempPosX[6] = oPosX[ring7members[i][6]];
            }
            xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4]+tempPosX[5]+tempPosX[6];
            xVal /= 7;

            tempPosY[0]=oPosY[ring7members[i][0]];
            diffVal = oPosY[ring7members[i][1]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring7members[i][1]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[1] = oPosY[ring7members[i][1]]+boxLength[1];
            } else {
                tempPosY[1] = oPosY[ring7members[i][1]];
            }
            diffVal = oPosY[ring7members[i][2]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring7members[i][2]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[2] = oPosY[ring7members[i][2]]+boxLength[1];
            } else {
                tempPosY[2] = oPosY[ring7members[i][2]];
            }
            diffVal = oPosY[ring7members[i][3]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring7members[i][3]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[3] = oPosY[ring7members[i][3]]+boxLength[1];
            } else {
                tempPosY[3] = oPosY[ring7members[i][3]];
            }
            diffVal = oPosY[ring7members[i][4]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring7members[i][4]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[4] = oPosY[ring7members[i][4]]+boxLength[1];
            } else {
                tempPosY[4] = oPosY[ring7members[i][4]];
            }
            diffVal = oPosY[ring7members[i][5]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[5] = oPosY[ring7members[i][5]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[5] = oPosY[ring7members[i][5]]+boxLength[1];
            } else {
                tempPosY[5] = oPosY[ring7members[i][5]];
            }
            diffVal = oPosY[ring7members[i][6]]-tempPosY[0];
            if (diffVal > 0.5*boxLength[1]){
                tempPosY[6] = oPosY[ring7members[i][6]]-boxLength[1];
            } else if (diffVal < -0.5*boxLength[1]){
                tempPosY[6] = oPosY[ring7members[i][6]]+boxLength[1];
            } else {
                tempPosY[6] = oPosY[ring7members[i][6]];
            }
            yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4]+tempPosY[5]+tempPosY[6];
            yVal /= 7;

            tempPosZ[0]=oPosZ[ring7members[i][0]];
            diffVal = oPosZ[ring7members[i][1]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring7members[i][1]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[1] = oPosZ[ring7members[i][1]]+boxLength[2];
            } else {
                tempPosZ[1] = oPosZ[ring7members[i][1]];
            }
            diffVal = oPosZ[ring7members[i][2]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring7members[i][2]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[2] = oPosZ[ring7members[i][2]]+boxLength[2];
            } else {
                tempPosZ[2] = oPosZ[ring7members[i][2]];
            }
            diffVal = oPosZ[ring7members[i][3]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring7members[i][3]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[3] = oPosZ[ring7members[i][3]]+boxLength[2];
            } else {
                tempPosZ[3] = oPosZ[ring7members[i][3]];
            }
            diffVal = oPosZ[ring7members[i][4]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring7members[i][4]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[4] = oPosZ[ring7members[i][4]]+boxLength[2];
            } else {
                tempPosZ[4] = oPosZ[ring7members[i][4]];
            }
            diffVal = oPosZ[ring7members[i][5]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[5] = oPosZ[ring7members[i][5]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[5] = oPosZ[ring7members[i][5]]+boxLength[2];
            } else {
                tempPosZ[5] = oPosZ[ring7members[i][5]];
            }
            diffVal = oPosZ[ring7members[i][6]]-tempPosZ[0];
            if (diffVal > 0.5*boxLength[2]){
                tempPosZ[6] = oPosZ[ring7members[i][6]]-boxLength[2];
            } else if (diffVal < -0.5*boxLength[2]){
                tempPosZ[6] = oPosZ[ring7members[i][6]]+boxLength[2];
            } else {
                tempPosZ[6] = oPosZ[ring7members[i][6]];
            }
            zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4]+tempPosZ[5]+tempPosZ[6];
            zVal /= 7;
            //outputer_xyz << setw(5) << "N" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
            outputer_xyz <<left<< setw(5) << "N" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";

            if (zVal <= slab_thickness && zVal >= -slab_thickness){
                zVal = -0.2;

                pov_out << "ring7(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";

                if (xVal < (-0.5*boxLength[0]+buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring7(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring7(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring7(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                    if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring7(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring7(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        pov_out << "ring7(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else {
                        pov_out << "ring7(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                    pov_out << "ring7(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                    pov_out << "ring7(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                }
            } 
        }

        x_frame = int(scale_factor * boxLength[0]);
        y_frame = int(scale_factor * boxLength[1]);
        xFrameInt.str("");
        yFrameInt.str("");
        xFrameInt << x_frame;
        yFrameInt << y_frame;
        xFrameString = xFrameInt.str();
        yFrameString = yFrameInt.str();
        outputer_pov << "povray -w" + xFrameString + " -h" + yFrameString + " +a0.1 -D " + povName2 + "_" + frameCountString + ".pov\n";
        frameCount++;

        // close the pov file handle
        pov_out.close();
        povDistOut.close();

        // clear reused vectors
        oPosX.clear();
        oPosY.clear(); 
        oPosZ.clear();
        hPosX.clear();
        hPosY.clear();
        hPosZ.clear();
        neighborListIndex.clear();
        neighborList.clear();
        hbondListIndex.clear();
        hbondList.clear();
        ring3members.clear(); 
        ring4members.clear();
        ring5members.clear();
        ring6members.clear();
        ring7members.clear();
        hbondCount = 0;
        lastRing3 = 0;
        lastRing4 = 0;
        lastRing5 = 0;
        lastRing6 = 0;
        lastRing7 = 0;
        lastRing8 = 0;
        totalRings = 0;

        // now we read in the next frame
        inputer.getline(inLine,999,'\n');
        inputer.getline(inLine,999,'\n');
        if (!inputer.eof()){
            token = strtok(inLine,delimit);
            strcpy(inValue,token);
            nAtoms = atoi(inValue);
            atomCount.push_back(nAtoms);
        }
    }

    //    delete calculator;
    //    delete povrayObjects;

    cout << "\nResults written to " << strungName << "\n\n";
    return 0;
}

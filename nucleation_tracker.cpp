/*
   Program to calculate the hbond rings in aqueous environments.
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
#include <iterator>
#include <algorithm>
#include "cmdline.h"

using namespace std;

// some global parameters
double distance_tol = 3.5;   // hbond distance tolerance (in Ångström)
double distance_tol2 = distance_tol*distance_tol; // hbond square distance tolerance
double distance_tol2_long = 2.0 * distance_tol*distance_tol; // long square distance tolerance
double angle_tol = 30;       // hbond angle tolerance (in degrees)
double angle_tol_rad = angle_tol * 3.1415926536 / 180.0; // hbond angle tolerance in radians
double OH_length = 0.926;    // OH bond length (in Ångström)
double OH_length_i = 1.0 / OH_length; // OH bond length (in Ångström)
double neighbor_tol = 0.5; // neighbor distance tolerance window triggering a H-bonding screening (in Ångström)
double rad_to_degree = 180.0 * 3.1415926536; // convert radians to degrees
double ring3_dot_size = 0.5; // central ring dot size
double ring4_dot_size = 0.6; // central ring dot size
double ring5_dot_size = 0.7; // central ring dot size
double ring6_dot_size = 0.8; // central ring dot size
double ring7_dot_size = 0.9; // central ring dot size
double ring8_dot_size = 1.0; // central ring dot size
double ring9_dot_size = 1.1; // central ring dot size
double ring10_dot_size = 1.2; // central ring dot size
double scale_factor = 20;    // scale the pixels!
double buffer_size = 2;      // the imaging buffer size
double slab_thickness = 100;  // 0.5*slab thickness
double transparency = 0.6;   // the transparency of ring dots
double outline = 0.04;
double axes_width = 0.08;

class Calcs {
    public:
        Calcs();
        ~Calcs();
        //int IsWat3HBond(double (&pos1)[9], double (&pos2)[9], double (&box)[3]);
        int IsWat3HBond(double (&pos1)[9], double (&pos2)[9], double (&hmatrix)[3][3], double (&ihmatrix)[3][3], double (&bondvec)[3], int &i, int &j);
        int IsWat5HBond(double (&pos1)[15], double (&pos2)[15], double (&hmatrix)[3][3], double (&ihmatrix)[3][3], double (&bondvec)[3], int &i, int &j);
    private:
        bool printFlag;
};

Calcs::Calcs() {
    printFlag = true;
}

Calcs::~Calcs() {
}

int Calcs::IsWat3HBond(double (&pos1)[9], double (&pos2)[9], double (&hmatrix)[3][3], double (&ihmatrix)[3][3], double (&bondvec)[3], int &i, int &j){
    // a function to identify if there is an HBond between 3 point water models
    int isHBond = 0;
    //double r, r2, ri, xVal, yVal, zVal, tempAngle;
    int k;
    double r, r2, ri, tempAngle;
    double posVec[3], vec1[3], vec2[3];
    double scaledVec[3];
    double dotProduct;

    // the 9 positions are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z)
    // first, we do the O-O vector
    posVec[0] = pos2[0]-pos1[0];
    posVec[1] = pos2[1]-pos1[1];
    posVec[2] = pos2[2]-pos1[2];

    // do vector wrapping of periodic boundary conditions
    for (k=0; k<3; k++){
        scaledVec[k] = posVec[k] * ihmatrix[k][k];
        scaledVec[k] -= round(scaledVec[k]);
        posVec[k] = scaledVec[k] * hmatrix[k][k];
    }

    // here, we take care of normalization
    // r2 = xVal*xVal + yVal*yVal + zVal*zVal;
    r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
    r = sqrt(r2);
    ri = 1.0 / r;

    vec1[0] = posVec[0] * ri;
    vec1[1] = posVec[1] * ri;
    vec1[2] = posVec[2] * ri;

    // the O-H1
    posVec[0] = pos1[3]-pos1[0];
    posVec[1] = pos1[4]-pos1[1];
    posVec[2] = pos1[5]-pos1[2];

    // do vector wrapping of periodic boundary conditions
    for (k=0; k<3; k++){
        scaledVec[k] = posVec[k] * ihmatrix[k][k];
        scaledVec[k] -= round(scaledVec[k]);
        posVec[k] = scaledVec[k] * hmatrix[k][k];
    }

    // here, we take care of normalization
    // r2 = xVal*xVal + yVal*yVal + zVal*zVal;
    r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
    r = sqrt(r2);
    ri = 1.0 / r;

    vec2[0] = posVec[0] * ri;
    vec2[1] = posVec[1] * ri;
    vec2[2] = posVec[2] * ri;

    dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
    if (dotProduct > 1.0){
        tempAngle = 0.0;
    } else {
        tempAngle = acos(dotProduct);
    }

    if (tempAngle <= angle_tol_rad){
        isHBond = 1;
    } else {
        // continue down the rabbit hole

        // the O-H2
        posVec[0] = pos1[6]-pos1[0];
        posVec[1] = pos1[7]-pos1[1];
        posVec[2] = pos1[8]-pos1[2];

        // do vector wrapping of periodic boundary conditions
        for (k=0; k<3; k++){
            scaledVec[k] = posVec[k] * ihmatrix[k][k];
            scaledVec[k] -= round(scaledVec[k]);
            posVec[k] = scaledVec[k] * hmatrix[k][k];
        }

        // here, we take care of normalization
        // r2 = xVal*xVal + yVal*yVal + zVal*zVal;
        r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
        r = sqrt(r2);
        ri = 1.0 / r;

        vec2[0] = posVec[0] * ri;
        vec2[1] = posVec[1] * ri;
        vec2[2] = posVec[2] * ri;

        dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
        if (dotProduct > 1.0){
            tempAngle = 0.0;
        } else {
            tempAngle = acos(dotProduct);
        }

        if (tempAngle <= angle_tol_rad){
            isHBond = 1;
        } else {
            // continue further down the rabbit hole

            // flip the O-O vector
            vec1[0] *= -1;
            vec1[1] *= -1;
            vec1[2] *= -1;

            // the O-H1
            posVec[0] = pos2[3]-pos2[0];
            posVec[1] = pos2[4]-pos2[1];
            posVec[2] = pos2[5]-pos2[2];

            // do vector wrapping of periodic boundary conditions
            for (k=0; k<3; k++){
                scaledVec[k] = posVec[k] * ihmatrix[k][k];
                scaledVec[k] -= round(scaledVec[k]);
                posVec[k] = scaledVec[k] * hmatrix[k][k];
            }

            // here, we take care of normalization
            // r2 = xVal*xVal + yVal*yVal + zVal*zVal;
            r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
            r = sqrt(r2);
            ri = 1.0 / r;

            vec2[0] = posVec[0] * ri;
            vec2[1] = posVec[1] * ri;
            vec2[2] = posVec[2] * ri;

            dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
            if (dotProduct > 1.0){
                tempAngle = 0.0;
            } else {
                tempAngle = acos(dotProduct);
            }

            if (tempAngle <= angle_tol_rad){
                isHBond = 1;
            } else {
                // last level... I promise...

                // the O-H2
                posVec[0] = pos2[6]-pos2[0];
                posVec[1] = pos2[7]-pos2[1];
                posVec[2] = pos2[8]-pos2[2];

                // do vector wrapping of periodic boundary conditions
                for (k=0; k<3; k++){
                    scaledVec[k] = posVec[k] * ihmatrix[k][k];
                    scaledVec[k] -= round(scaledVec[k]);
                    posVec[k] = scaledVec[k] * hmatrix[k][k];
                }

                // here, we take care of normalization
                // r2 = xVal*xVal + yVal*yVal + zVal*zVal;
                r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
                r = sqrt(r2);
                ri = 1.0 / r;

                vec2[0] = posVec[0] * ri;
                vec2[1] = posVec[1] * ri;
                vec2[2] = posVec[2] * ri;

                dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
                if (dotProduct > 1.0){
                    tempAngle = 0.0;
                } else {
                    tempAngle = acos(dotProduct);
                }

                if (tempAngle <= angle_tol_rad){
                    isHBond = 1;
                }
            }
        }
    }

    if (isHBond){
        bondvec[0] = vec1[0];
        bondvec[1] = vec1[1];
        bondvec[2] = vec1[2];
    }

    return isHBond;
}

int Calcs::IsWat5HBond(double (&pos1)[15], double (&pos2)[15], double (&hmatrix)[3][3], double (&ihmatrix)[3][3], double (&bondvec)[3], int &i, int &j){
    // a function to identify if there is an HBond between 3 point water models
    int isHBond = 0;
    int isHBondStep1 = 0;
    int isHBondStep2 = 0;
    //double r, r2, ri, xVal, yVal, zVal, tempAngle;
    int k;
    double r, r2, ri, tempAngle;
    double posVec[3], vec1[3], vec2[3];
    double scaledVec[3];
    double dotProduct;

    // the 15 positions are sorted by O(x,y,z), H1(x,y,z), H2(x,y,z), H3(x,y,z), H4(x,y,z)
    // first, we do the O-O vector
    posVec[0] = pos2[0]-pos1[0];
    posVec[1] = pos2[1]-pos1[1];
    posVec[2] = pos2[2]-pos1[2];

    // do vector wrapping of periodic boundary conditions
    for (k=0; k<3; k++){
        scaledVec[k] = posVec[k] * ihmatrix[k][k];
        scaledVec[k] -= round(scaledVec[k]);
        posVec[k] = scaledVec[k] * hmatrix[k][k];
    }

    // here, we take care of normalization
    r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
    r = sqrt(r2);
    ri = 1.0 / r;

    vec1[0] = posVec[0] * ri;
    vec1[1] = posVec[1] * ri;
    vec1[2] = posVec[2] * ri;

    // the O-H1
    posVec[0] = pos1[3]-pos1[0];
    posVec[1] = pos1[4]-pos1[1];
    posVec[2] = pos1[5]-pos1[2];

    // do vector wrapping of periodic boundary conditions
    for (k=0; k<3; k++){
        scaledVec[k] = posVec[k] * ihmatrix[k][k];
        scaledVec[k] -= round(scaledVec[k]);
        posVec[k] = scaledVec[k] * hmatrix[k][k];
    }

    // here, we take care of normalization
    r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
    r = sqrt(r2);
    ri = 1.0 / r;

    vec2[0] = posVec[0] * ri;
    vec2[1] = posVec[1] * ri;
    vec2[2] = posVec[2] * ri;

    dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
    if (dotProduct > 1.0){
        tempAngle = 0.0;
    } else {
        tempAngle = acos(dotProduct);
    }

    if (tempAngle <= angle_tol_rad){
        isHBondStep1 = 1;
    } else {
        // continue down the rabbit hole

        // the O-H2
        posVec[0] = pos1[6]-pos1[0];
        posVec[1] = pos1[7]-pos1[1];
        posVec[2] = pos1[8]-pos1[2];

        // do vector wrapping of periodic boundary conditions
        for (k=0; k<3; k++){
            scaledVec[k] = posVec[k] * ihmatrix[k][k];
            scaledVec[k] -= round(scaledVec[k]);
            posVec[k] = scaledVec[k] * hmatrix[k][k];
        }

        // here, we take care of normalization
        r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
        r = sqrt(r2);
        ri = 1.0 / r;

        vec2[0] = posVec[0] * ri;
        vec2[1] = posVec[1] * ri;
        vec2[2] = posVec[2] * ri;

        dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
        if (dotProduct > 1.0){
            tempAngle = 0.0;
        } else {
            tempAngle = acos(dotProduct);
        }

        if (tempAngle <= angle_tol_rad){
            isHBondStep1 = 1;
        } else {
            // continue down the rabbit hole #2

            // the O-H3
            posVec[0] = pos1[9]-pos1[0];
            posVec[1] = pos1[10]-pos1[1];
            posVec[2] = pos1[11]-pos1[2];

            // do vector wrapping of periodic boundary conditions
            for (k=0; k<3; k++){
                scaledVec[k] = posVec[k] * ihmatrix[k][k];
                scaledVec[k] -= round(scaledVec[k]);
                posVec[k] = scaledVec[k] * hmatrix[k][k];
            }

            // here, we take care of normalization
            r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
            r = sqrt(r2);
            ri = 1.0 / r;

            vec2[0] = posVec[0] * ri;
            vec2[1] = posVec[1] * ri;
            vec2[2] = posVec[2] * ri;

            dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
            if (dotProduct > 1.0){
                tempAngle = 0.0;
            } else {
                tempAngle = acos(dotProduct);
            }

            if (tempAngle <= angle_tol_rad){
                isHBondStep1 = 1;
            } else {
                // continue down the rabbit hole #3

                // the O-H4
                posVec[0] = pos1[12]-pos1[0];
                posVec[1] = pos1[13]-pos1[1];
                posVec[2] = pos1[14]-pos1[2];

                // do vector wrapping of periodic boundary conditions
                for (k=0; k<3; k++){
                    scaledVec[k] = posVec[k] * ihmatrix[k][k];
                    scaledVec[k] -= round(scaledVec[k]);
                    posVec[k] = scaledVec[k] * hmatrix[k][k];
                }

                // here, we take care of normalization
                r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
                r = sqrt(r2);
                ri = 1.0 / r;

                vec2[0] = posVec[0] * ri;
                vec2[1] = posVec[1] * ri;
                vec2[2] = posVec[2] * ri;

                dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
                if (dotProduct > 1.0){
                    tempAngle = 0.0;
                } else {
                    tempAngle = acos(dotProduct);
                }

                if (tempAngle <= angle_tol_rad){
                    isHBondStep1 = 1;
                }
            }
        }
    }

    if (isHBondStep1 == 1){
        // now check the other water orientation
        // flip the O-O vector
        vec1[0] *= -1;
        vec1[1] *= -1;
        vec1[2] *= -1;

        // the O-H1
        posVec[0] = pos2[3]-pos2[0];
        posVec[1] = pos2[4]-pos2[1];
        posVec[2] = pos2[5]-pos2[2];

        // do vector wrapping of periodic boundary conditions
        for (k=0; k<3; k++){
            scaledVec[k] = posVec[k] * ihmatrix[k][k];
            scaledVec[k] -= round(scaledVec[k]);
            posVec[k] = scaledVec[k] * hmatrix[k][k];
        }

        // here, we take care of normalization
        r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
        r = sqrt(r2);
        ri = 1.0 / r;

        vec2[0] = posVec[0] * ri;
        vec2[1] = posVec[1] * ri;
        vec2[2] = posVec[2] * ri;

        dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
        if (dotProduct > 1.0){
            tempAngle = 0.0;
        } else {
            tempAngle = acos(dotProduct);
        }

        if (tempAngle <= angle_tol_rad){
            isHBondStep2 = 1;
        } else {
            // last level... I promise...

            // the O-H2
            posVec[0] = pos2[6]-pos2[0];
            posVec[1] = pos2[7]-pos2[1];
            posVec[2] = pos2[8]-pos2[2];

            // do vector wrapping of periodic boundary conditions
            for (k=0; k<3; k++){
                scaledVec[k] = posVec[k] * ihmatrix[k][k];
                scaledVec[k] -= round(scaledVec[k]);
                posVec[k] = scaledVec[k] * hmatrix[k][k];
            }

            // here, we take care of normalization
            r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
            r = sqrt(r2);
            ri = 1.0 / r;

            vec2[0] = posVec[0] * ri;
            vec2[1] = posVec[1] * ri;
            vec2[2] = posVec[2] * ri;

            dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
            if (dotProduct > 1.0){
                tempAngle = 0.0;
            } else {
                tempAngle = acos(dotProduct);
            }

            if (tempAngle <= angle_tol_rad){
                isHBondStep2 = 1;
            } else {
                // another level...

                // the O-H3
                posVec[0] = pos2[9]-pos2[0];
                posVec[1] = pos2[10]-pos2[1];
                posVec[2] = pos2[11]-pos2[2];

                // do vector wrapping of periodic boundary conditions
                for (k=0; k<3; k++){
                    scaledVec[k] = posVec[k] * ihmatrix[k][k];
                    scaledVec[k] -= round(scaledVec[k]);
                    posVec[k] = scaledVec[k] * hmatrix[k][k];
                }

                // here, we take care of normalization
                r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
                r = sqrt(r2);
                ri = 1.0 / r;

                vec2[0] = posVec[0] * ri;
                vec2[1] = posVec[1] * ri;
                vec2[2] = posVec[2] * ri;

                dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
                if (dotProduct > 1.0){
                    tempAngle = 0.0;
                } else {
                    tempAngle = acos(dotProduct);
                }

                if (tempAngle <= angle_tol_rad){
                    isHBondStep2 = 1;
                } else {
                    // the last level...

                    // the O-H4
                    posVec[0] = pos2[12]-pos2[0];
                    posVec[1] = pos2[13]-pos2[1];
                    posVec[2] = pos2[14]-pos2[2];

                    // do vector wrapping of periodic boundary conditions
                    for (k=0; k<3; k++){
                        scaledVec[k] = posVec[k] * ihmatrix[k][k];
                        scaledVec[k] -= round(scaledVec[k]);
                        posVec[k] = scaledVec[k] * hmatrix[k][k];
                    }

                    // here, we take care of normalization
                    r2 = posVec[0]*posVec[0] + posVec[1]*posVec[1] + posVec[2]*posVec[2];
                    r = sqrt(r2);
                    ri = 1.0 / r;

                    vec2[0] = posVec[0] * ri;
                    vec2[1] = posVec[1] * ri;
                    vec2[2] = posVec[2] * ri;

                    dotProduct = ((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]));
                    if (dotProduct > 1.0){
                        tempAngle = 0.0;
                    } else {
                        tempAngle = acos(dotProduct);
                    }

                    if (tempAngle <= angle_tol_rad){
                        isHBondStep2 = 1;
                    }
                }
            }
        }

        // final determination of if we have a HBond
        if (isHBondStep2 == 1){
            isHBond = 1;
        }
    }

    if (isHBond){
        bondvec[0] = vec1[0];
        bondvec[1] = vec1[1];
        bondvec[2] = vec1[2];
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
        void printRing8(ofstream& outputFile, double (&boxLength)[3]);
	void printRing9(ofstream& outputFile, double (&boxLength)[3]);
        void printRing10(ofstream& outputFile, double (&boxLength)[3]);
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

void PovObjects::printRing8(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring8 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n "
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring8_dot_size << "\n"
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

void PovObjects::printRing9(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring9 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring9_dot_size << "\n"
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

void PovObjects::printRing10(ofstream& outputFile, double (&boxLength)[3]){
    outputFile << "#macro ring10 (center_x, center_y, center_z, outerRed, outerGreen, outerBlue, dot_transparency)\n " 
        "  intersection{\n"
        "    sphere{\n"
        "      < center_x, center_y, center_z >,\n    " << ring10_dot_size << "\n"
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

struct Neighbor{
    int id;
    double distance;
    double vectorX;
    double vectorY;
    double vectorZ;
};

bool compareByDistance(const Neighbor &a, const Neighbor &b)
{
    return a.distance < b.distance;
};

int main(int argc, char **argv) {
    bool trigger;
    bool isPos;
    bool fiveAtomWater = false;
    char povDistFileName[200], povFileName[200], trajFileName[200], fileName[200], inLine[1000], inValue[200], tag[100], tetraPDBFileName[200];
    char *token;
    const char *file;
    const char *delimit = " \t\n";
    const char *period = ".";
    const char *OWAtom = "   OW";
    const char *HWAtom1 = "  HW1";
    const char *HWAtom2 = "  HW2";
    const char *HWAtom3 = "  HW3";
    const char *HWAtom4 = "  HW4";
    const char *OAtom = "    O";
    const char *HAtom1 = "   H1";
    const char *HAtom2 = "   H2";
    const char *HAtom3 = "   H3";
    const char *HAtom4 = "   H4";
    const char *HAtom = "    H";
    int count = 0;
    int tmp_count, loaded;
    int i, j, k, l, m, n, o, p, q, r, s, nAtoms, waterModel, waterNumber;
    int oCount, hCount, mCount, totalWaterCount, hbondCount, neighborCount, frameCount, largeNeighborCount;
    int index2, index3, index4, index5, index6, index7, index8, index9, index10, index11;
    int lastRing3, lastRing4, lastRing5, lastRing6, lastRing7, lastRing8, lastRing9, lastRing10;
    int x_frame, y_frame;
    int ringOutOpt;
    int maxRingOpt;
    int directionalityOpt = 0;
    int povrayOutOpt = 0;
    int ringTrajOutOpt = 0;
    int tetraPDBOutOpt = 0;
    int totalRings;
    int index, test_slot;
    int temp_index;
    int is_donor = 0;
    double cosine_sum, dotProduct, inv_nearestNeighborMag1, inv_nearestNeighborMag2, unitVec1x, unitVec1y, unitVec1z, unitVec2x, unitVec2y, unitVec2z;
    double tetrahedralParam, avg_tetrahedrality, tet_square, stdev_tetrahedrality;
    double xVal, yVal, zVal;
    double temp_dist_inv, temp_dist, temp_dist2;
    double temp_value;
    double temp_val_x, temp_val_y, temp_val_z;
    double ionPosition[3], position[3];
    double water1[9], water2[9];
    double water15[15], water25[15];
    double pos[3], scaled[3];
    double boxLength[3];
    double hmat[3][3], invhmat[3][3];
    double determinant, invdeterminant;
    double tempPosX[10], tempPosY[10], tempPosZ[10];
    double tempHBVec[3], tempHBVec2[3], tempScaledHBVec[3];
    double diffVal;
    double ringProbabilities[11];
    double occupancy;
    Neighbor closest_neighbor, other_neighbor;
    vector<int> atomCount;
    vector<int> neighborListIndex;
    vector<int> largeNeighborListIndex;
    vector<int> hbondListIndex;
    vector<int> tempVec;
    vector<int> bondSequence;
    vector<int> angleSequence;
    vector<int> dihedralSequence;
    vector<double> oPosX, oPosY, oPosZ;
    vector<double> hPosX, hPosY, hPosZ;
    vector<double> tempList, tempVecX, tempVecY, tempVecZ;
    vector<double> tetrahedrality;
    vector<vector<int> > neighborList;
    vector<vector<Neighbor> > largeNeighborList;
    vector<Neighbor> singleSort;
    vector<vector<Neighbor> > nearestNeighborList;
    vector<vector<int> > hbondList;
    vector<vector<int> > ring3members;
    vector<vector<int> > ring4members;
    vector<vector<int> > ring5members;
    vector<vector<int> > ring6members;
    vector<vector<int> > ring7members;
    vector<vector<int> > ring8members;
    vector<vector<int> > ring9members;
    vector<vector<int> > ring10members;
    vector<vector<double> > hbondVecX;
    vector<vector<double> > hbondVecY;
    vector<vector<double> > hbondVecZ;
    vector<vector<double> > nearestNeighborVecX;
    vector<vector<double> > nearestNeighborVecY;
    vector<vector<double> > nearestNeighborVecZ;
    vector<vector<double> > nearestNeighborMag;
    string strungName, strungName2, strungName3, strungName4;
    string povName;
    string povName2;
    string tetraPDBDistName;
    string povDirName;
    string povDistName;
    string red_val, green_val, blue_val;
    string lineString;
    string shortString;
    string frameCountString, xFrameString, yFrameString;
    stringstream frameInt, xFrameInt, yFrameInt;
    ofstream pov_out;
    ofstream outputer_pov;
    ofstream povDistOut;
    ofstream outputer_xyz;
    ofstream outputer_PDB;

    gengetopt_args_info args_info;

    /* let's call our cmdline parser */
    if (cmdline_parser (argc, argv, &args_info) != 0)
        exit(1) ;

    if (args_info.directionality_flag == 1){
        directionalityOpt = 1;
    }
    if (args_info.povray_flag == 1){
        povrayOutOpt = 1;
    }
    if (args_info.ring_trajectory_flag == 1){
        ringTrajOutOpt = 1;
    }
    if (args_info.tetra_pdb_flag == 1){
        tetraPDBOutOpt = 1;
    }
    maxRingOpt = args_info.max_ring_arg;
    ringOutOpt = args_info.closure_method_arg;

    // Now try opening the file
    ifstream prayer(args_info.input_file_arg);

    // Make sure the file exists
    if (!prayer) {
        cout << "Unable to open " << args_info.input_file_arg << " for reading.\n";
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
    //strungName = argv[2];
    strungName = args_info.input_file_arg;
    file = strungName.c_str();
    strcpy(fileName, file);
    token = strtok(fileName, period);

    // Build readers and writers as necessary given options
    // ifstream inputer(args_info.input_file_arg);
    // ofstream outputer(fileName);

    if (ringTrajOutOpt){ 
        strcpy(trajFileName, token);
        strcat(trajFileName,"_ringtrj.xyz");
        strungName2 = trajFileName;
        outputer_xyz.open(trajFileName);
    }

    if (povrayOutOpt){
        strcpy(povFileName, token);
        strcpy(povDistFileName, token);
        povName2 = povFileName;
        povDistName = povDistFileName;
        strcat(povFileName,"_pov.txt");
        strungName = fileName;
        povDirName = "pov_files/";
        strungName3 = povDirName + povFileName;
        outputer_pov.open(strungName3);
    }

    if (tetraPDBOutOpt){
        strcpy(tetraPDBFileName, token);
        strcat(tetraPDBFileName,"_tetra.pdb");
        strungName4 = tetraPDBFileName;
		outputer_PDB.open(tetraPDBFileName);
    }

    // finally appending "_nuc_info.txt" to a file
    strcpy(fileName, token);
    strcat(fileName,"_nuc_info.txt");


    // Build readers and writers as necessary given options
    ifstream inputer(args_info.input_file_arg);
    ofstream outputer(fileName);
 
    cmdline_parser_free (&args_info); /* release gotten options allocated memory */

    // initialize our calculator
    Calcs *calculator = new Calcs();


    // Read the .gro file and load the positions
    cout << "\nLoading and processing trajectory...\n";

    outputer << setw(8) << "# frame" <<  setw(8) << "HBonds"  << setw(8) << "<q>" << setw(8) << "stdev" << setw(8) << "3";
    if (maxRingOpt > 3){
        outputer << setw(8) << "4";
        if (maxRingOpt > 4){
            outputer << setw(8) << "5";
            if (maxRingOpt > 5){
                outputer << setw(8) << "6";
                if (maxRingOpt > 6){
                    outputer << setw(8) << "7";
                    if (maxRingOpt > 7){
                        outputer << setw(8) << "8";
                    	if (maxRingOpt > 8){
                            outputer << setw(8) << "9";
                    	    if (maxRingOpt > 9){
                    	    	outputer << setw(8) << "10";
                    	    }
                   	}
                    }
                }
            }
        }
    }
    outputer << "\n";

    frameCount = 1;
    totalWaterCount = 0;
    hbondCount = 0;
    neighborCount = 0;
    largeNeighborCount = 0;
    lastRing3 = 0;
    lastRing4 = 0;
    lastRing5 = 0;
    lastRing6 = 0;
    lastRing7 = 0;
    lastRing8 = 0;
    lastRing9 = 0;
    lastRing10 = 0;
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
                totalWaterCount++;
                // test if it is a hydrogen atom and load in the appropriate vector
            } else if (!strcmp(HAtom1, token) || !strcmp(HWAtom1, token) || !strcmp(HAtom, token)){
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
            } else if (!strcmp(HAtom3, token) || !strcmp(HWAtom3, token)){
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
            } else if (!strcmp(HAtom4, token) || !strcmp(HWAtom4, token)){
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

        hmat[0][0] = boxLength[0];
        hmat[1][1] = boxLength[1];
        hmat[2][2] = boxLength[2];

        // and load the rest of the box cell matrix (hmat) using gromacs form.
        token = strtok(NULL,delimit);
        if (token != NULL){
            strcpy(inValue,token);
            hmat[0][1] = atof(inValue) * 10.0;
            token = strtok(NULL,delimit);
            strcpy(inValue,token);
            hmat[0][2] = atof(inValue) * 10.0;
            token = strtok(NULL,delimit);
            strcpy(inValue,token);
            hmat[1][0] = atof(inValue) * 10.0;
            token = strtok(NULL,delimit);
            strcpy(inValue,token);
            hmat[1][2] = atof(inValue) * 10.0;
            token = strtok(NULL,delimit);
            strcpy(inValue,token);
            hmat[2][0] = atof(inValue) * 10.0;
            token = strtok(NULL,delimit);
            strcpy(inValue,token);
            hmat[2][1] = atof(inValue) * 10.0;
        } else {
            hmat[0][1] = 0;
            hmat[0][2] = 0;
            hmat[1][0] = 0;
            hmat[1][2] = 0;
            hmat[2][0] = 0;
            hmat[2][1] = 0;
        }

        // let's determine the inverse of the hmat
        // first the determinant...
        determinant = 0;
        for(i=0; i<3; i++){
            determinant = determinant + (hmat[0][i] * (hmat[1][(i+1)%3] * hmat[2][(i+2)%3] - hmat[1][(i+2)%3] * hmat[2][(i+1)%3]));
        }
        invdeterminant = 1.0 / determinant;
        // now the double loop inverse...
        for(i=0; i<3; i++){
            for(j=0; j<3; j++){
                invhmat[i][j] = ((hmat[(j+1)%3][(i+1)%3] * hmat[(j+2)%3][(i+2)%3]) - (hmat[(j+1)%3][(i+2)%3] * hmat[(j+2)%3][(i+1)%3])) * invdeterminant;
            }
        }

        //cerr << "\n" << oPosX.size() << " : number of water molecules identified\n";
        //cerr << "\n" << hPosX.size() << " : number of water H atoms identified\n";
        if (4*oPosX.size() == hPosX.size()){
            fiveAtomWater = true;
        } else {
            if (oPosX.size() <= 0 || 2*oPosX.size() != hPosX.size()){
                cerr << "Error: Cannot identify any water molecules.  Best to use atomtypes   of 'OW', 'HW1', and 'HW2'.\n";
                return 0;
            }
        }
        // okay, let's do some calculations on this frame

        // we loop over the waters in the frame to determine if they
        // are neighbors
        for (i=0; i<3; i++){
            tempHBVec[i] = 0;
        }

        Neighbor blank;
        blank.id = -1;
        blank.distance = 1000000;
        blank.vectorX = 0;
        blank.vectorY = 0;
        blank.vectorZ = 0;

        for (i=0; i<oPosX.size(); i++){
            // can have up to 50 neighbors
            // vec2D.push_back(std::vector<int>(4, 11));
            neighborListIndex.push_back(0);
            neighborList.push_back(vector<int>(50,-1));
            largeNeighborListIndex.push_back(0);
            largeNeighborList.push_back(vector<Neighbor>(100,blank));
            nearestNeighborList.push_back(vector<Neighbor>(4,blank));

            // can have up to 50 HBonds
            hbondListIndex.push_back(0);
            hbondList.push_back(vector<int>(50,-1));
            hbondVecX.push_back(vector<double>(50,0));
            hbondVecY.push_back(vector<double>(50,0));
            hbondVecZ.push_back(vector<double>(50,0));
        }

        for (i=0; i<oPosX.size()-1; i++){
            // just use o positions
            // load the reference water
            water1[0] = oPosX[i];
            water1[1] = oPosY[i];
            water1[2] = oPosZ[i];

            for (j=i+1; j<oPosX.size(); j++){
                // load the water
                water2[0] = oPosX[j];
                water2[1] = oPosY[j];
                water2[2] = oPosZ[j];

                // calculate if they are neighbors
                pos[0] = water2[0]-water1[0];
                pos[1] = water2[1]-water1[1];
                pos[2] = water2[2]-water1[2];

                // do vector wrapping of periodic boundary conditions
                for (k=0; k<3; k++){
                    scaled[k] = pos[k] * invhmat[k][k];
                    scaled[k] -= round(scaled[k]);
                    pos[k] = scaled[k] * hmat[k][k];
                }

                temp_dist2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
                if(temp_dist2 <= distance_tol2_long){
                    temp_dist = sqrt(temp_dist2);
                    largeNeighborList[i][largeNeighborListIndex[i]].id = j;
                    largeNeighborList[i][largeNeighborListIndex[i]].distance = temp_dist;
                    largeNeighborList[i][largeNeighborListIndex[i]].vectorX = pos[0];
                    largeNeighborList[i][largeNeighborListIndex[i]].vectorY = pos[1];
                    largeNeighborList[i][largeNeighborListIndex[i]].vectorZ = pos[2];
                    largeNeighborListIndex[i]++;
                    largeNeighborList[j][largeNeighborListIndex[j]].id = i;
                    largeNeighborList[j][largeNeighborListIndex[j]].distance = temp_dist;
                    largeNeighborList[j][largeNeighborListIndex[j]].vectorX = -pos[0];
                    largeNeighborList[j][largeNeighborListIndex[j]].vectorY = -pos[1];
                    largeNeighborList[j][largeNeighborListIndex[j]].vectorZ = -pos[2];
                    largeNeighborListIndex[j]++;
                    largeNeighborCount++;
                    if(temp_dist2 <= distance_tol2){
                        neighborList[i][neighborListIndex[i]] = j;
                        neighborListIndex[i]++;
                        neighborList[j][neighborListIndex[j]] = i;
                        neighborListIndex[j]++;
                        neighborCount++;
                    }
                }
            }
        }

        // now build our hbond lists
        for (i=0; i<neighborList.size(); i++){
            // load the reference water
            if (fiveAtomWater){
                water15[0] = oPosX[i];
                water15[1] = oPosY[i];
                water15[2] = oPosZ[i];
                water15[3] = hPosX[4*i];
                water15[4] = hPosY[4*i];
                water15[5] = hPosZ[4*i];
                water15[6] = hPosX[4*i+1];
                water15[7] = hPosY[4*i+1];
                water15[8] = hPosZ[4*i+1];
                water15[9] = hPosX[4*i+2];
                water15[10] = hPosY[4*i+2];
                water15[11] = hPosZ[4*i+2];
                water15[12] = hPosX[4*i+3];
                water15[13] = hPosY[4*i+3];
                water15[14] = hPosZ[4*i+3];
            } else {
                water1[0] = oPosX[i];
                water1[1] = oPosY[i];
                water1[2] = oPosZ[i];
                water1[3] = hPosX[2*i];
                water1[4] = hPosY[2*i];
                water1[5] = hPosZ[2*i];
                water1[6] = hPosX[2*i+1];
                water1[7] = hPosY[2*i+1];
                water1[8] = hPosZ[2*i+1];
            }

            for (j=0; j<neighborListIndex[i]; j++){
                if (neighborList[i][j] > i){
                    // load the water
                    if (fiveAtomWater){
                        water25[0] = oPosX[neighborList[i][j]];
                        water25[1] = oPosY[neighborList[i][j]];
                        water25[2] = oPosZ[neighborList[i][j]];
                        water25[3] = hPosX[4*neighborList[i][j]];
                        water25[4] = hPosY[4*neighborList[i][j]];
                        water25[5] = hPosZ[4*neighborList[i][j]];
                        water25[6] = hPosX[4*neighborList[i][j]+1];
                        water25[7] = hPosY[4*neighborList[i][j]+1];
                        water25[8] = hPosZ[4*neighborList[i][j]+1];
                        water25[9] = hPosX[4*neighborList[i][j]+2];
                        water25[10] = hPosY[4*neighborList[i][j]+2];
                        water25[11] = hPosZ[4*neighborList[i][j]+2];
                        water25[12] = hPosX[4*neighborList[i][j]+3];
                        water25[13] = hPosY[4*neighborList[i][j]+3];
                        water25[14] = hPosZ[4*neighborList[i][j]+3];
                    } else {
                        water2[0] = oPosX[neighborList[i][j]];
                        water2[1] = oPosY[neighborList[i][j]];
                        water2[2] = oPosZ[neighborList[i][j]];
                        water2[3] = hPosX[2*neighborList[i][j]];
                        water2[4] = hPosY[2*neighborList[i][j]];
                        water2[5] = hPosZ[2*neighborList[i][j]];
                        water2[6] = hPosX[2*neighborList[i][j]+1];
                        water2[7] = hPosY[2*neighborList[i][j]+1];
                        water2[8] = hPosZ[2*neighborList[i][j]+1];
                    }

                    if(fiveAtomWater){
                        if(calculator->IsWat5HBond(water15, water25, hmat, invhmat, tempHBVec, i, neighborList[i][j])){
                            hbondList[i][hbondListIndex[i]] = neighborList[i][j];
                            hbondVecX[i][hbondListIndex[i]] = 0.5*tempHBVec[0];
                            hbondVecY[i][hbondListIndex[i]] = 0.5*tempHBVec[1];
                            hbondVecZ[i][hbondListIndex[i]] = 0.5*tempHBVec[2];
                            hbondListIndex[i]++;
                            hbondList[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = i;
                            hbondVecX[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[0];
                            hbondVecY[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[1];
                            hbondVecZ[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[2];
                            hbondListIndex[neighborList[i][j]]++;
                            hbondCount++;
                        }
                    } else {
                        if(calculator->IsWat3HBond(water1, water2, hmat, invhmat, tempHBVec, i, neighborList[i][j])){
                            if (directionalityOpt == 1){
                                tempHBVec2[0] = water2[0]-water1[0];
                                tempHBVec2[1] = water2[1]-water1[1];
                                tempHBVec2[2] = water2[2]-water1[2];
                                // do vector wrapping of periodic boundary conditions
                                for (k=0; k<3; k++){
                                    tempScaledHBVec[k] = tempHBVec2[k] * invhmat[k][k];
                                    tempScaledHBVec[k] -= round(tempScaledHBVec[k]);
                                    tempHBVec2[k] = tempScaledHBVec[k] * hmat[k][k];
                                }
                                temp_dist2 = tempHBVec2[0]*tempHBVec2[0] + tempHBVec2[1]*tempHBVec2[1] + tempHBVec2[2]*tempHBVec2[2]; 
                                temp_dist_inv = 1.0 / sqrt(temp_dist2);
                                tempHBVec2[0] *= temp_dist_inv;
                                tempHBVec2[1] *= temp_dist_inv;
                                tempHBVec2[2] *= temp_dist_inv;
                                // check if it is a donating HBond
                                if ((tempHBVec[0]*tempHBVec2[0] + tempHBVec[1]*tempHBVec2[1] + tempHBVec[2]*tempHBVec2[2]) > 0){
                                    is_donor = 1;
                                } else {
                                    is_donor = 0;
                                }

                                if (is_donor == 1){
                                    hbondList[i][hbondListIndex[i]] = neighborList[i][j];
                                    hbondVecX[i][hbondListIndex[i]] = 0.5*tempHBVec[0];
                                    hbondVecY[i][hbondListIndex[i]] = 0.5*tempHBVec[1];
                                    hbondVecZ[i][hbondListIndex[i]] = 0.5*tempHBVec[2];
                                    hbondListIndex[i]++;
                                } else {
                                    hbondList[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = i;
                                    hbondVecX[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[0];
                                    hbondVecY[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[1];
                                    hbondVecZ[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[2];
                                    hbondListIndex[neighborList[i][j]]++;
                                }
                            } else {
                                hbondList[i][hbondListIndex[i]] = neighborList[i][j];
                                hbondVecX[i][hbondListIndex[i]] = 0.5*tempHBVec[0];
                                hbondVecY[i][hbondListIndex[i]] = 0.5*tempHBVec[1];
                                hbondVecZ[i][hbondListIndex[i]] = 0.5*tempHBVec[2];
                                hbondListIndex[i]++;
                                hbondList[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = i;
                                hbondVecX[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[0];
                                hbondVecY[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[1];
                                hbondVecZ[neighborList[i][j]][hbondListIndex[neighborList[i][j]]] = 0.5*tempHBVec[2];
                                hbondListIndex[neighborList[i][j]]++;
                            }
                            hbondCount++;

                        }
                    }
                }
            }
        }

        // build the nearest neighbor list and calculate tetrahedrality
        for (i = 0; i < largeNeighborList.size(); i++){
            // let's sort our nearest neighbor list!
            singleSort = largeNeighborList[i]; 

            // sort based on distances
            sort(singleSort.begin(), singleSort.end(), compareByDistance);

            // check neighbors to see if we need to sort with H-bonding as a
            // priority...
            closest_neighbor = singleSort[0];
            count = 0;
            for (j=1; j<singleSort.size(); j++){
                temp_dist = singleSort[j].distance - closest_neighbor.distance;
                if (temp_dist < neighbor_tol) count++;
            }

            if ( count > 4 ){
                // load the 4 neighbors with a H-bonding+distance preference, 
                // then distance only up to the 4 total neighbors
                tmp_count = 0;
                for (j=0; j<hbondListIndex[i]; j++){
                    for (k=0; k<count; k++){
                        if (singleSort[k].id == hbondList[i][j]){ 
                            nearestNeighborList[i][tmp_count] = singleSort[k];
                            tmp_count++;
                            break;
                        }
                    }
                    if (tmp_count > 3) break;
                }
                if (tmp_count < 4){
                    // load in closest neighbors that haven't already been
                    // selected as H-bonding
                    for (j=0; j<count; j++){
                        loaded = 0; 
                        for (k=0; k<tmp_count; k++){
                            if (nearestNeighborList[i][k].id == singleSort[j].id) {
                                loaded = 1;
                            }
                        }
                        if (loaded == 0){
                            nearestNeighborList[i][tmp_count] = singleSort[j];
                            tmp_count++;
                        }
                        if (tmp_count > 3) break;
                    }
                }
            } else {
                // load the first 4 neighbors in the nearest neighbor lists
                for (j = 0; j < 4; j++){
                    nearestNeighborList[i][j] = singleSort[j];
                }
            }

            // load the sorted list of neighbors back in
            largeNeighborList[i] = singleSort;

            // tetrahedrality double loop over the nearest neighbors of each molecule
            cosine_sum = 0;
            for (j = 0; j < 3; j++){
                inv_nearestNeighborMag1 = 1.0/nearestNeighborList[i][j].distance;
                unitVec1x = nearestNeighborList[i][j].vectorX * inv_nearestNeighborMag1;
                unitVec1y = nearestNeighborList[i][j].vectorY * inv_nearestNeighborMag1;
                unitVec1z = nearestNeighborList[i][j].vectorZ * inv_nearestNeighborMag1;

                for (k = j + 1; k < 4; k++){
                    inv_nearestNeighborMag2 = 1.0/nearestNeighborList[i][k].distance;
                    unitVec2x = nearestNeighborList[i][k].vectorX * inv_nearestNeighborMag2;
                    unitVec2y = nearestNeighborList[i][k].vectorY * inv_nearestNeighborMag2;
                    unitVec2z = nearestNeighborList[i][k].vectorZ * inv_nearestNeighborMag2;

                    dotProduct = (unitVec1x * unitVec2x) + (unitVec1y * unitVec2y) + (unitVec1z * unitVec2z);
                    cosine_sum += pow((dotProduct + 1.0/3.0), 2);
                }
            }
            tetrahedralParam = 1 - (0.375 * cosine_sum);
            tetrahedrality.push_back(tetrahedralParam);
        }

        // calculate the tetrahedrality for the frame

        avg_tetrahedrality = 0;
        tet_square = 0;
        for(i=0; i<tetrahedrality.size(); i++){
            avg_tetrahedrality += tetrahedrality[i];
            tet_square += tetrahedrality[i]*tetrahedrality[i];
        }
        avg_tetrahedrality /= tetrahedrality.size();
        tet_square /= tetrahedrality.size();
        stdev_tetrahedrality = sqrt(tet_square - avg_tetrahedrality*avg_tetrahedrality);


	// output pdb file with tetrahedrality and use it to color the structure
	if (tetraPDBOutOpt){
		outputer_PDB << "COMPND" << "    " << tetraPDBFileName << "\n";
		outputer_PDB << "AUTHOR" << "    " << "PDB file generated by nucleation_tracker program" << "\n";
		outputer_PDB << "CRYST1" << right << setw(9) << fixed << setprecision(3) << boxLength[0] << setw(9) << boxLength[1] << setw(9) << boxLength[2] << setw(7) << fixed << setprecision(2) << 90.00 << setw(7) << 90.00 << setw(7) << 90.00 << " " << "P" << " " << 1 << setw(12) << 1 << "\n";
		m = 1;
		n = 0;
		o = 1;
		occupancy = 1;

		if (2*oPosX.size() == hPosX.size()){
			for (i = 0; i < oPosX.size(); i++){
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "OW" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << oPosX[i] << setw(8) << oPosY[i] << setw(8) << oPosZ[i] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "O" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW1" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n] << setw(8) << hPosY[n] << setw(8) << hPosZ[n] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW2" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n+1] << setw(8) << hPosY[n+1] << setw(8) << hPosZ[n+1] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				n += 2;
				m++;
				o++;
			}
		}

		else {
			for (i = 0; i < oPosX.size(); i++){
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "OW" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << oPosX[i] << setw(8) << oPosY[i] << setw(8) << oPosZ[i] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "O" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW1" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n] << setw(8) << hPosY[n] << setw(8) << hPosZ[n] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW2" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n+1] << setw(8) << hPosY[n+1] << setw(8) << hPosZ[n+1] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW3" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n+2] << setw(8) << hPosY[n+2] << setw(8) << hPosZ[n+2] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				m++;
				outputer_PDB << left << setw(6) << "ATOM" << right << setw(5) << m << setw(5) << "HW4" << setw(4) << "WAT" << setw(2) << "A" << setw(4) << o << fixed << setprecision(3) <<  setw(12) << hPosX[n+3] << setw(8) << hPosY[n+3] << setw(8) << hPosZ[n+3] << fixed << setprecision(2) << setw(6) << occupancy << setw(6) << tetrahedrality[i] << setw(12) << "H" << setw(2) << 0 << "\n";
				n += 4;
				m++;
				o++;
			}
		}
/*
		m = 1;
		if (2*oPosX.size() == hPosX.size()){
			for(j = 0; j < 3*oPosX.size(); j++){
				outputer_PDB << "CONECT" << right << setw(5) << m << "\n";
				m++;
			}
		}
		else{
			for(j = 0; j < 5*oPosX.size(); j++){
				outputer_PDB << "CONECT" << right << setw(5) << m << "\n";
				m++;
			}
		}
*/
		if (2*oPosX.size() == hPosX.size()){
			outputer_PDB << "MASTER" << right << setw(9) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 3*oPosX.size() << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << "\n" << left << setw(6) << "END";
		}
		else{
			outputer_PDB << "MASTER" << right << setw(9) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(5) << 5*oPosX.size() << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << "\n" << left << setw(6) << "END";
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
                                if (maxRingOpt > 3 && index4 != index3 && index4 != index2){
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
                                        if (maxRingOpt > 4 && index5 != index4 && index5 != index3 && index5 != index2){
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
                                                if (maxRingOpt > 5 && index6 != index5 && index6 != index4 && index6 != index3 && index6 != index2){
                                                    for(o=0; o<hbondListIndex[index6]; o++){
                                                        index7 = hbondList[index6][o];
                                                        // test if a 6 member ring
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
                                                        if (maxRingOpt > 6 && index7 != index6 && index7 != index5 && index7 != index4 && index7 != index3 && index7 != index2){
                                                            for(p=0; p<hbondListIndex[index7]; p++){
                                                                index8 = hbondList[index7][p];
                                                                // test if a 7 member ring
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

                                                                if (maxRingOpt > 7 && index8 != index7 && index8 != index6 && index8 != index5 && index8 != index4 && index8 != index3 && index8 != index2){
                                                                    for(q=0; q<hbondListIndex[index8]; q++){
                                                                        index9 = hbondList[index8][q];
                                                                        // test if a 8 member ring
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

                                                                            if (maxRingOpt > 8 && index9 != index8 && index9 != index7 && index9 != index6 && index9 != index5 && index9 != index4 && index9 != index3 && index9 != index2){
                                                                                for(r=0; r<hbondListIndex[index9]; r++){
                                                                                    index10 = hbondList[index9][r];
                                                                                    // test if a 9 member ring
                                                                                    if (index10 == i){
                                                                                        tempVec.push_back(i); 
                                                                                        tempVec.push_back(index2); 
                                                                                        tempVec.push_back(index3);
                                                                                        tempVec.push_back(index4);
                                                                                        tempVec.push_back(index5);
                                                                                        tempVec.push_back(index6);
                                                                                        tempVec.push_back(index7);
                                                                                        tempVec.push_back(index8);
                                                                                        tempVec.push_back(index9);
                                                                                        sort(tempVec.begin(),tempVec.end());
                                                                                        ring9members.push_back(tempVec);
                                                                                        tempVec.clear();
                                                                                        break;
                                                                                    }

                                                                                    if (maxRingOpt > 9 && index10 != index9 && index10 != index8 && index10 != index7 && index10 != index6 && index10 != index5 && index10 != index4 && index10 != index3 && index10 != index2){
                                                                                        for(s=0; s<hbondListIndex[index10]; s++){
                                                                                            index11 = hbondList[index10][s];
                                                                                            // test if a 10 member ring
                                                                                            if (index11 == i){
                                                                                                tempVec.push_back(i); 
                                                                                                tempVec.push_back(index2); 
                                                                                                tempVec.push_back(index3);
                                                                                                tempVec.push_back(index4);
                                                                                                tempVec.push_back(index5);
                                                                                                tempVec.push_back(index6);
                                                                                                tempVec.push_back(index7);
                                                                                                tempVec.push_back(index8);
                                                                                                tempVec.push_back(index9);
                                                                                                tempVec.push_back(index10);
                                                                                                sort(tempVec.begin(),tempVec.end());
                                                                                                ring10members.push_back(tempVec);
                                                                                                tempVec.clear();
                                                                                                break;
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
				    }
				}
			    }
			}
		    }
		}
	    }
//            // purge atom i from hbond lists
//            // ...this is only valid if we use common vertex be default
//            for(j=0; j<hbondListIndex[i]; j++){
//                index2 = hbondList[i][j];
//                for(k=0; k<hbondListIndex[index2]; k++){
//                    hbondList[index2][k] = hbondList[index2][k+1];
//                }
//                hbondListIndex[index2]--;
//            }
        }
        
        // cerr << "initial ring set\n";
        // cerr << ring3members.size() << " is the 3 member ring count\n";
        // cerr << ring4members.size() << " is the 4 member ring count\n";
        // cerr << ring5members.size() << " is the 5 member ring count\n";
        // cerr << ring6members.size() << " is the 6 member ring count\n";
        // cerr << ring7members.size() << " is the 7 member ring count\n";
        // cerr << ring8members.size() << " is the 8 member ring count\n";
        // for(i=0; i<ring3members.size(); i++){
        //     cerr << "ring3: " << ring3members[i][0] << " " << ring3members[i][1] << " " << ring3members[i][2] << "\n";
        // }

        // prune the subset information
        // ...remove duplicate self ring lists first...
        for(i=ring3members.size()-1; i>0; i--){
            for (j=i-1; j>=0; j--){
                if(includes(ring3members[j].begin(),ring3members[j].end(), ring3members[i].begin(), ring3members[i].end())){
                    ring3members.erase(ring3members.begin() + i);
                    break;
                }
            }
        }
        if (maxRingOpt > 3){
            for(i=ring4members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring4members[j].begin(),ring4members[j].end(), ring4members[i].begin(), ring4members[i].end())){
                        ring4members.erase(ring4members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 4){
            for(i=ring5members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring5members[j].begin(),ring5members[j].end(), ring5members[i].begin(), ring5members[i].end())){
                        ring5members.erase(ring5members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 5){
            for(i=ring6members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring6members[j].begin(),ring6members[j].end(), ring6members[i].begin(), ring6members[i].end())){
                        ring6members.erase(ring6members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 6){
            for(i=ring7members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring7members[j].begin(),ring7members[j].end(), ring7members[i].begin(), ring7members[i].end())){
                        ring7members.erase(ring7members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 7){
            for(i=ring8members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring8members[j].begin(),ring8members[j].end(), ring8members[i].begin(), ring8members[i].end())){
                        ring8members.erase(ring8members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 8){
            for(i=ring9members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring9members[j].begin(),ring9members[j].end(), ring9members[i].begin(), ring9members[i].end())){
                        ring9members.erase(ring9members.begin() + i);
                        break;
                    }
                }
            }
        }
        if (maxRingOpt > 9){
            for(i=ring10members.size()-1; i>0; i--){
                for (j=i-1; j>=0; j--){
                    if(includes(ring10members[j].begin(),ring10members[j].end(), ring10members[i].begin(), ring10members[i].end())){
                        ring10members.erase(ring10members.begin() + i);
                        break;
                    }
                }
            }
        }

        // cerr << "\n...after pruning self rings:\n";
        // cerr << ring3members.size() << " is the 3 member ring count\n";
        // cerr << ring4members.size() << " is the 4 member ring count\n";
        // cerr << ring5members.size() << " is the 5 member ring count\n";
        // cerr << ring6members.size() << " is the 6 member ring count\n";
        // cerr << ring7members.size() << " is the 7 member ring count\n";
        // cerr << ring8members.size() << " is the 8 member ring count\n";
        if (ringOutOpt == 0){
            // Sort through the rings to insure only minimal path rings are included for
            // each water molecule.
            //
            // This is the common vertex algorithm.
            //
            for(i=0; i<oPosX.size(); i++){
                bool doEliminateRings = false;
                for (j=ring3members.size()-1; j>=0; j--){
                    for (k=0; k<3; k++){
                        if (ring3members[j][k] == i){
                            doEliminateRings = true;
                            break;
                        }
                    }
                }

                if (maxRingOpt > 3){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring4members.size()-1; j>=0; j--){
                            for (k=0; k<4; k++){
                                if (ring4members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring4members.size()-1; j>=0; j--){
                            for (k=0; k<4; k++){
                                if (ring4members[j][k] == i){
                                    ring4members.erase(ring4members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 4){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring5members.size()-1; j>=0; j--){
                            for (k=0; k<5; k++){
                                if (ring5members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring5members.size()-1; j>=0; j--){
                            for (k=0; k<5; k++){
                                if (ring5members[j][k] == i){
                                    ring5members.erase(ring5members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 5){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring6members.size()-1; j>=0; j--){
                            for (k=0; k<6; k++){
                                if (ring6members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring6members.size()-1; j>=0; j--){
                            for (k=0; k<6; k++){
                                if (ring6members[j][k] == i){
                                    ring6members.erase(ring6members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 6){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring7members.size()-1; j>=0; j--){
                            for (k=0; k<7; k++){
                                if (ring7members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring7members.size()-1; j>=0; j--){
                            for (k=0; k<7; k++){
                                if (ring7members[j][k] == i){
                                    ring7members.erase(ring7members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 7){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring8members.size()-1; j>=0; j--){
                            for (k=0; k<8; k++){
                                if (ring8members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring8members.size()-1; j>=0; j--){
                            for (k=0; k<8; k++){
                                if (ring8members[j][k] == i){
                                    ring8members.erase(ring8members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 8){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring9members.size()-1; j>=0; j--){
                            for (k=0; k<9; k++){
                                if (ring9members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring9members.size()-1; j>=0; j--){
                            for (k=0; k<9; k++){
                                if (ring9members[j][k] == i){
                                    ring9members.erase(ring9members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxRingOpt > 9){
                    if (doEliminateRings == false){
                        // continue our search for a minimum ring size for this water molecule
                        for (j=ring10members.size()-1; j>=0; j--){
                            for (k=0; k<10; k++){
                                if (ring10members[j][k] == i){
                                    doEliminateRings = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        // eliminate all rings containing this water molecule
                        for (j=ring10members.size()-1; j>=0; j--){
                            for (k=0; k<10; k++){
                                if (ring10members[j][k] == i){
                                    ring10members.erase(ring10members.begin() + j);
                                    break;
        			}
                            }
                        }
                    }
                }
            }
        } else if(ringOutOpt == 1){
           // common bond ring closure/pruning
           // run over possible angle sequences by making a temporary 2 member
           // sequences to probe all remaining larger rings.
            for (i=0; i<ring3members.size(); i++){
                for (k=0; k<ring3members[i].size()-2; k++){
                    bondSequence.push_back(ring3members[i][k]);
                    bondSequence.push_back(ring3members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 3){
                        for (j=ring4members.size()-1; j>=0; j--){
                            if(includes(ring4members[j].begin(),ring4members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring4members.erase(ring4members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 4){
                            for (j=ring5members.size()-1; j>=0; j--){
                                if(includes(ring5members[j].begin(),ring5members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring5members.erase(ring5members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 5){
                                for (j=ring6members.size()-1; j>=0; j--){
                                    if(includes(ring6members[j].begin(),ring6members[j].end(), bondSequence.begin(), bondSequence.end())){
                                        ring6members.erase(ring6members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 6){
                                    for (j=ring7members.size()-1; j>=0; j--){
                                        if(includes(ring7members[j].begin(),ring7members[j].end(), bondSequence.begin(), bondSequence.end())){
                                            ring7members.erase(ring7members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 7){
                                        for (j=ring8members.size()-1; j>=0; j--){
                                            if(includes(ring8members[j].begin(),ring8members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                ring8members.erase(ring8members.begin() + j);
                                                break;
                                            }
                                        }
                                        if (maxRingOpt > 8){
                                            for (j=ring9members.size()-1; j>=0; j--){
                                                if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                    ring9members.erase(ring9members.begin() + j);
                                                    break;
                                                }
                                            }
                                            if (maxRingOpt > 9){
                                                for (j=ring10members.size()-1; j>=0; j--){
                                                    if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                        ring10members.erase(ring10members.begin() + j);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring4members.size(); i++){
                for (k=0; k<ring4members[i].size()-2; k++){
                    bondSequence.push_back(ring4members[i][k]);
                    bondSequence.push_back(ring4members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 4){
                        for (j=ring5members.size()-1; j>=0; j--){
                            if(includes(ring5members[j].begin(),ring5members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring5members.erase(ring5members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 5){
                            for (j=ring6members.size()-1; j>=0; j--){
                                if(includes(ring6members[j].begin(),ring6members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring6members.erase(ring6members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 6){
                                for (j=ring7members.size()-1; j>=0; j--){
                                    if(includes(ring7members[j].begin(),ring7members[j].end(), bondSequence.begin(), bondSequence.end())){
                                        ring7members.erase(ring7members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 7){
                                    for (j=ring8members.size()-1; j>=0; j--){
                                        if(includes(ring8members[j].begin(),ring8members[j].end(), bondSequence.begin(), bondSequence.end())){
                                            ring8members.erase(ring8members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 8){
                                        for (j=ring9members.size()-1; j>=0; j--){
                                            if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                ring9members.erase(ring9members.begin() + j);
                                                break;
                                            }
                                        }
                                        if (maxRingOpt > 9){
                                            for (j=ring10members.size()-1; j>=0; j--){
                                                if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                    ring10members.erase(ring10members.begin() + j);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring5members.size(); i++){
                for (k=0; k<ring5members[i].size()-2; k++){
                    bondSequence.push_back(ring5members[i][k]);
                    bondSequence.push_back(ring5members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 5){
                        for (j=ring6members.size()-1; j>=0; j--){
                            if(includes(ring6members[j].begin(),ring6members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring6members.erase(ring6members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 6){
                            for (j=ring7members.size()-1; j>=0; j--){
                                if(includes(ring7members[j].begin(),ring7members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring7members.erase(ring7members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 7){
                                for (j=ring8members.size()-1; j>=0; j--){
                                    if(includes(ring8members[j].begin(),ring8members[j].end(), bondSequence.begin(), bondSequence.end())){
                                        ring8members.erase(ring8members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 8){
                                    for (j=ring9members.size()-1; j>=0; j--){
                                        if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                            ring9members.erase(ring9members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 9){
                                        for (j=ring10members.size()-1; j>=0; j--){
                                            if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                                ring10members.erase(ring10members.begin() + j);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring6members.size(); i++){
                for (k=0; k<ring6members[i].size()-2; k++){
                    bondSequence.push_back(ring6members[i][k]);
                    bondSequence.push_back(ring6members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 6){
                        for (j=ring7members.size()-1; j>=0; j--){
                            if(includes(ring7members[j].begin(),ring7members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring7members.erase(ring7members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 7){
                            for (j=ring8members.size()-1; j>=0; j--){
                                if(includes(ring8members[j].begin(),ring8members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring8members.erase(ring8members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 8){
                                for (j=ring9members.size()-1; j>=0; j--){
                                    if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                        ring9members.erase(ring9members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 9){
                                    for (j=ring10members.size()-1; j>=0; j--){
                                        if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                            ring10members.erase(ring10members.begin() + j);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring7members.size(); i++){
                for (k=0; k<ring7members[i].size()-2; k++){
                    bondSequence.push_back(ring7members[i][k]);
                    bondSequence.push_back(ring7members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 7){
                        for (j=ring8members.size()-1; j>=0; j--){
                            if(includes(ring8members[j].begin(),ring8members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring8members.erase(ring8members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 8){
                            for (j=ring9members.size()-1; j>=0; j--){
                                if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring9members.erase(ring9members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 9){
                                for (j=ring10members.size()-1; j>=0; j--){
                                    if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                        ring10members.erase(ring10members.begin() + j);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring8members.size(); i++){
                for (k=0; k<ring8members[i].size()-2; k++){
                    bondSequence.push_back(ring8members[i][k]);
                    bondSequence.push_back(ring8members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 8){
                        for (j=ring9members.size()-1; j>=0; j--){
                            if(includes(ring9members[j].begin(),ring9members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring9members.erase(ring9members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 9){
                            for (j=ring10members.size()-1; j>=0; j--){
                                if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                    ring10members.erase(ring10members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
            for (i=0; i<ring9members.size(); i++){
                for (k=0; k<ring9members[i].size()-2; k++){
                    bondSequence.push_back(ring9members[i][k]);
                    bondSequence.push_back(ring9members[i][k+1]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 9){
                        for (j=ring10members.size()-1; j>=0; j--){
                            if(includes(ring10members[j].begin(),ring10members[j].end(), bondSequence.begin(), bondSequence.end())){
                                ring10members.erase(ring10members.begin() + j);
                                break;
                            }
                        }
                    }
                    bondSequence.clear();
                }
            }
        } else if(ringOutOpt == 2){
            // common angle ring closure/pruning
            //
            // start by removing any ring that contains one of the 3 member rings these are common angle cases
            for (i=0; i<ring3members.size(); i++){
                angleSequence = ring3members[i];
                if (maxRingOpt > 3){
                    // ...3s in 4s...
                    for (j=ring4members.size()-1; j>=0; j--){
                        if(includes(ring4members[j].begin(),ring4members[j].end(), angleSequence.begin(), angleSequence.end())){
                            ring4members.erase(ring4members.begin() + j);
                            break;
                        }
                    }
                    if (maxRingOpt > 4){
                        // ...3s in 5s...
                        for (j=ring5members.size()-1; j>=0; j--){
                            if(includes(ring5members[j].begin(),ring5members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring5members.erase(ring5members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 5){
                            // ...3s in 6s...
                            for (j=ring6members.size()-1; j>=0; j--){
                                if(includes(ring6members[j].begin(),ring6members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring6members.erase(ring6members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 6){
                                // ...3s in 7s...
                                for (j=ring7members.size()-1; j>=0; j--){
                                    if(includes(ring7members[j].begin(),ring7members[j].end(), angleSequence.begin(), angleSequence.end())){
                                        ring7members.erase(ring7members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 7){
                                    // ...3s in 8s...
                                    for (j=ring8members.size()-1; j>=0; j--){
                                        if(includes(ring8members[j].begin(),ring8members[j].end(), angleSequence.begin(), angleSequence.end())){
                                            ring8members.erase(ring8members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 8){
                                        // ...3s in 9s...
                                        for (j=ring9members.size()-1; j>=0; j--){
                                            if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                                ring9members.erase(ring9members.begin() + j);
                                                break;
                                            }
                                        }
                                        if (maxRingOpt > 9){
                                            // ...3s in 10s...
                                            for (j=ring10members.size()-1; j>=0; j--){
                                                if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                                    ring10members.erase(ring10members.begin() + j);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                angleSequence.clear();
            }
           // run over possible angle sequences by making a temporary 3 member
           // sequences to probe all remaining larger rings.
            for (i=0; i<ring4members.size(); i++){
                for (k=0; k<ring4members[i].size()-3; k++){
                    angleSequence.push_back(ring4members[i][k]);
                    angleSequence.push_back(ring4members[i][k+1]);
                    angleSequence.push_back(ring4members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 4){
                        for (j=ring5members.size()-1; j>=0; j--){
                            if(includes(ring5members[j].begin(),ring5members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring5members.erase(ring5members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 5){
                            for (j=ring6members.size()-1; j>=0; j--){
                                if(includes(ring6members[j].begin(),ring6members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring6members.erase(ring6members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 6){
                                for (j=ring7members.size()-1; j>=0; j--){
                                    if(includes(ring7members[j].begin(),ring7members[j].end(), angleSequence.begin(), angleSequence.end())){
                                        ring7members.erase(ring7members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 7){
                                    for (j=ring8members.size()-1; j>=0; j--){
                                        if(includes(ring8members[j].begin(),ring8members[j].end(), angleSequence.begin(), angleSequence.end())){
                                            ring8members.erase(ring8members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 8){
                                        for (j=ring9members.size()-1; j>=0; j--){
                                            if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                                ring9members.erase(ring9members.begin() + j);
                                                break;
                                            }
                                        }
                                        if (maxRingOpt > 9){
                                            for (j=ring10members.size()-1; j>=0; j--){
                                                if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                                    ring10members.erase(ring10members.begin() + j);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
            for (i=0; i<ring5members.size(); i++){
                for (k=0; k<ring5members[i].size()-3; k++){
                    angleSequence.push_back(ring5members[i][k]);
                    angleSequence.push_back(ring5members[i][k+1]);
                    angleSequence.push_back(ring5members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 5){
                        for (j=ring6members.size()-1; j>=0; j--){
                            if(includes(ring6members[j].begin(),ring6members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring6members.erase(ring6members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 6){
                            for (j=ring7members.size()-1; j>=0; j--){
                                if(includes(ring7members[j].begin(),ring7members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring7members.erase(ring7members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 7){
                                for (j=ring8members.size()-1; j>=0; j--){
                                    if(includes(ring8members[j].begin(),ring8members[j].end(), angleSequence.begin(), angleSequence.end())){
                                        ring8members.erase(ring8members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 8){
                                    for (j=ring9members.size()-1; j>=0; j--){
                                        if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                            ring9members.erase(ring9members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 9){
                                        for (j=ring10members.size()-1; j>=0; j--){
                                            if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                                ring10members.erase(ring10members.begin() + j);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
            for (i=0; i<ring6members.size(); i++){
                for (k=0; k<ring6members[i].size()-3; k++){
                    angleSequence.push_back(ring6members[i][k]);
                    angleSequence.push_back(ring6members[i][k+1]);
                    angleSequence.push_back(ring6members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 6){
                        for (j=ring7members.size()-1; j>=0; j--){
                            if(includes(ring7members[j].begin(),ring7members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring7members.erase(ring7members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 7){
                            for (j=ring8members.size()-1; j>=0; j--){
                                if(includes(ring8members[j].begin(),ring8members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring8members.erase(ring8members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 8){
                                for (j=ring9members.size()-1; j>=0; j--){
                                    if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                        ring9members.erase(ring9members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 9){
                                    for (j=ring10members.size()-1; j>=0; j--){
                                        if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                            ring10members.erase(ring10members.begin() + j);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
            for (i=0; i<ring7members.size(); i++){
                for (k=0; k<ring7members[i].size()-3; k++){
                    angleSequence.push_back(ring7members[i][k]);
                    angleSequence.push_back(ring7members[i][k+1]);
                    angleSequence.push_back(ring7members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 7){
                        for (j=ring8members.size()-1; j>=0; j--){
                            if(includes(ring8members[j].begin(),ring8members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring8members.erase(ring8members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 8){
                            for (j=ring9members.size()-1; j>=0; j--){
                                if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring9members.erase(ring9members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 9){
                                for (j=ring10members.size()-1; j>=0; j--){
                                    if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                        ring10members.erase(ring10members.begin() + j);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
            for (i=0; i<ring8members.size(); i++){
                for (k=0; k<ring8members[i].size()-3; k++){
                    angleSequence.push_back(ring8members[i][k]);
                    angleSequence.push_back(ring8members[i][k+1]);
                    angleSequence.push_back(ring8members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 8){
                        for (j=ring9members.size()-1; j>=0; j--){
                            if(includes(ring9members[j].begin(),ring9members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring9members.erase(ring9members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 9){
                            for (j=ring10members.size()-1; j>=0; j--){
                                if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                    ring10members.erase(ring10members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
            for (i=0; i<ring9members.size(); i++){
                for (k=0; k<ring9members[i].size()-3; k++){
                    angleSequence.push_back(ring9members[i][k]);
                    angleSequence.push_back(ring9members[i][k+1]);
                    angleSequence.push_back(ring9members[i][k+2]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 9){
                        for (j=ring10members.size()-1; j>=0; j--){
                            if(includes(ring10members[j].begin(),ring10members[j].end(), angleSequence.begin(), angleSequence.end())){
                                ring10members.erase(ring10members.begin() + j);
                                break;
                            }
                        }
                    }
                    angleSequence.clear();
                }
            }
        } else if(ringOutOpt == 3){
            // common dihedral ring closure/pruning
           // run over possible dihedral sequences by making a temporary 4 member
           // sequences to probe all remaining larger rings.
            for (i=0; i<ring4members.size(); i++){
                for (k=0; k<ring4members[i].size()-4; k++){
                    dihedralSequence.push_back(ring4members[i][k]);
                    dihedralSequence.push_back(ring4members[i][k+1]);
                    dihedralSequence.push_back(ring4members[i][k+2]);
                    dihedralSequence.push_back(ring4members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 4){
                        for (j=ring5members.size()-1; j>=0; j--){
                            if(includes(ring5members[j].begin(),ring5members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring5members.erase(ring5members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 5){
                            for (j=ring6members.size()-1; j>=0; j--){
                                if(includes(ring6members[j].begin(),ring6members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                    ring6members.erase(ring6members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 6){
                                for (j=ring7members.size()-1; j>=0; j--){
                                    if(includes(ring7members[j].begin(),ring7members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                        ring7members.erase(ring7members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 7){
                                    for (j=ring8members.size()-1; j>=0; j--){
                                        if(includes(ring8members[j].begin(),ring8members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                            ring8members.erase(ring8members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 8){
                                        for (j=ring9members.size()-1; j>=0; j--){
                                            if(includes(ring9members[j].begin(),ring9members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                                ring9members.erase(ring9members.begin() + j);
                                                break;
                                            }
                                        }
                                        if (maxRingOpt > 9){
                                            for (j=ring10members.size()-1; j>=0; j--){
                                                if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                                    ring10members.erase(ring10members.begin() + j);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
            for (i=0; i<ring5members.size(); i++){
                for (k=0; k<ring5members[i].size()-4; k++){
                    dihedralSequence.push_back(ring5members[i][k]);
                    dihedralSequence.push_back(ring5members[i][k+1]);
                    dihedralSequence.push_back(ring5members[i][k+2]);
                    dihedralSequence.push_back(ring5members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 5){
                        for (j=ring6members.size()-1; j>=0; j--){
                            if(includes(ring6members[j].begin(),ring6members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring6members.erase(ring6members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 6){
                            for (j=ring7members.size()-1; j>=0; j--){
                                if(includes(ring7members[j].begin(),ring7members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                    ring7members.erase(ring7members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 7){
                                for (j=ring8members.size()-1; j>=0; j--){
                                    if(includes(ring8members[j].begin(),ring8members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                        ring8members.erase(ring8members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 8){
                                    for (j=ring9members.size()-1; j>=0; j--){
                                        if(includes(ring9members[j].begin(),ring9members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                            ring9members.erase(ring9members.begin() + j);
                                            break;
                                        }
                                    }
                                    if (maxRingOpt > 9){
                                        for (j=ring10members.size()-1; j>=0; j--){
                                            if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                                ring10members.erase(ring10members.begin() + j);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
            for (i=0; i<ring6members.size(); i++){
                for (k=0; k<ring6members[i].size()-4; k++){
                    dihedralSequence.push_back(ring6members[i][k]);
                    dihedralSequence.push_back(ring6members[i][k+1]);
                    dihedralSequence.push_back(ring6members[i][k+2]);
                    dihedralSequence.push_back(ring6members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 6){
                        for (j=ring7members.size()-1; j>=0; j--){
                            if(includes(ring7members[j].begin(),ring7members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring7members.erase(ring7members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 7){
                            for (j=ring8members.size()-1; j>=0; j--){
                                if(includes(ring8members[j].begin(),ring8members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                    ring8members.erase(ring8members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 8){
                                for (j=ring9members.size()-1; j>=0; j--){
                                    if(includes(ring9members[j].begin(),ring9members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                        ring9members.erase(ring9members.begin() + j);
                                        break;
                                    }
                                }
                                if (maxRingOpt > 9){
                                    for (j=ring10members.size()-1; j>=0; j--){
                                        if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                            ring10members.erase(ring10members.begin() + j);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
            for (i=0; i<ring7members.size(); i++){
                for (k=0; k<ring7members[i].size()-4; k++){
                    dihedralSequence.push_back(ring7members[i][k]);
                    dihedralSequence.push_back(ring7members[i][k+1]);
                    dihedralSequence.push_back(ring7members[i][k+2]);
                    dihedralSequence.push_back(ring7members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 7){
                        for (j=ring8members.size()-1; j>=0; j--){
                            if(includes(ring8members[j].begin(),ring8members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring8members.erase(ring8members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 8){
                            for (j=ring9members.size()-1; j>=0; j--){
                                if(includes(ring9members[j].begin(),ring9members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                    ring9members.erase(ring9members.begin() + j);
                                    break;
                                }
                            }
                            if (maxRingOpt > 9){
                                for (j=ring10members.size()-1; j>=0; j--){
                                    if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                        ring10members.erase(ring10members.begin() + j);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
            for (i=0; i<ring8members.size(); i++){
                for (k=0; k<ring8members[i].size()-4; k++){
                    dihedralSequence.push_back(ring8members[i][k]);
                    dihedralSequence.push_back(ring8members[i][k+1]);
                    dihedralSequence.push_back(ring8members[i][k+2]);
                    dihedralSequence.push_back(ring8members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 8){
                        for (j=ring9members.size()-1; j>=0; j--){
                            if(includes(ring9members[j].begin(),ring9members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring9members.erase(ring9members.begin() + j);
                                break;
                            }
                        }
                        if (maxRingOpt > 9){
                            for (j=ring10members.size()-1; j>=0; j--){
                                if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                    ring10members.erase(ring10members.begin() + j);
                                    break;
                                }
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
            for (i=0; i<ring9members.size(); i++){
                for (k=0; k<ring9members[i].size()-4; k++){
                    dihedralSequence.push_back(ring9members[i][k]);
                    dihedralSequence.push_back(ring9members[i][k+1]);
                    dihedralSequence.push_back(ring9members[i][k+2]);
                    dihedralSequence.push_back(ring9members[i][k+3]);

                    // now loop to eliminate from larger rings
                    if (maxRingOpt > 9){
                        for (j=ring10members.size()-1; j>=0; j--){
                            if(includes(ring10members[j].begin(),ring10members[j].end(), dihedralSequence.begin(), dihedralSequence.end())){
                                ring10members.erase(ring10members.begin() + j);
                                break;
                            }
                        }
                    }
                    dihedralSequence.clear();
                }
            }
        } else {
            // all rings
        }
        //             cerr << "\n...after elimination:\n";
        // cerr << ring3members.size() << " is the 3 member ring count\n";
        // cerr << ring4members.size() << " is the 4 member ring count\n";
        // cerr << ring5members.size() << " is the 5 member ring count\n";
        // cerr << ring6members.size() << " is the 6 member ring count\n";
        // cerr << ring7members.size() << " is the 7 member ring count\n";
        // cerr << ring8members.size() << " is the 8 member ring count\n";

        // check for "rings" crossing periodic boundary images
         

        // output ring data
        lastRing3 = ring3members.size()-lastRing3;
        lastRing4 = ring4members.size()-lastRing4;
        lastRing5 = ring5members.size()-lastRing5;
        lastRing6 = ring6members.size()-lastRing6;
        lastRing7 = ring7members.size()-lastRing7;
        lastRing8 = ring8members.size()-lastRing8;
        lastRing9 = ring9members.size()-lastRing9;
        lastRing10 = ring10members.size()-lastRing10;
        outputer << setw(8) << frameCount ;
        outputer << setw(8) << hbondCount ;
        outputer << setw(8) << setprecision(4) << avg_tetrahedrality;
        outputer << setw(8) << setprecision(1) << stdev_tetrahedrality;
        outputer << setw(8) << lastRing3 ;
        totalRings = lastRing3;
        if (maxRingOpt > 3){
            outputer << setw(8) << lastRing4;
            totalRings += lastRing4;
            if (maxRingOpt > 4){
                outputer << setw(8) << lastRing5;
                totalRings += lastRing5;
                if (maxRingOpt > 5){
                    outputer << setw(8) << lastRing6;
                    totalRings += lastRing6;
                    if (maxRingOpt > 6){
                        outputer << setw(8) << lastRing7;
                        totalRings += lastRing7;
                        if (maxRingOpt > 7){
                            outputer << setw(8) << lastRing8; 
                            totalRings += lastRing8;
                            if (maxRingOpt > 8){
                            	outputer << setw(9) << lastRing9; 
                            	totalRings += lastRing9;
                        	if (maxRingOpt > 9){
                            	    outputer << setw(10) << lastRing10; 
                            	    totalRings += lastRing10;
                        }
                    }
                }
            }
        }
    }
}
        outputer << "\n";

        // establish our ring probabilities
        if (totalRings == 0){
            cerr << "Error: No closed ring paths found!\n\n";
            return 0;
        }
        ringProbabilities[0] = 0;
        ringProbabilities[1] = 0;
        ringProbabilities[2] = 0;
        ringProbabilities[3] = (float)lastRing3/totalRings;
        ringProbabilities[4] = (float)lastRing4/totalRings;
        ringProbabilities[5] = (float)lastRing5/totalRings;
        ringProbabilities[6] = (float)lastRing6/totalRings;
        ringProbabilities[7] = (float)lastRing7/totalRings;
        ringProbabilities[8] = (float)lastRing8/totalRings;
        ringProbabilities[9] = (float)lastRing9/totalRings;
        ringProbabilities[10] = (float)lastRing10/totalRings;

        if (ringTrajOutOpt){
            outputer_xyz << totalRings << "\n";
            outputer_xyz << 0 << "; " << fixed<<setprecision(5)<<setw(8) << boxLength[0] << " " << hmat[0][1] << " " << hmat[0][2] << "; " << hmat[1][0] << " " << boxLength[1] << " " << hmat[1][2] << "; " << hmat[2][0] << " " << hmat[2][1] << " " << boxLength[2] << "\n";
        }

        // outputer << setw(8) << lastRing7 ;
        // outputer << setw(8) << lastRing8 << "\n";

        if (povrayOutOpt){
            // open a pov outputer
            frameInt.str("");
            frameInt << frameCount;
            frameCountString = frameInt.str();
            povName = "pov_files/" + povName2 + "_" + frameCountString + ".pov";
            povDistName = "pov_files/" + povName2 + "_dist_" + frameCountString + ".pov";
            pov_out.open(povName.c_str());
            povDistOut.open(povDistName.c_str());

            // output ring info. in pov files
            PovObjects *povrayObjects = new PovObjects();

            // print out the header and object info. to the .pov file
            povrayObjects->printHeader(pov_out, boxLength);
            povrayObjects->printRing3(pov_out, boxLength);
            povrayObjects->printRing4(pov_out, boxLength);
            povrayObjects->printRing5(pov_out, boxLength);
            povrayObjects->printRing6(pov_out, boxLength);
            povrayObjects->printRing7(pov_out, boxLength);
            povrayObjects->printRing8(pov_out, boxLength);
            povrayObjects->printRing9(pov_out, boxLength);
            povrayObjects->printRing10(pov_out, boxLength);

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
            povDistOut << "bar(2.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[8] << ", 0.7, 0.5, 0.5, 1, 0.0)\n";
            povDistOut << "bar(-3.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[9] << ", 0.7, 1.0, 0.5, 0.5, 0.0)\n";
            povDistOut << "bar(3.0, " <<setprecision(8)<<setw(10)<< 4.0*ringProbabilities[10] << ", 0.7, 1.0, 0.75, 0.5, 0.0)\n";
        }

        if ((povrayOutOpt || ringTrajOutOpt) && maxRingOpt > 2){
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

                //outputer_xyz << setw(5) << "Li " << setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
                outputer_xyz <<left<< setw(5) << "Li" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
            }

                if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                    red_val = "1.0";
                    green_val = "0.4";
                    blue_val = "1.0";
                    zVal = -0.05;

                    pov_out << "ring3(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";

                    if (xVal < (-0.5*boxLength[0]+buffer_size)){
                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                            pov_out << "ring3(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                            pov_out << "ring3(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else {
                            pov_out << "ring3(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        }
                    } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                            pov_out << "ring3(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                            pov_out << "ring3(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring3(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else {
                            pov_out << "ring3(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        }
                    } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring3(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring3(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                }

        }

        if ((povrayOutOpt || ringTrajOutOpt) && maxRingOpt > 3) {
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
                //outputer_xyz << setw(5) << "Be" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
                if (ringTrajOutOpt){
                    outputer_xyz <<left<< setw(5) << "Be" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                }


                if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                    red_val = "0.5";
                    green_val = "0.0";
                    blue_val = "1.0";
                    zVal = -0.15;

                    pov_out << "ring4(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";

                    if (xVal < (-0.5*boxLength[0]+buffer_size)){
                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                            pov_out << "ring4(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                            pov_out << "ring4(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else {
                            pov_out << "ring4(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        }
                    } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                            pov_out << "ring4(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                            pov_out << "ring4(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                            pov_out << "ring4(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        } else {
                            pov_out << "ring4(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                        }
                    } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                        pov_out << "ring4(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                        pov_out << "ring4(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                    }
                }
            }


            if (maxRingOpt > 4){
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
                    if (ringTrajOutOpt){
                        outputer_xyz <<left<< setw(5) << "B" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                    }


                    if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                        red_val = "0.0";
                        green_val = "0.5";
                        blue_val = "1.0";
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

                if (maxRingOpt > 5){
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
                        if (ringTrajOutOpt){
                            outputer_xyz <<left<< setw(5) << "C" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                        }

                        if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                            red_val = "1.0";
                            green_val = "0.0";
                            blue_val = "0.0";
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

                    if (maxRingOpt > 6){
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
                            if (ringTrajOutOpt){
                                outputer_xyz <<left<< setw(5) << "N" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                            }

                            if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                                red_val = "1.0";
                                green_val = "0.5";
                                blue_val = "0.0";
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

                        if (maxRingOpt > 7){
                            for (i=0; i<ring8members.size(); i++){
                                tempPosX[0]=oPosX[ring8members[i][0]];
                                diffVal = oPosX[ring8members[i][1]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring8members[i][1]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring8members[i][1]]+boxLength[0];
                                } else {
                                    tempPosX[1] = oPosX[ring8members[i][1]];
                                }
                                diffVal = oPosX[ring8members[i][2]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring8members[i][2]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring8members[i][2]]+boxLength[0];
                                } else {
                                    tempPosX[2] = oPosX[ring8members[i][2]];
                                }
                                diffVal = oPosX[ring8members[i][3]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring8members[i][3]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring8members[i][3]]+boxLength[0];
                                } else {
                                    tempPosX[3] = oPosX[ring8members[i][3]];
                                }
                                diffVal = oPosX[ring8members[i][4]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring8members[i][4]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring8members[i][4]]+boxLength[0];
                                } else {
                                    tempPosX[4] = oPosX[ring8members[i][4]];
                                }
                                diffVal = oPosX[ring8members[i][5]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring8members[i][5]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring8members[i][5]]+boxLength[0];
                                } else {
                                    tempPosX[5] = oPosX[ring8members[i][5]];
                                }
                                diffVal = oPosX[ring8members[i][6]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring8members[i][6]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring8members[i][6]]+boxLength[0];
                                } else {
                                    tempPosX[6] = oPosX[ring8members[i][6]];
                                }
                                diffVal = oPosX[ring8members[i][7]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring8members[i][7]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring8members[i][7]]+boxLength[0];
                                } else {
                                    tempPosX[7] = oPosX[ring8members[i][7]];
                                }
                                xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4]+tempPosX[5]+tempPosX[6]+tempPosX[7];
                                xVal /= 8;

                                tempPosY[0]=oPosY[ring8members[i][0]];
                                diffVal = oPosY[ring8members[i][1]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring8members[i][1]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring8members[i][1]]+boxLength[1];
                                } else {
                                    tempPosY[1] = oPosY[ring8members[i][1]];
                                }
                                diffVal = oPosY[ring8members[i][2]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring8members[i][2]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring8members[i][2]]+boxLength[1];
                                } else {
                                    tempPosY[2] = oPosY[ring8members[i][2]];
                                }
                                diffVal = oPosY[ring8members[i][3]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring8members[i][3]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring8members[i][3]]+boxLength[1];
                                } else {
                                    tempPosY[3] = oPosY[ring8members[i][3]];
                                }
                                diffVal = oPosY[ring8members[i][4]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring8members[i][4]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring8members[i][4]]+boxLength[1];
                                } else {
                                    tempPosY[4] = oPosY[ring8members[i][4]];
                                }
                                diffVal = oPosY[ring8members[i][5]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring8members[i][5]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring8members[i][5]]+boxLength[1];
                                } else {
                                    tempPosY[5] = oPosY[ring8members[i][5]];
                                }
                                diffVal = oPosY[ring8members[i][6]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring8members[i][6]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring8members[i][6]]+boxLength[1];
                                } else {
                                    tempPosY[6] = oPosY[ring8members[i][6]];
                                }
                                diffVal = oPosY[ring8members[i][7]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring8members[i][7]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring8members[i][7]]+boxLength[1];
                                } else {
                                    tempPosY[7] = oPosY[ring8members[i][7]];
                                }
                                yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4]+tempPosY[5]+tempPosY[6]+tempPosY[7];
                                yVal /= 8;

                                tempPosZ[0]=oPosZ[ring8members[i][0]];
                                diffVal = oPosZ[ring8members[i][1]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring8members[i][1]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring8members[i][1]]+boxLength[2];
                                } else {
                                    tempPosZ[1] = oPosZ[ring8members[i][1]];
                                }
                                diffVal = oPosZ[ring8members[i][2]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring8members[i][2]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring8members[i][2]]+boxLength[2];
                                } else {
                                    tempPosZ[2] = oPosZ[ring8members[i][2]];
                                }
                                diffVal = oPosZ[ring8members[i][3]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring8members[i][3]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring8members[i][3]]+boxLength[2];
                                } else {
                                    tempPosZ[3] = oPosZ[ring8members[i][3]];
                                }
                                diffVal = oPosZ[ring8members[i][4]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring8members[i][4]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring8members[i][4]]+boxLength[2];
                                } else {
                                    tempPosZ[4] = oPosZ[ring8members[i][4]];
                                }
                                diffVal = oPosZ[ring8members[i][5]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring8members[i][5]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring8members[i][5]]+boxLength[2];
                                } else {
                                    tempPosZ[5] = oPosZ[ring8members[i][5]];
                                }
                                diffVal = oPosZ[ring8members[i][6]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring8members[i][6]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring8members[i][6]]+boxLength[2];
                                } else {
                                    tempPosZ[6] = oPosZ[ring8members[i][6]];
                                }
                                diffVal = oPosZ[ring8members[i][7]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring8members[i][7]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring8members[i][7]]+boxLength[2];
                                } else {
                                    tempPosZ[7] = oPosZ[ring8members[i][7]];
                                }
                                zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4]+tempPosZ[5]+tempPosZ[6]+tempPosZ[7];
                                zVal /= 8;
                                //outputer_xyz << setw(5) << "O" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
                                if (ringTrajOutOpt){
                                    outputer_xyz <<left<< setw(5) << "O" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                                }

                                if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                                    red_val = "0.0";
                                    green_val = "1.0";
                                    blue_val = "0.5";
                                    zVal = -0.25;
                                    pov_out << "ring8(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    if (xVal < (-0.5*boxLength[0]+buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring8(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring8(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring8(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring8(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring8(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring8(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring8(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                        pov_out << "ring8(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                        pov_out << "ring8(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    }
                                }
                            }
                        if (maxRingOpt > 8){
                            for (i=0; i<ring9members.size(); i++){
                                tempPosX[0]=oPosX[ring9members[i][0]];
                                diffVal = oPosX[ring9members[i][1]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring9members[i][1]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring9members[i][1]]+boxLength[0];
                                } else {
                                    tempPosX[1] = oPosX[ring9members[i][1]];
                                }
                                diffVal = oPosX[ring9members[i][2]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring9members[i][2]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring9members[i][2]]+boxLength[0];
                                } else {
                                    tempPosX[2] = oPosX[ring9members[i][2]];
                                }
                                diffVal = oPosX[ring9members[i][3]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring9members[i][3]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring9members[i][3]]+boxLength[0];
                                } else {
                                    tempPosX[3] = oPosX[ring9members[i][3]];
                                }
                                diffVal = oPosX[ring9members[i][4]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring9members[i][4]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring9members[i][4]]+boxLength[0];
                                } else {
                                    tempPosX[4] = oPosX[ring9members[i][4]];
                                }
                                diffVal = oPosX[ring9members[i][5]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring9members[i][5]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring9members[i][5]]+boxLength[0];
                                } else {
                                    tempPosX[5] = oPosX[ring9members[i][5]];
                                }
                                diffVal = oPosX[ring9members[i][6]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring9members[i][6]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring9members[i][6]]+boxLength[0];
                                } else {
                                    tempPosX[6] = oPosX[ring9members[i][6]];
                                }
                                diffVal = oPosX[ring9members[i][7]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring9members[i][7]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring9members[i][7]]+boxLength[0];
                                } else {
                                    tempPosX[7] = oPosX[ring9members[i][7]];
                                }
                                diffVal = oPosX[ring9members[i][8]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[8] = oPosX[ring9members[i][8]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[8] = oPosX[ring9members[i][8]]+boxLength[0];
                                } else {
                                    tempPosX[8] = oPosX[ring9members[i][8]];
                                }
                                xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4]+tempPosX[5]+tempPosX[6]+tempPosX[7]+tempPosX[8];
                                xVal /= 9;

                                tempPosY[0]=oPosY[ring9members[i][0]];
                                diffVal = oPosY[ring9members[i][1]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring9members[i][1]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring9members[i][1]]+boxLength[1];
                                } else {
                                    tempPosY[1] = oPosY[ring9members[i][1]];
                                }
                                diffVal = oPosY[ring9members[i][2]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring9members[i][2]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring9members[i][2]]+boxLength[1];
                                } else {
                                    tempPosY[2] = oPosY[ring9members[i][2]];
                                }
                                diffVal = oPosY[ring9members[i][3]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring9members[i][3]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring9members[i][3]]+boxLength[1];
                                } else {
                                    tempPosY[3] = oPosY[ring9members[i][3]];
                                }
                                diffVal = oPosY[ring9members[i][4]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring9members[i][4]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring9members[i][4]]+boxLength[1];
                                } else {
                                    tempPosY[4] = oPosY[ring9members[i][4]];
                                }
                                diffVal = oPosY[ring9members[i][5]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring9members[i][5]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring9members[i][5]]+boxLength[1];
                                } else {
                                    tempPosY[5] = oPosY[ring9members[i][5]];
                                }
                                diffVal = oPosY[ring9members[i][6]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring9members[i][6]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring9members[i][6]]+boxLength[1];
                                } else {
                                    tempPosY[6] = oPosY[ring9members[i][6]];
                                }
                                diffVal = oPosY[ring9members[i][7]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring9members[i][7]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring9members[i][7]]+boxLength[1];
                                } else {
                                    tempPosY[7] = oPosY[ring9members[i][7]];
                                }
                                diffVal = oPosY[ring9members[i][8]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[8] = oPosY[ring9members[i][8]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[8] = oPosY[ring9members[i][8]]+boxLength[1];
                                } else {
                                    tempPosY[8] = oPosY[ring9members[i][8]];
                                }
                                yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4]+tempPosY[5]+tempPosY[6]+tempPosY[7]+tempPosY[8];
                                yVal /= 9;

                                tempPosZ[0]=oPosZ[ring9members[i][0]];
                                diffVal = oPosZ[ring9members[i][1]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring9members[i][1]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring9members[i][1]]+boxLength[2];
                                } else {
                                    tempPosZ[1] = oPosZ[ring9members[i][1]];
                                }
                                diffVal = oPosZ[ring9members[i][2]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring9members[i][2]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring9members[i][2]]+boxLength[2];
                                } else {
                                    tempPosZ[2] = oPosZ[ring9members[i][2]];
                                }
                                diffVal = oPosZ[ring9members[i][3]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring9members[i][3]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring9members[i][3]]+boxLength[2];
                                } else {
                                    tempPosZ[3] = oPosZ[ring9members[i][3]];
                                }
                                diffVal = oPosZ[ring9members[i][4]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring9members[i][4]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring9members[i][4]]+boxLength[2];
                                } else {
                                    tempPosZ[4] = oPosZ[ring9members[i][4]];
                                }
                                diffVal = oPosZ[ring9members[i][5]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring9members[i][5]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring9members[i][5]]+boxLength[2];
                                } else {
                                    tempPosZ[5] = oPosZ[ring9members[i][5]];
                                }
                                diffVal = oPosZ[ring9members[i][6]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring9members[i][6]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring9members[i][6]]+boxLength[2];
                                } else {
                                    tempPosZ[6] = oPosZ[ring9members[i][6]];
                                }
                                diffVal = oPosZ[ring9members[i][7]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring9members[i][7]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring9members[i][7]]+boxLength[2];
                                } else {
                                    tempPosZ[7] = oPosZ[ring9members[i][7]];
                                }
                                diffVal = oPosZ[ring9members[i][8]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[8] = oPosZ[ring9members[i][8]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[8] = oPosZ[ring9members[i][8]]+boxLength[2];
                                } else {
                                    tempPosZ[8] = oPosZ[ring9members[i][8]];
                                }
                                zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4]+tempPosZ[5]+tempPosZ[6]+tempPosZ[7]+tempPosZ[8];
                                zVal /= 9;
                                //outputer_xyz << setw(5) << "O" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
                                if (ringTrajOutOpt){
                                    outputer_xyz <<left<< setw(5) << "F" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                                }

                                if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                                    red_val = "0.5";
                                    green_val = "1.0";
                                    blue_val = "0.0";
                                    zVal = -0.3;
                                    pov_out << "ring9(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    if (xVal < (-0.5*boxLength[0]+buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring9(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring9(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring9(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring9(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring9(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring9(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring9(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                        pov_out << "ring9(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                        pov_out << "ring9(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    }
                                }
                            }
                        if (maxRingOpt > 9){
                            for (i=0; i<ring10members.size(); i++){
                                tempPosX[0]=oPosX[ring10members[i][0]];
                                diffVal = oPosX[ring10members[i][1]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring10members[i][1]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[1] = oPosX[ring10members[i][1]]+boxLength[0];
                                } else {
                                    tempPosX[1] = oPosX[ring10members[i][1]];
                                }
                                diffVal = oPosX[ring10members[i][2]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring10members[i][2]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[2] = oPosX[ring10members[i][2]]+boxLength[0];
                                } else {
                                    tempPosX[2] = oPosX[ring10members[i][2]];
                                }
                                diffVal = oPosX[ring10members[i][3]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring10members[i][3]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[3] = oPosX[ring10members[i][3]]+boxLength[0];
                                } else {
                                    tempPosX[3] = oPosX[ring10members[i][3]];
                                }
                                diffVal = oPosX[ring10members[i][4]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring10members[i][4]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[4] = oPosX[ring10members[i][4]]+boxLength[0];
                                } else {
                                    tempPosX[4] = oPosX[ring10members[i][4]];
                                }
                                diffVal = oPosX[ring10members[i][5]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring10members[i][5]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[5] = oPosX[ring10members[i][5]]+boxLength[0];
                                } else {
                                    tempPosX[5] = oPosX[ring10members[i][5]];
                                }
                                diffVal = oPosX[ring10members[i][6]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring10members[i][6]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[6] = oPosX[ring10members[i][6]]+boxLength[0];
                                } else {
                                    tempPosX[6] = oPosX[ring10members[i][6]];
                                }
                                diffVal = oPosX[ring10members[i][7]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring10members[i][7]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[7] = oPosX[ring10members[i][7]]+boxLength[0];
                                } else {
                                    tempPosX[7] = oPosX[ring10members[i][7]];
                                }
                                diffVal = oPosX[ring10members[i][8]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[8] = oPosX[ring10members[i][8]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[8] = oPosX[ring10members[i][8]]+boxLength[0];
                                } else {
                                    tempPosX[8] = oPosX[ring10members[i][8]];
                                }
                                diffVal = oPosX[ring10members[i][9]]-tempPosX[0];
                                if (diffVal > 0.5*boxLength[0]){
                                    tempPosX[9] = oPosX[ring10members[i][9]]-boxLength[0];
                                } else if (diffVal < -0.5*boxLength[0]){
                                    tempPosX[9] = oPosX[ring10members[i][9]]+boxLength[0];
                                } else {
                                    tempPosX[9] = oPosX[ring10members[i][9]];
                                }
                                xVal = tempPosX[0]+tempPosX[1]+tempPosX[2]+tempPosX[3]+tempPosX[4]+tempPosX[5]+tempPosX[6]+tempPosX[7]+tempPosX[8]+tempPosX[9];
                                xVal /= 10;

                                tempPosY[0]=oPosY[ring10members[i][0]];
                                diffVal = oPosY[ring10members[i][1]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring10members[i][1]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[1] = oPosY[ring10members[i][1]]+boxLength[1];
                                } else {
                                    tempPosY[1] = oPosY[ring10members[i][1]];
                                }
                                diffVal = oPosY[ring10members[i][2]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring10members[i][2]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[2] = oPosY[ring10members[i][2]]+boxLength[1];
                                } else {
                                    tempPosY[2] = oPosY[ring10members[i][2]];
                                }
                                diffVal = oPosY[ring10members[i][3]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring10members[i][3]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[3] = oPosY[ring10members[i][3]]+boxLength[1];
                                } else {
                                    tempPosY[3] = oPosY[ring10members[i][3]];
                                }
                                diffVal = oPosY[ring10members[i][4]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring10members[i][4]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[4] = oPosY[ring10members[i][4]]+boxLength[1];
                                } else {
                                    tempPosY[4] = oPosY[ring10members[i][4]];
                                }
                                diffVal = oPosY[ring10members[i][5]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring10members[i][5]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[5] = oPosY[ring10members[i][5]]+boxLength[1];
                                } else {
                                    tempPosY[5] = oPosY[ring10members[i][5]];
                                }
                                diffVal = oPosY[ring10members[i][6]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring10members[i][6]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[6] = oPosY[ring10members[i][6]]+boxLength[1];
                                } else {
                                    tempPosY[6] = oPosY[ring10members[i][6]];
                                }
                                diffVal = oPosY[ring10members[i][7]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring10members[i][7]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[7] = oPosY[ring10members[i][7]]+boxLength[1];
                                } else {
                                    tempPosY[7] = oPosY[ring10members[i][7]];
                                }
                                diffVal = oPosY[ring10members[i][8]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[8] = oPosY[ring10members[i][8]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[8] = oPosY[ring10members[i][8]]+boxLength[1];
                                } else {
                                    tempPosY[8] = oPosY[ring10members[i][8]];
                                }
                                diffVal = oPosY[ring10members[i][9]]-tempPosY[0];
                                if (diffVal > 0.5*boxLength[1]){
                                    tempPosY[9] = oPosY[ring10members[i][9]]-boxLength[1];
                                } else if (diffVal < -0.5*boxLength[1]){
                                    tempPosY[9] = oPosY[ring10members[i][9]]+boxLength[1];
                                } else {
                                    tempPosY[9] = oPosY[ring10members[i][9]];
                                }
                                yVal = tempPosY[0]+tempPosY[1]+tempPosY[2]+tempPosY[3]+tempPosY[4]+tempPosY[5]+tempPosY[6]+tempPosY[7]+tempPosY[8]+tempPosY[9];
                                yVal /= 10;

                                tempPosZ[0]=oPosZ[ring10members[i][0]];
                                diffVal = oPosZ[ring10members[i][1]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring10members[i][1]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[1] = oPosZ[ring10members[i][1]]+boxLength[2];
                                } else {
                                    tempPosZ[1] = oPosZ[ring10members[i][1]];
                                }
                                diffVal = oPosZ[ring10members[i][2]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring10members[i][2]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[2] = oPosZ[ring10members[i][2]]+boxLength[2];
                                } else {
                                    tempPosZ[2] = oPosZ[ring10members[i][2]];
                                }
                                diffVal = oPosZ[ring10members[i][3]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring10members[i][3]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[3] = oPosZ[ring10members[i][3]]+boxLength[2];
                                } else {
                                    tempPosZ[3] = oPosZ[ring10members[i][3]];
                                }
                                diffVal = oPosZ[ring10members[i][4]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring10members[i][4]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[4] = oPosZ[ring10members[i][4]]+boxLength[2];
                                } else {
                                    tempPosZ[4] = oPosZ[ring10members[i][4]];
                                }
                                diffVal = oPosZ[ring10members[i][5]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring10members[i][5]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[5] = oPosZ[ring10members[i][5]]+boxLength[2];
                                } else {
                                    tempPosZ[5] = oPosZ[ring10members[i][5]];
                                }
                                diffVal = oPosZ[ring10members[i][6]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring10members[i][6]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[6] = oPosZ[ring10members[i][6]]+boxLength[2];
                                } else {
                                    tempPosZ[6] = oPosZ[ring10members[i][6]];
                                }
                                diffVal = oPosZ[ring10members[i][7]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring10members[i][7]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[7] = oPosZ[ring10members[i][7]]+boxLength[2];
                                } else {
                                    tempPosZ[7] = oPosZ[ring10members[i][7]];
                                }
                                diffVal = oPosZ[ring10members[i][8]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[8] = oPosZ[ring10members[i][8]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[8] = oPosZ[ring10members[i][8]]+boxLength[2];
                                } else {
                                    tempPosZ[8] = oPosZ[ring10members[i][8]];
                                }
                                diffVal = oPosZ[ring10members[i][9]]-tempPosZ[0];
                                if (diffVal > 0.5*boxLength[2]){
                                    tempPosZ[9] = oPosZ[ring10members[i][9]]-boxLength[2];
                                } else if (diffVal < -0.5*boxLength[2]){
                                    tempPosZ[9] = oPosZ[ring10members[i][9]]+boxLength[2];
                                } else {
                                    tempPosZ[9] = oPosZ[ring10members[i][9]];
                                }
                                zVal = tempPosZ[0]+tempPosZ[1]+tempPosZ[2]+tempPosZ[3]+tempPosZ[4]+tempPosZ[5]+tempPosZ[6]+tempPosZ[7]+tempPosZ[8]+tempPosZ[9];
                                zVal /= 10;
                                //outputer_xyz << setw(5) << "O" <<" "<< setw(12) << xVal <<" "<< setw(12) << yVal <<" "<< setw(12) << zVal << "\n";
                                if (ringTrajOutOpt){
                                    outputer_xyz <<left<< setw(5) << "Ne" <<" "<<right<<fixed<<setprecision(6)<<setw(10) << xVal <<" "<<right<<setw(10) << yVal <<" "<<right<<setw(10) << zVal << "\n";
                                }

                                if (povrayOutOpt && zVal <= slab_thickness && zVal >= -slab_thickness){
                                    red_val = "1.0";
                                    green_val = "0.0";
                                    blue_val = "0.5";
                                    zVal = -0.35;

                                    pov_out << "ring10(" << xVal << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    if (xVal < (-0.5*boxLength[0]+buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring10(" << xVal+boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring10(" << xVal+boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring10(" << xVal+boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (xVal > (0.5*boxLength[0]-buffer_size)){
                                        if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                            pov_out << "ring10(" << xVal-boxLength[0] << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                            pov_out << "ring10(" << xVal-boxLength[0] << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                            pov_out << "ring10(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        } else {
                                            pov_out << "ring10(" << xVal-boxLength[0] << ", " << yVal << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                        }
                                    } else if (yVal < (-0.5*boxLength[1]+buffer_size)){
                                        pov_out << "ring10(" << xVal << ", " << yVal+boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    } else if (yVal > (0.5*boxLength[1]-buffer_size)){
                                        pov_out << "ring10(" << xVal << ", " << yVal-boxLength[1] << ", " << zVal << ", " << red_val << ", " << green_val << ", " << blue_val << ", " << transparency << ")\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
	}
    }

            if (povrayOutOpt){
                x_frame = int(scale_factor * boxLength[0]);
                y_frame = int(scale_factor * boxLength[1]);
                xFrameInt.str("");
                yFrameInt.str("");
                xFrameInt << x_frame;
                yFrameInt << y_frame;
                xFrameString = xFrameInt.str();
                yFrameString = yFrameInt.str();
                outputer_pov << "povray -w" + xFrameString + " -h" + yFrameString + " +a0.1 -D " + povName2 + "_" + frameCountString + ".pov\n";
            }
            frameCount++;

            // close the pov file handle
            pov_out.close();
            povDistOut.close();
        }

        // clear reused vectors
        oPosX.clear();
        oPosY.clear();
        oPosZ.clear();
        hPosX.clear();
        hPosY.clear();
        hPosZ.clear();
        neighborListIndex.clear();
        neighborList.clear();
        largeNeighborListIndex.clear();
        largeNeighborList.clear();
        nearestNeighborList.clear();
        nearestNeighborMag.clear();
        nearestNeighborVecX.clear();
        nearestNeighborVecY.clear();
        nearestNeighborVecZ.clear();
        tetrahedrality.clear();
        hbondListIndex.clear();
        hbondList.clear();
        hbondVecX.clear();
        hbondVecY.clear();
        hbondVecZ.clear();
        ring3members.clear();
        ring4members.clear();
        ring5members.clear();
        ring6members.clear();
        ring7members.clear();
        ring8members.clear();
        ring9members.clear();
        ring10members.clear();
        hbondCount = 0;
        lastRing3 = 0;
        lastRing4 = 0;
        lastRing5 = 0;
        lastRing6 = 0;
        lastRing7 = 0;
        lastRing8 = 0;
        lastRing9 = 0;
        lastRing10 = 0;
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

    cout << "\nRing distribution results written to " << fileName << "\n\n";
    if (tetraPDBOutOpt){
	cout << "A pdb file has been generated as " << strungName4 << "\n";
    }
    if (povrayOutOpt){
        cout << "POVRay rendering command written to " << strungName3 << "\n";
        cout << "\t...potentially useful for rendering files in the NEWLY made pov_files directory\n\n";
    }
    if (ringTrajOutOpt){
        cout << "Finally, a ring trajectory file was written to " << strungName2 << "\n\n";
    }
    return 0;
}


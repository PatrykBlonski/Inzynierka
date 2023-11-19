/*
  ==============================================================================

    Calculate.c
    Created: 18 Nov 2023 9:04:37pm
    Author:  igosie

  ==============================================================================
*/

#include "CalculateXY.h"

#include <math.h>

#define PI 3.14159265358979323846

// Function to convert degrees to radians
float toRadians(float degrees) 
{
    return degrees * (PI / 180.0);
}

void calculateCoordinates(float distance, float azimuth, float* x, float* y) //we get x and y
{
    float azimuthRadians = toRadians(azimuth);
    *x = distance * cos(azimuthRadians);
    *y = distance * sin(azimuthRadians);
}
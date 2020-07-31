/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "geos95.h"
#include <netcdf.h>
#include <math.h>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int calc_adjacent_vector_indices_h(int from_index[], int to_index[], int adjacent_signs_h[], int adjacent_vector_indices_h[])
{
    int trouble_detected = 0;
    int counter;
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            if (from_index[j] == i || to_index[j] == i)
            {
                if (from_index[j] == to_index[j])
				{
                    printf("It is from_index == to_index at some point.\n");
					exit(1);
				}
                adjacent_vector_indices_h[6*i + counter] = j;
                if (from_index[j] == i)
                    adjacent_signs_h[6*i + counter] = 1;
                if (to_index[j] == i)
                    adjacent_signs_h[6*i + counter] = -1;
                ++counter;
            }
        }
        if (counter != 6)
        {
            trouble_detected = 1;
            if (counter == 5 && i < NO_OF_PENTAGONS)
                trouble_detected = 0;
        }
        if (trouble_detected == 1)
		{
            printf("Trouble detected, place 1.\n");
			exit(1);
		}
        if (i < NO_OF_PENTAGONS)
        {
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
        }
    }
    int no_of_edges, double_check, sign_sum_check;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        counter = 0;
        sign_sum_check = 0;
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            no_of_edges = 6;
            if (j < NO_OF_PENTAGONS)
                no_of_edges = 5;
            double_check = 0;
            for (int k = 0; k < no_of_edges; ++k)
            {
                if (adjacent_vector_indices_h[6*j + k] == i)
                {
                    ++counter;
                    ++double_check;
                    sign_sum_check += adjacent_signs_h[6*j + k];
                }
            }
            if (double_check > 1)
			{
                printf("Same vector twice in adjacent_vector_indices_h of same grid cell.\n");
				exit(1);
			}
        }
        if (sign_sum_check != 0)
            printf("Problem with adjacent_signs_h.\n");
        if (counter != 2)
            printf("Problem with adjacent_vector_indices_h.\n");
    }
    return 0;
}
    
int set_horizontal_curl_indices(double direction_dual[], double direction[], int h_curl_indices[], int from_index[], int to_index[], double ORTH_CRITERION_DEG, int h_curl_signs[])
{
	int sign, counter;
	double direction_change;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        sign = 1;
        find_angle_change(direction_dual[i], direction[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            sign = -1;
        h_curl_indices[4*i + 0] = i + NO_OF_VECTORS_PER_LAYER;
        h_curl_signs[4*i + 0] = sign;
        if (sign == 1)
            h_curl_indices[4*i + 1] = to_index[i];
        else
            h_curl_indices[4*i + 1] = from_index[i];
        h_curl_signs[4*i + 1] = 1;
        h_curl_indices[4*i + 2] = i;
        h_curl_signs[4*i + 2] = -sign;
        if (sign == 1)
            h_curl_indices[4*i + 3] = from_index[i];
        else
            h_curl_indices[4*i + 3] = to_index[i];
        h_curl_signs[4*i + 3] = -1;
    }
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NO_OF_VECTORS_H; ++j)
    	{
    		if (h_curl_indices[4*j + 1] == i || h_curl_indices[4*j + 3] == i)
    			++counter;
    	}
    	if (i < NO_OF_PENTAGONS && counter != 5)
    	{
    		printf("Error in h_curl_indices, position 0.\n");
    		exit(1);
		}
    	if (i >= NO_OF_PENTAGONS && counter != 6)
    	{
    		printf("Error in h_curl_indices, position 1.\n");
    		exit(1);
		}
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NO_OF_VECTORS_H; ++j)
    	{
    		if (h_curl_indices[4*j + 0] == i + NO_OF_VECTORS_PER_LAYER || h_curl_indices[4*j + 2] == i)
    		{
    			++counter;
    			if (h_curl_indices[4*j + 0] == i + NO_OF_VECTORS_PER_LAYER && h_curl_indices[4*j + 2] == i)
    				++counter;
			}
    	}
    	if (counter != 2)
    	{
    		printf("Error in h_curl_indices, position 2.\n");
    		exit(1);
		}
    }
	return 0;
}
    
int set_vorticity_indices(int vorticity_indices_pre[], int vorticity_signs_pre[], int from_index_dual[], int to_index_dual[], int vorticity_indices[], int vorticity_signs[])
{
	int counter;
    int indices_list_pre[6];
    int signs_list_pre[6];
    int indices_list[4];
    int signs_list[4];
    int double_indices[2];
    double_indices[0] = -1;
    double_indices[1] = -1;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[j] = vorticity_indices_pre[3*to_index_dual[i] + j];
			signs_list_pre[j] = vorticity_signs_pre[3*to_index_dual[i] + j];
    	}
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[3 + j] = vorticity_indices_pre[3*from_index_dual[i] + j];
			signs_list_pre[3 + j] = vorticity_signs_pre[3*from_index_dual[i] + j];
    	}
		for (int j = 0; j < 6; ++j)
		{
			for (int k = j + 1; k < 6; ++k)
			{
				if (indices_list_pre[j] == indices_list_pre[k])
				{
					double_indices[0] = j;
					double_indices[1] = k;
				}
			}
		}
		counter = 0;
    	for (int j = 0; j < 6; ++j)
    	{
    		if (j != double_indices[0] && j != double_indices[1])
    		{
				indices_list[counter] = indices_list_pre[j];
				signs_list[counter] = signs_list_pre[j];
				counter++;
			}
		}
		if (counter != 4)
		{
			printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 1.\n");
			exit(1);
		}
    	for (int j = 0; j < 4; ++j)
    	{
			vorticity_indices[4*i + j] = indices_list[j];
			vorticity_signs[4*i + j] = signs_list[j];
			if (vorticity_signs[4*i + j] != 1 && vorticity_signs[4*i + j] != -1)
			{
				printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 2.");
				exit(1);
			}
			if (vorticity_indices[4*i + j] >= NO_OF_VECTORS_H || vorticity_indices[4*i + j] < 0)
			{
				printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 3.");
				exit(1);
			}
		}
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NO_OF_VECTORS_H; ++j)
    	{
    		for (int k = 0; k < 4; ++k)
    		{
    			if (vorticity_indices[4*j + k] == i)
    				++counter;
    		}
    	}
    	if (counter != 4)
    	{
    		printf("Error in vorticity_indices, position 0.\n");
    		exit(1);
    	}
    }
	return 0;
}
    
    
    
    
    
    
    

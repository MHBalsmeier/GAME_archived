/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains helper functions concerned with simple algebraic operations on vectors.
*/

#include <math.h>
#include "grid_generator.h"

int flinspace(double in_vector[], int vector_length, double min, double max)
{
	/*
	This function separates the interval between a minimum and a maximum into equidistant sections.
	*/
    double delta = (max - min)/(vector_length - 1.0);
    for (int i = 0; i < vector_length; i++)
    {
        in_vector[i] = min + i*delta;
        if(i == vector_length - 1)
        {
            in_vector[i] = max;
    	}
    }
	return 0;
}

int fwhere_double(double vector[], int vector_length, double wanted_value, double epsilon)
{
	/*
	This function returns the index where a vector of doubles (almost) assumes a certain value.
	*/
    int index = -1;
    int continue_bool = 1;
    int i = 0;
    while(continue_bool == 1)
    {
        if(fabs(vector[i] - wanted_value) < epsilon)
        {
            index = i;
            continue_bool = 0;
        }
        i++;
        if(i == vector_length)
        {
            continue_bool = 0;
    	}
    }
    return index;
}

int find_min_index(double vector[], int vector_length)
{
	/*
	This function returns the index where a vector has its minimum.
	*/
    int result = 0;
    double current_min = vector[0];
    for (int i = 1; i < vector_length; ++i)
    {
        if (vector[i] < current_min)
        {
            current_min = vector[i];
            result = i;
        }
    }
    return result;
}

int find_max_index(double vector[], int vector_length)
{
	/*
	This function returns the index where a vector has its maximum.
	*/
    int result = 0;
    double current_max = vector[0];
    for (int i = 1; i < vector_length; i++)
    {
        if(vector[i] > current_max)
        {
            current_max = vector[i];
            result = i;
        }
    }
    return result;
}

int find_min_index_exclude(double vector[], int vector_length, int exclude_indices_vector[], int exclude_indices_vector_length)
{
	/*
	This function finds the index where a vector has its minimum, excluding the elements of another vector.
	*/
    int result = 0;
    double current_min = vector[0];
    for (int i = 1; i < vector_length; ++i)
    {
        if (vector[i] < current_min)
        {
        	if (in_bool_calculator(i, exclude_indices_vector, exclude_indices_vector_length) == 0)
        	{
            	current_min = vector[i];
            	result = i;
        	}
        }
    }
    return result;
}

int fmake_abs(double vector[], int vector_length)
{
	/*
	This function turns each element of a vector of doubles into its absolute value.
	*/
    for (int i = 0; i < vector_length; i++)
    {
        vector[i] = fabs(vector[i]);
	}
	return 0;
}

int fvec_scal(double result[], double vector[], double scalar, int vector_length)
{
	/*
	This function adds a scalar value to a vector.
	*/
    for (int i = 0; i < vector_length; i++)
    {
        result[i] = vector[i] + scalar;
	}
	return 0;
}

int find_index_n_equal(double time_vector_a[], double time_vector_b[], int n_time_a, int n_time_b)
{
	/*
	This function detects how many values of two vectors of doubles are equal.
	*/
    int n_overlap = 0;
    for (int i = 0; i < n_time_a; ++i)
    {
        for (int j = 0; j < n_time_b; ++j)
        {
            if (time_vector_a[i] == time_vector_b[j])
            {
                n_overlap = n_overlap + 1;
        	}
        }
    }
    return n_overlap;
}

int find_index_where_equal_a(double time_vector_a[], double time_vector_b[], int n_time_a, int n_time_b, int time_indices_overlap_a[])
{
	/*
	This function detects the indices of vector_a where it has a value that is an element of vector_b.
	*/
    int index = 0;
    for (int i = 0; i < n_time_a; ++i)
    {
        for (int j = 0; j < n_time_b; ++j)
        {
            if (time_vector_a[i] == time_vector_b[j])
            {
                time_indices_overlap_a[index] = i;
                ++index;
            }
        }
    }
    return 0;
}

int find_n_between_points(double vector[], int n_values, double min_value, double max_value)
{
	/*
	This function counts the elements of a vector that are between two given values.
	*/
    int n_points_between = 0;
    for (int i = 0; i < n_values; ++i)
    {
        if (vector[i] > min_value && vector[i] < max_value)
        {
            n_points_between = n_points_between + 1;
    	}
    }
    return n_points_between;
}

double fsum(double vector[], int vector_length)
{
	/*
	This function calculates the sum of a vector.
	*/
	double result = 0.0;
	for (int i = 0; i < vector_length; ++i)
	{
		result = result + vector[i];
	}
	return result;
}

int freverse(double vector_in[], int vector_length, double vector_out[])
{
	/*
	This function reverses a vector of doubles.
	*/
	for (int i = 0; i < vector_length; ++i)
	{
		vector_out[i] = vector_in[((int) vector_length) - 1 - i];
	}
	return 0;
}

int freverse_int(int vector_in[], int vector_length, int vector_out[])
{
	/*
	This function reverses a vector of integers.
	*/
	for (int i = 0; i < vector_length; ++i)
	{
		vector_out[i] = vector_in[((int) vector_length) - 1 - i];
	}
	return 0;
}

double double_sum_gen(double vector[], int vector_length, int first_index, int second_index)
{
	/*
	This function calculates the sum of all elements of a vector of doubles between two indices.
	*/
	double result = 0.0;
	if (first_index <= second_index)
	{
		for (int i = first_index ; i <= second_index; ++i)
		{
			result += vector[i];
		}
	}
	else
	{
		for (int i = first_index; i < vector_length; ++i)
		{
			result += vector[i];
		}
		for (int i = 0; i <= second_index; ++i)
		{
			result += vector[i];
		}	
	}
	return result;
}

int in_bool_calculator(int value, int check_array[], int array_length)
{
	/*
	This function checks if a vector of integers contains a certain value.
	*/
    int result = 0;
    for (int i = 0; i < array_length; i++)
    {
        if(check_array[i] == value)
        {
            result = 1;
            break;
    	}
    }
    return result;
}










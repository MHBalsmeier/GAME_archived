#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "io.h"

void write_out(State state_write_out, double t_init, double t_write, int write_out_index, char output_foldername[])
{
	Scalar_field pressure, temperature;
	pressure_diagnostics(state_write_out.pot_temp,state_write_out.density,pressure);
	temperature_diagnostics(state_write_out.pot_temp,pressure,temperature);
	FILE *output_file;
	int number_of_digits = 1;
	if(write_out_index > 9)
	{
		number_of_digits = 2;
	}
	if(write_out_index > 99)
		number_of_digits = 3;
	if(write_out_index > 999)
		number_of_digits = 4;
	if(write_out_index > 9999)
		number_of_digits = 5;
	char number[number_of_digits];
	itoa(write_out_index, number, 10);
	char output_folder_pre[] = "C:/Users/Max/Desktop/Programme/C/ESS_1_0_0/output/";
	char output_filename[strlen(output_folder_pre)+strlen(output_foldername)+1+number_of_digits+1];
	for (int i = 0; i < strlen(output_folder_pre)+strlen(output_foldername)+1+number_of_digits+1; ++i)
	{
		if(i<strlen(output_folder_pre))
			output_filename[i] = output_folder_pre[i];
		else if(i<strlen(output_folder_pre)+strlen(output_foldername))
			output_filename[i] = output_foldername[i-strlen(output_folder_pre)];
		else if(i<strlen(output_folder_pre)+strlen(output_foldername)+1)
			output_filename[i] = '/';
		else
			output_filename[i] = number[i-(strlen(output_folder_pre)+strlen(output_foldername)+1)];
	}
	output_filename[strlen(output_folder_pre)+strlen(output_foldername)+1+number_of_digits] = '\0';
	output_file = fopen(output_filename, "w");
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		fprintf(output_file,"%d %f %f %f\n",i,pressure[i],temperature[i],state_write_out.density[i]);
	}
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		fprintf(output_file,"%d %f\n",i,state_write_out.wind[i]);
	}
	fclose(output_file);
}

#include <math.h>
#include "../enum_and_typedefs.h"
#include <stdio.h>
#include <string.h>
#include "io.h"

int get_line_length(int);
int get_add_length(int);

State initializer(char input_filename[])
{
	extern double R_d, c_p, p_0;
	FILE *input_file;
	State init_state;
	char input_folder_pre[] = "C:/Users/Max/Desktop/Programme/C/ESS_1_0_0/input/";
	char input_filename_with_path[strlen(input_folder_pre)+strlen(input_filename)];
	for (int i = 0; i < strlen(input_folder_pre)+strlen(input_filename); ++i)
	{
		if(i<strlen(input_folder_pre))
			input_filename_with_path[i] = input_folder_pre[i];
		else if(i<strlen(input_folder_pre)+strlen(input_filename))
			input_filename_with_path[i] = input_filename[i-strlen(input_folder_pre)];
	}
	input_file = fopen(input_filename_with_path, "r");
	for (int i = 0 ;i < NUMBER_OF_SCALARS+NUMBER_OF_VECTORS; ++i)
	{
		int number_of_characters = get_line_length(i)+1;
		char line[number_of_characters];
		if(fgets(line, 61, input_file))
		{
			int line_checker;
			if (i<NUMBER_OF_SCALARS)
			{
				int line_checker;
				double pressure;
				double temp;
				double density;
				sscanf(line,"%ld %lf %lf %lf\n\0",&line_checker,&pressure,&temp,&density);
				init_state.density[i] = density;
				init_state.pot_temp[i] = temp*pow(p_0/pressure,R_d/c_p);
			}
			if (i >= NUMBER_OF_SCALARS)
			{
				int line_checker;
				double wind;
				sscanf(line,"%ld %lf\n\0",&line_checker,&wind);
				init_state.wind[i-(NUMBER_OF_SCALARS+NUMBER_OF_SCALARS_H)] = wind;
			}
		}
	}
	fclose(input_file);
	return init_state;
}

int get_line_length(int j)
{
	int line_length, base_value, add_value = 0;
		if (j<NUMBER_OF_SCALARS)
		{
			base_value = 19+1;
			int add_value = get_add_length(j);
		}
		if (j >= NUMBER_OF_SCALARS)
		{
			base_value = 10+1;
			int add_value = get_add_length(j-(NUMBER_OF_SCALARS+NUMBER_OF_SCALARS_H));
		}
	line_length = base_value + add_value;
	return line_length;
}

int get_add_length(int k)
{
	int add_value = 0;
	if( k > (1e1)-1)
		add_value = 1;
	if( k > (1e2)-1)
		add_value = 2;
	if( k > (1e3)-1)
		add_value = 3;
	if( k > (1e4)-1)
		add_value = 4;
	if( k > (1e5)-1)
		add_value = 5;
	if( k > (1e6)-1)
		add_value = 6;
	if( k > (1e7)-1)
		add_value = 7;
	if( k > (1e8)-1)
		add_value = 8;
	return add_value;
}

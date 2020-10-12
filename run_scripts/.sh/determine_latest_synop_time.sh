#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

now=$(date +%s)
year=$(date --utc -d @$now +%Y)
month=$(date --utc -d @$now +%m)
day=$(date --utc -d @$now +%d)
hour=$(date --utc -d @$now +%H)
now_hr=$(date --utc -d "$year$month$day $hour:00:00" +%s)
hr_substract=$(python -c "print(int('$hour') - 6*int(int('$hour')/6))")
wanted_timestamp=$(($now_hr - 3600*$hr_substract))
year=$(date --utc -d @$wanted_timestamp +%Y)
month=$(date --utc -d @$wanted_timestamp +%m)
day=$(date --utc -d @$wanted_timestamp +%d)
hour=$(date --utc -d @$wanted_timestamp +%H)

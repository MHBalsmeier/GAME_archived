old_string="Geophysical"
new_string="Geophysical"
grep -rl $old_string . | xargs sed -i 's/'$old_string'/'$new_string'/g'

rm /home/max/figs/game_output/movie_elements/*
python /home/max/my_code/game_plotting/create_movie_elements.py $run_span $write_out_interval $disp_level $disp_shortname $grid_props_file $fig_save_path $core_path $output_dir
if [ -f movie.mp4 ]
then
rm movie.mp4
fi
/home/max/my_code/game_plotting/create_movie.sh

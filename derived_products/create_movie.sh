ffmpeg -framerate 25 -i /home/max/figs/game_output/movie_elements/file-%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p derived_products/movie.mp4

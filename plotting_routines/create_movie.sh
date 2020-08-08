echo creating movie ...
ffmpeg -y -hide_banner -loglevel warning -framerate 50 -i $fig_save_path/file-%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $fig_save_path/movie.mp4
if [ $? -ne 0 ]
then
echo -e ${RED}Creaing movie failed.$NC
else
echo Movie successfully created.
fi

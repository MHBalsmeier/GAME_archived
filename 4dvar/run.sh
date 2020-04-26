if  [ -f  ]
then
rm -r obs/$run_name
fi
mkdir obs/$run_name
source sh/downloader.sh

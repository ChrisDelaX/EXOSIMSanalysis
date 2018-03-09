python keepout_graphics.py -s 65000 -l 2.0 -d 5.0 -m $HOME/Desktop/keepout.mp4 ./sampleScript_coron_v2.json
python keepout_graphics.py -s 60900 -l 1.0 -d 10.0 -m $HOME/Desktop/keepout.mp4 ./sampleScript_coron_v2.json
python keepout_graphics.py -s 60900 -l 1.0 -d 100 -m $HOME/Desktop/keepout.mp4 ./sampleScript_coron_v2.json


Where:

-s is the start-time of the movie (MJD)
-l is the length in years
-d is the delta-t between frames in days
-m is the name of the movie
[optionally, -f is the directory name for .png output]
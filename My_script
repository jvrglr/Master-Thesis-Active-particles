#!/bin/bash
echo "Compile active matter..."
gfortran 2D_soft_core.f dranxor2.f Subroutines.f -o exe1
echo "Execute active matter: creating data files..."
time ./exe1
echo "Creating .png files with snapshots..."
#gnuplot multiple.plt #Create snapshots with gnuplot
python Multiplot.py  #Create snapshots with python
echo "Creating .mp4 video from snapshots..."
#ffmpeg -r 14 -i image.%d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4 #for gnuplot
ffmpeg -r 14 -i image.%d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf scale=1600:800 out.mp4 #for python

#This last lines will remove .dat and .png files after creating the video
#WARNING!!!!!!!! all .png and .dat files in the directory will be erased
echo "Deleting .dat files..."
find . -type f -name '*.dat' -delete
echo "Deleting .png files..."
find . -type f -name '*.png' -delete



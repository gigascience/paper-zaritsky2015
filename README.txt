Matlab script from Data Note by Zaritsky et al., (2015)

The Matlab script, speedKymograph.m, creates a kymograph per experiment that displays the spatiotemporal dynamics of cell speed as function of time and distance from the monolayer edge.

The script accepts as input the main directory name, goes over the .tif/.lsm data (each file holds raw data of a single time-lapse experiment) and uses the intermediately processed data in the corresponding directories to generate the kymograph. One copy of the kymograph is placed in the experimental directory and the other in a directory named speedKymographs at the main directory. Two files are generated in each of the two directories: a .mat file with the kymograph data and a .bmp image for visualization.

Please refer to the documentation in the script for further information.

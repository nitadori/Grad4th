#!/bin/sh
g++-mp-5 -Wall -O3 -std=c++11 -fopenmp -DINTGRT=DKD_2nd main.cpp -o dkd
g++-mp-5 -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDK_2nd main.cpp -o kdk
g++-mp-5 -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_2nd main.cpp -o kdkdk2
g++-mp-5 -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_4th main.cpp -o kdkdk4

#!/bin/bash

# assumes that the folder name is passed as a directory

pwd
cd ../OUT

if [ -d $1 ] 
then
    echo "Creating new indexed output directory." 

    currName=$1
    i="0"

    while [ -d $1 ]
    do
        if [ -d $currName$i ]
        then
            i=$[$i+1]
        else
            mv $1 $currName$i
        fi
    done

    mkdir $1

else
    mkdir $1
fi

cp ../src/main_prog.f90 $1

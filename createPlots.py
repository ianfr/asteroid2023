# ffmpeg -i %10d.png -c:v libx264 -vf "fps=25,format=yuv420p" out.mp4

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

folder = 'bbox'
directory = './OUT/' + folder

BOUND = 2.5

try:
    os.mkdir(directory + '_FIGURES')
except:
    print()

for filename in os.scandir(directory):
    if filename.is_file():
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        df = pd.read_csv(filename)
        # print(df.head())
        ax.scatter(df.iloc[:,1], df.iloc[:,2], df.iloc[:,3])
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_xlim([-BOUND, BOUND])
        ax.set_ylim([-BOUND, BOUND])
        ax.set_zlim([-BOUND, BOUND])
        fileNum = filename.name.replace('ast.','').replace('csv.','')

        plt.savefig(directory + '_FIGURES/' + fileNum.zfill(10) + '.png')
        plt.close()
from __future__ import division
import sys
import numpy as np
from numpy import pi
import random
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv

def update_particle(num):
    v1.set_data(PlotX[:,num], PlotY[:,num])
    v1.set_markersize(1)
    v1.set_markerfacecolor([0.2,0.7,0.7])
    
if __name__ == "__main__":
    global PlotX, PlotY
    PlotX = np.zeros((49, 200),float)
    PlotY = np.zeros((49, 200),float)
    
    f= open('location_x.csv','r')

    n=0
    for line in csv.reader(f):
        PlotX[:,n] = [float(i) for i in line[:-2]]
        n += 1
    f.close()


    g= open('location_y.csv','r')

    n=0
    for line in csv.reader(g):
        PlotY[:,n] = [float(i) for i in line[:-2]]
        n += 1
    g.close()

    l= open('tubeloc.csv','r')

    n=0
    BLOC=np.zeros((3,41),float)
    for line in csv.reader(l):
        BLOC[:,n] = line
        n += 1
    l.close()
    
    fig1 = plt.figure(figsize=(8.2, 8))
    
    l, = plt.plot([], [], 'r-')

    v1, = plt.plot([], [], 'o',alpha=0.8,lw = 0);
    fig1.gca().add_artist(v1);

    for i in range(BLOC.shape[1]):
        s = plt.Circle((BLOC[0,i],BLOC[1,i]), BLOC[2,i], alpha=0.2)
        fig1.gca().add_artist(s)       
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Steam Generator')
    #line_ani = animation.FuncAnimation(fig1, update_line, int(Steps/20), fargs=(Vapour, l), interval=10, blit=False)
    particle_ani = animation.FuncAnimation(fig1, update_particle, PlotY.shape[1], interval=10, blit=False)
    plt.show()

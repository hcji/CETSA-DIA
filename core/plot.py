# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:52:12 2020

@author: hcji
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_ms(mz, intensity, red_mz=None, blue_mz=None):
    if (red_mz is not None) or (blue_mz is not None):
        alpha = 0.5
    else:
        alpha = 0.9
    plt.vlines(mz, 0, intensity, alpha = alpha)
    plt.axhline(0, color='black')
    plt.xlim(100, 600)
    plt.xlabel('m/z')
    plt.ylabel('intensity')
    
    if red_mz is not None:
        red = []
        for i in range(len(mz)):
            if np.min(np.abs(red_mz - mz[i])) < 0.01:
                red.append(i)
        plt.vlines(mz[red], 0, intensity[red], color='red')
    
    if blue_mz is not None:
        blue = []
        for i in range(len(mz)):
            if np.min(np.abs(blue_mz - mz[i])) < 0.01:
                blue.append(i)
        plt.vlines(mz[blue], 0, intensity[blue], color='blue')
    plt.show()
    
    
def plot_chroms(rt, eic, red_eics=None, blue_eics=None):
    plt.xlabel('retention time')
    plt.ylabel('intensity')

    plt.plot(rt, eic, color='black')
    if red_eics is not None:
        for i in range(len(red_eics)):
            plt.plot(rt, red_eics[i,:], color='red')
    if blue_eics is not None:
        for i in range(len(blue_eics)):
            plt.plot(rt, blue_eics[i,:], color='blue')
    plt.show()
    
    
    
    
    
    
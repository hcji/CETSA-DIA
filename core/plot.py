# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:52:12 2020

@author: hcji
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_ms(spectrum):
    plt.figure(figsize=(6, 4), dpi=300)
    plt.vlines(spectrum['mz'], np.zeros(spectrum.shape[0]), np.array(spectrum['intensity']), 'red')
    plt.axhline(0, color='black')
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity')
    plt.show()


def plot_anno_ms(spectrum, mzs, annotations, tol=0.1):
    spectrum = pd.DataFrame(spectrum)
    spectrum.columns = ['mz', 'intensity']
    spectrum['intensity'] /= max(spectrum['intensity'])
    c_mz, c_int, c_ann = [], [], []
    for i, mz in enumerate(mzs):
        diffs = abs(spectrum['mz'] - mz)
        if min(diffs) < tol:
            c_mz.append(spectrum['mz'][np.argmin(diffs)])
            c_int.append(spectrum['intensity'][np.argmin(diffs)])
            c_ann.append(annotations[i])
    plt.figure(figsize=(6, 4), dpi=300)
    plt.vlines(spectrum['mz'], np.zeros(spectrum.shape[0]), np.array(spectrum['intensity']), 'gray')
    plt.vlines(c_mz, np.zeros(len(c_mz)), c_int, 'red')
    for i in range(len(c_mz)):
        plt.text(c_mz[i], c_int[i] + 0.025, c_ann[i], color='red', rotation='vertical')
    plt.axhline(0, color='black')
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity')
    plt.show()


def plot_compare_ms(spectrum1, spectrum2, tol=0.5):
    spectrum1 = pd.DataFrame(spectrum1)
    spectrum2 = pd.DataFrame(spectrum2)
    spectrum1.columns = ['mz', 'intensity']
    spectrum2.columns = ['mz', 'intensity']
    spectrum1['intensity'] /= np.max(spectrum1['intensity'])
    spectrum2['intensity'] /= np.max(spectrum2['intensity'])
    c_mz = []
    c_int = []
    for i in spectrum1.index:
        diffs = abs(spectrum2['mz'] - spectrum1['mz'][i])
        if min(diffs) < tol:
            c_mz.append(spectrum1['mz'][i])
            c_mz.append(spectrum2['mz'][np.argmin(diffs)])
            c_int.append(spectrum1['intensity'][i])
            c_int.append(-spectrum2['intensity'][np.argmin(diffs)])
    c_spec = pd.DataFrame({'mz': c_mz, 'intensity': c_int})
    plt.figure(figsize=(6, 6), dpi=300)
    plt.vlines(spectrum1['mz'], np.zeros(spectrum1.shape[0]), np.array(spectrum1['intensity']), 'gray')
    plt.axhline(0, color='black')
    plt.vlines(spectrum2['mz'], np.zeros(spectrum2.shape[0]), -np.array(spectrum2['intensity']), 'gray')
    plt.vlines(c_spec['mz'], np.zeros(c_spec.shape[0]), c_spec['intensity'], 'red')
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity')
    plt.show()
    
    
    
    
    
    
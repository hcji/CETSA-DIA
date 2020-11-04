# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 09:14:11 2019

@author: hcji
"""

import pymzml
import numpy as np
import pandas as pd
import tensorflow.compat.v1 as tf
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
from scipy import signal
from bisect import bisect_left, bisect_right
from prosit import model
from prosit import prediction
from prosit import constants
from prosit import tensorize
from pyopenms import AASequence, EmpiricalFormula, CoarseIsotopePatternGenerator, MSSpectrum, Deisotoper
from seq_encode.ms import ms_to_vec


def load_data(mzml, swath):
    """
    Load raw ms file and SWATH window.
    
    Parameters
    ----------
    mzml : str
        The path of the ms file. Only *.mzml is support, other format can be
        converted by ProteoWizard. Note: the peaks should be in centroid mode.
    swath : str
        The path of the SWTAH window file. Shold be *.csv format. Include columns: 
        'start_mz', 'end_mz', 'center'.
    """
    
    window = pd.read_csv(swath)
    precursors = np.array([-1.0] + list(window['center']))
    reader = pymzml.run.Reader(mzml)
    all_peaks = []
    all_rts = []
    for p in reader:
        all_rts.append(p.scan_time[0])
        all_peaks.append(p)
    all_rts = np.array(all_rts)
    all_peaks = np.array(all_peaks)
    return {'rts': all_rts,
            'precursors': precursors,
            'window': window,
            'peaks': all_peaks}


def calc_isotope(pep, charge=2):
    """
    calculate isotope distribution of a peptide.
    
    Parameters
    ----------
    pep : str
        The peptide sequence.
    charge : int
        The charge of the peptide precursor.
    """

    seq = AASequence.fromString(pep)
    seq_formula = seq.getFormula() + EmpiricalFormula("H" + str(charge))
    isotopes = seq_formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(4))

    pattern = {'isotope': [], 'abundance': []}
    s = MSSpectrum()
    for iso in isotopes.getContainer():
        iso.setMZ(iso.getMZ() / charge)
        s.push_back(iso)
        pattern['isotope'].append(iso.getMZ())
        pattern['abundance'].append(iso.getIntensity())
    return pd.DataFrame(pattern)


def predict_library(peplist):
    """
    Predict iRT and MS/MS by prosit package:
        https://github.com/kusterlab/prosit
    
    Parameters
    ----------
    peplist : str
        The path of the peptide list file. Shold be *.csv format. Include columns: 
        'modified_sequence', 'collision_energy', 'precursor_charge'.
    """
    pep = pd.read_csv(peplist)
    pept = tensorize.csv(pep)
    d_spectra = {}
    d_irt = {}

    d_spectra["graph"] = tf.Graph()
    with d_spectra["graph"].as_default():
        d_spectra["session"] = tf.Session()
        with d_spectra["session"].as_default():
            d_spectra["model"], d_spectra["config"] = model.load(
                'prosit/model/fragmentation_prediction',
                trained=True)
            d_spectra["model"].compile(optimizer="adam", loss="mse")
            
    d_irt["graph"] = tf.Graph()
    with d_irt["graph"].as_default():
        d_irt["session"] = tf.Session()
        with d_irt["session"].as_default():
            d_irt["model"], d_irt["config"] = model.load(
                'prosit/model/irt_prediction',
                trained=True)
            d_irt["model"].compile(optimizer="adam", loss="mse")
            
    pept = prediction.predict(pept, d_irt)
    pept = prediction.predict(pept, d_spectra)
    
    output = {}
    for i in range(len(pep)):
        n = pep['modified_sequence'][i] + '_'
        n += str(pep['collision_energy'][i]) + 'NCE_'
        n += str(pep['precursor_charge'][i]) + '+'
        isot = calc_isotope(pep['modified_sequence'][i], pep['precursor_charge'][i])
        irt = pept['iRT'][i][0]
        mz = pept['masses_pred'][i]
        abund = pept['intensities_pred'][i]
        keep = np.where(abund > 0.01)[0]
        mz = mz[keep]
        abund = abund[keep]
        od = np.argsort(mz)
        mz = mz[od]
        abund = abund[od]
        
        seq = pept['sequence_integer'][i]
        mw = 0
        for s in seq:
            if s == 0:
                break
            aa = list(constants.ALPHABET.keys())[list(constants.ALPHABET.values()).index(s)]
            mw += constants.AMINO_ACID[aa]
        precursor = (mw + constants.H2O + constants.PROTON * pep['precursor_charge'][i]) / pep['precursor_charge'][i]
        
        output[n] = {'irt': irt, 'mz': mz, 'intensity': abund, 'precursor': precursor, 'isotope': isot}
    return output
    

def precursor_eic(data, library, mztol=0.05, rtlength=5):
    """
    Get extracted ion profile of the targeted peptide precursor
    
    Parameters
    ----------
    data : dict
        The results returned by the *load_data* function.
    library : dict
        The library generated by the *predict_library* function.
    """
    
    ind = np.arange(0, len(data['peaks']), len(data['precursors']))
    peaks = data['peaks'][ind]
    eics = {}
    for i in library.keys():
        exmz = library[i]['precursor']
        mzrange = [exmz - mztol, exmz + mztol]
        if 'corrected_rt' in library[i]:
            exrt = library[i]['corrected_rt']
            rtrange = [exrt - rtlength, exrt + rtlength]
        else:
            rtrange = [0, float('inf')]
            
        rts, abunds = [], []
        for p in peaks:
            if p.scan_time[0] > rtrange[1]:
                break
            if (p.scan_time[0] >= rtrange[0]) and (p.scan_time[0] <= rtrange[1]) and (p.ms_level==1):
                rts.append(p.scan_time[0])
                mzs, intensities = p.centroidedPeaks[:,0], p.centroidedPeaks[:,1]
                sel_peaks = np.arange(bisect_left(mzs, mzrange[0]), bisect_right(mzs, mzrange[1]))
                abunds.append(np.sum(intensities[sel_peaks]))
            # plt.plot(rts, abunds)
        eics[i] = {'rt': rts, 'intensity': abunds}
    return eics


def get_ms2(data, exmz, exrt):
    """
    Get the ms/ms spectrum of a targeted scan and SWATH window.
    
    Parameters
    ----------
    data : dict
        The results returned by the *load_data* function.
    exmz : float
        The targeted precursor located in a SWATH window.
    exrt : list
        A list of the targeted retention time.
    """
    window = data['window']
    wid = np.where(np.logical_and(window['start_mz'] < exmz, window['end_mz'] > exmz))[0][0] + 1
    ind1 = np.arange(wid, len(data['peaks']), len(data['precursors']))
    
    output = {}
    rts = data['rts'][ind1]
    for r in exrt:    
        rid = np.argmin(np.abs(rts - r))
        peak = data['peaks'][ind1[rid]].centroidedPeaks
        output[str(r)] = peak
    return output


def get_ms1(data, exrt):
    """
    Get the ms1 spectrum of a targeted scan.

    Parameters
    ----------
    data : dict
        The results returned by the *load_data* function.
    exrt : list
        The targeted retention time.
    """

    ind = np.arange(0, len(data['peaks']), len(data['precursors']))
    output = {}
    rts = data['rts'][ind]
    for r in exrt:
        rid = np.argmin(np.abs(rts - r))
        peak = data['peaks'][ind[rid]].centroidedPeaks
        output[str(r)] = peak
    return output


def dpscore(v1, v2):
    return np.dot(v1, v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))


def irt_curve(data, irtlib, irteic, mscor=0.85):
    """
    Get the calibration curve based on the iRT peptides
    
    Parameters
    ----------
    data : dict
        The results returned by the *load_data* function.
    irtlib : dict
        The iRT peptide library generated by the *predict_library* function.
    mscor : float
        The threshold of ms/ms similairty between the ms in library and the ms in raw data.
    """
    allrt, allirt, allscore = [], [], []
    mtv = ms_to_vec(max_mz = 2000)
    for i in range(len(irteic)):
        k = list(irteic.keys())[i]
        eic = irteic[k]
        rts = np.array(eic['rt'])
        intensities = np.array(eic['intensity'])
        pks = signal.find_peaks(intensities, width = 2, height = np.percentile(intensities, 90))[0]
        exrt = rts[pks]
        exmz = irtlib[k]['precursor']
        fragms = np.array([irtlib[k]['mz'], irtlib[k]['intensity']]).transpose()
        precms = np.array(irtlib[k]['isotope'])
        fragvec = mtv.transform(fragms)
        precvec = mtv.transform(precms)
        candms1 = get_ms1(data, exrt)
        candms2 = get_ms2(data, exmz, exrt)

        compress1 = np.where(precvec == 0)[0]
        compress2 = np.where(fragvec == 0)[0]
        rt, score = [], []
        for r, ms2 in candms2.items():
            ms1 = candms1[r]
            candvec1 = mtv.transform(ms1)
            candvec2 = mtv.transform(ms2)
            candvec1[compress1] = 0
            candvec1 /= max(candvec1)
            candvec2[compress2] = 0
            candvec2 /= max(candvec2)
            sc1 = dpscore(precvec, candvec1)
            sc2 = dpscore(fragvec, candvec2)
            sc = (sc1 + sc2) / 2
            if sc >= mscor:
                rt.append(float(r))
                score.append(sc)
        allrt.append(rt)
        allirt.append(irtlib[k]['irt'])
        allscore.append(score)

    irtl, rtl = [], []
    for i in range(len(allscore)):
        k = list(irtlib.keys())[i]
        scores = allscore[i]
        irt = irtlib[k]['irt']
        if len(scores) == 0:
            continue
        else:
            rt = allrt[i][np.argmax(allscore[i])]
        irtl.append(irt)
        rtl.append(rt)

    rtl = np.array(rtl)
    irtl = np.array(irtl).reshape(1,-1).transpose()
    lm = LinearRegression()
    lm.fit(irtl, rtl)
    rtp = lm.predict(np.array(allirt).reshape(1,-1).transpose())

    irtl, rtl = [], []
    for i in range(len(allscore)):
        k = list(irtlib.keys())[i]
        scores = allscore[i] + (1 - np.abs(allrt[i] - rtp[i]) / rtp[i]) * 0.5
        irt = irtlib[k]['irt']
        if len(scores) == 0:
            continue
        else:
            rt = allrt[i][np.argmax(scores)]
        irtl.append(irt)
        rtl.append(rt)

    rtl = np.array(rtl)
    irtl = np.array(irtl).reshape(1,-1).transpose()
    lm = LinearRegression()
    lm.fit(irtl, rtl)
    b0 = lm.intercept_
    b1 = lm.coef_[0]
    r2 = lm.score(irtl, rtl)
    return b0, b1, r2


def rt_correct(peplib, irtcurve):
    """
    Get the corrected rt of peptides.

    Parameters
    ----------
    peplib : dict
        The peptide library generated by the *predict_library* function.
    irtcurve : tuple
        The results returned by *irt_curve* function.
    """
    for i in range(len(peplib)):
        k = list(peplib.keys())[i]
        irt = peplib[k]['irt']
        rt = irtcurve[1] * irt + irtcurve[0]
        peplib[k]['corrected_rt'] = rt
    return peplib


def peptide_feature(data, peplib, irtcurve, mscor=0.85, rttol=0.3, mztol=0.05):
    """
    Get the features for classification peptides

    Parameters
    ----------
    data : dict
        The results returned by the *load_data* function.
    peplib : dict
        The peptide library generated by the *predict_library* function.
    mscor : float
        The threshold of ms/ms similairty between the ms in library and the ms in raw data.
    rttol : float
        The threshold of difference of the retention time in library and the retention time in raw data.
    """

    peplib = rt_correct(peplib, irtcurve)
    pepeic = precursor_eic(data, peplib, mztol=mztol, rtlength=5)
    for i in range(len(peplib)):
        k = list(peplib.keys())[i]

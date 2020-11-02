# -*- coding: utf-8 -*-
#
# Copyright 2014 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  When publishing work that uses OrthoFinder please cite:
#      Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically 
#      improves orthogroup inference accuracy, Genome Biology 16:157
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com

import os
import glob
import numpy as np
try:
    import cPickle as pic
except ImportError:
    import pickle as pic

from . import util, files

def DumpMatrix(name, m, iSpecies, jSpecies, d_pickle):
    with open(d_pickle + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'wb') as picFile:
        pic.dump(m, picFile, protocol=util.picProtocol)
    
def DumpMatrixArray(name, matrixArray, iSpecies, d_pickle):
    for jSpecies, m in enumerate(matrixArray):
        DumpMatrix(name, m, iSpecies, jSpecies, d_pickle)

def LoadMatrix(name, iSpecies, jSpecies, d_pickle): 
    with open(d_pickle + "%s%d_%d.pic" % (name, iSpecies, jSpecies), 'rb') as picFile:  
        M = pic.load(picFile)
    return M
        
def LoadMatrixArray(name, seqsInfo, iSpecies, d_pickle, row=True):
    matrixArray = []
    for jSpecies in range(seqsInfo.nSpecies):
        if row == True:
            matrixArray.append(LoadMatrix(name, iSpecies, jSpecies, d_pickle))
        else:
            matrixArray.append(LoadMatrix(name, jSpecies, iSpecies, d_pickle))
    return matrixArray
              
def MatricesAnd_s(Xarr, Yarr):
    Zarr = []
    for x, y in zip(Xarr, Yarr):
        Zarr.append(x.multiply(y))
    return Zarr
                
def MatricesAndTr_s(Xarr, Yarr):
    Zarr = []
    for x, y in zip(Xarr, Yarr):
        Zarr.append(x.multiply(y.transpose()))
    return Zarr   
    
def DeleteMatrices(baseName, d_pickle):
    for f in glob.glob(d_pickle + baseName + "*_*.pic"):
        if os.path.exists(f): os.remove(f)

def sparse_max_row(csr_mat):
    ret = np.zeros(csr_mat.shape[0])
    ret[np.diff(csr_mat.indptr) != 0] = np.maximum.reduceat(csr_mat.data,csr_mat.indptr[:-1][np.diff(csr_mat.indptr)>0])
    return ret
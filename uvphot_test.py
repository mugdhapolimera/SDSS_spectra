# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 06:49:24 2022

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

uvphotfull = pd.read_csv('RESOLVE_UVids.txt', sep = '\s+', header = None,
                  names = ['name', 'rad_50_g_a',
                'maxuvrad_a', 'NUV_outre', 'K_outre', 'eNUV_outre', 'eK_outre',
                'UVBflag', 'rad_uvsf_a', 'NUV_Xoutre', 'K_Xoutre', 'eNUV_Xoutre',
                'eK_Xoutre', 'limk_Xoutre', 'T1sigcode', 'XUVflag'])
uvphotfull.index = uvphotfull.name

uvphot = uvphotfull[uvphotfull.NUV_outre > 0]
nouv_ruvsf = uvphotfull[uvphotfull.NUV_outre == -777]
uv_within_r50g = uvphotfull[uvphotfull.NUV_outre == -888]
fail = uvphotfull[uvphotfull.NUV_outre == -999]

lowsfoutliers = np.array(['rs0892', 'rs0851', 'rs0749', 'rs0801', 'rs0633', 'rs0608', 
                 'rs0875', 'rs0619', 'rs0685', 'rs0662', 'rs0694', 'rs0632', 
                 'rs0705', 'rs0825', 'rs0031', 'rs0648', 'rs0658', 'rs0714', 
                 'rs0647', 'rs0669', 'rs0682', 'rs0835', 'rs0709', 'rs0653', 
                 'rs1231', 'rs0310', 'rs0511', 'rs0951', 'rs0033', 'rs0205', 
                 'rs0091', 'rs0131', 'rs0383', 'rs0279', 'rs0561', 'rs0399', 
                 'rs0954', 'rs0280', 'rs1288', 'rs1279', 'rs0830', 'rs1292', 
                 'rs0972', 'rs0622', 'rs1201', 'rs0610', 'rs1007', 'rs1097', 
                 'rs1174', 'rf0423', 'rf0211', 'rf0212', 'rf0509', 'rf0220', 
                 'rf0198', 'rf0512', 'rf0498', 'rf0021', 'rf0406', 'rf0492', 
                 'rf0484', 'rf0422', 'rf0080', 'rf0196', 'rf0213', 'rf0184', 
                 'rf0188', 'rf0200', 'rf0215', 'rf0081', 'rf0219', 'rf0597', 
                 'rf0307', 'rf0230', 'rf0139', 'rf0622', 'rf0293', 'rf0114', 
                 'rf0115', 'rf0116', 'rf0082', 'rf0149', 'rf0153', 'rf0161', 
                 'rf0056', 'rf0222', 'rf0281', 'rf0214', 'rf0277', 'rf0288', 
                 'rs1158', 'rs0595', 'rs0757', 'rf0334', 'rf0350', 'rf0345', 
                 'rf0355', 'rf0447', 'rf0382'])

lowsf_xuv = np.intersect1d(lowsfoutliers, uvphot.name[uvphot.XUVflag == 1])
lowsf_nouv = np.intersect1d(lowsfoutliers, nouv_ruvsf.name[nouv_ruvsf.XUVflag == 1])
lowsf_uvrad50 = np.intersect1d(lowsfoutliers, uv_within_r50g.name[uv_within_r50g.XUVflag == 1])
lowsf_fail = np.intersect1d(lowsfoutliers, fail.name[fail.XUVflag == 1])

#highsfoutliers_df = pd.read_csv('highsf_outliers.txt', sep = '\s+')
#highsfoutliers = list(highsfoutliers_df.name)
highsfoutliers = np.array(['rs0053', 'rs0087', 'rs0133', 'rs0135', 'rs0148', 'rs0240', 
                  'rs0243', 'rs0244', 'rs0270', 'rs0461', 'rs0466', 'rs0503',
                  'rs0219', 'rs0531', 'rs0547', 'rs0625', 'rs1337', 'rs0839'])
highsf_xuv = np.intersect1d(highsfoutliers, uvphotfull.name[uvphotfull.XUVflag == 1])
highsf_nouv = np.intersect1d(highsfoutliers, nouv_ruvsf.name)#[nouv_ruvsf.XUVflag == 1])
highsf_uvrad50 = np.intersect1d(highsfoutliers, uv_within_r50g.name)#[uv_within_r50g.XUVflag == 1])
highsf_fail = np.intersect1d(highsfoutliers, fail.name)#[fail.XUVflag == 1])


rescat = readsav('resolvecatalog.dat')
rescatphot= readsav('resolvecatalogphot.dat')
nuvmag = rescatphot.nuvmag
rmag = rescatphot.rmag
good = (nuvmag > 0) & (rmag > 0)
nuvmag = nuvmag[good]
nuv90 = rescatphot.nuv90[good]
rmag = rmag[good]
resname = rescat.name[good]
nuvext = rescatphot.deextrestnuvmag[good] - rescatphot.smoothrestnuvmag[good]
lowoutres, lowndx, reslowndx = np.intersect1d(lowsfoutliers, resname, return_indices = True)
highoutres, highndx, reshighndx = np.intersect1d(highsfoutliers, resname, return_indices = True)

plt.figure()
plt.scatter(rmag, nuvmag)
plt.scatter(rmag[reslowndx], nuvmag[reslowndx], c= 'k')
plt.scatter(rmag[reshighndx], nuvmag[reshighndx], c= 'r')
plt.xlim(21,11)
plt.ylim(24,14)
plt.xlabel('RESOLVE r mag')
plt.ylabel('RESOLVE NUV mag')

plt.figure()
plt.scatter(nuvmag, nuvmag - nuv90)
plt.scatter(nuvmag[reslowndx], nuvmag[reslowndx] - nuv90[reslowndx], c= 'k')
plt.scatter(nuvmag[reshighndx], nuvmag[reshighndx] - nuv90[reshighndx], c= 'r')

plt.figure()
plt.scatter(nuvmag, nuvext)
plt.scatter(nuvmag[reslowndx], nuvext[reslowndx], c= 'k')
plt.scatter(nuvmag[reshighndx], nuvext[reshighndx], c= 'r')

agn = pd.read_csv('ECO+RESOLVE_AGN_list.csv')
agn.index = agn.name
highsfout_agn = np.intersect1d(highsfoutliers, agn.name)#[fail.XUVflag == 1])


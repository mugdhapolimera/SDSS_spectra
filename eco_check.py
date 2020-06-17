# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 12:24:12 2020

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import scipy.stats.kde as kde
import matplotlib.pyplot as plt
import statsmodels.api as sm

eco = pd.read_csv('ECO_snr5_master_new.csv')
res = pd.read_csv('RESOLVE_snr5_master_new.csv')

ecoflags = pd.read_csv('eco_emlineclass_full_snr5_master_new.csv')
resflags = pd.read_csv('resolve_emlineclass_full_snr5_master_new.csv')

ecodwarf = (eco.logmstar < 9.5) & ((eco.logmgas - eco.logmstar) > 0)
resdwarf = (res.logmstar < 9.5) & ((res.logmgas - res.logmstar) > 0)

ptype = 'kde'
def ppty_kde(arr, label, xlabel = ' ',plttype = 'kde', **kwargs):
    arr = np.array(arr, '<d')
    hist = np.histogram(arr, bins = 'fd', normed = True)
    #n, bins, patches = plt.hist(arr, bins = 'fd', normed = True, alpha = 0.3)
    try:
        bin_width = kwargs.pop('bin_width')
    except KeyError:
        bin_width = np.diff(hist[1])[0]
    xaxis = np.arange(min(arr)-0.5, max(arr)+0.5, 0.1)
#    xaxis = np.arange(min(arr), max(arr), 0.1)
    pptykde = sm.nonparametric.KDEUnivariate(arr)
    try:                                                                  
        pptykde.fit(bw = bin_width)
    except ValueError:
        arr = np.array(arr,dtype = np.double)
        newarr = arr.byteswap().newbyteorder()
        pptykde = sm.nonparametric.KDEUnivariate(newarr)
        pptykde.fit(bw = bin_width)
    try: 
        lw = kwargs.pop('lw')
        if plttype == 'kde':
            plt.plot(pptykde.support, pptykde.density, label = label, lw = lw)
        elif plttype == 'hist':
            plt.hist(arr, bins = 'fd', histtype = 'step', lw = lw, normed = True,
                     label = label)
            #plt.plot(hist[1][1:],hist[0]/float(len(arr)),label=label, lw= lw)
    except KeyError:
        if plttype == 'kde':
            plt.plot(pptykde.support, pptykde.density, label = label)
        elif plttype == 'hist':
            plt.hist(arr, bins = 'fd', histtype = 'step', normed = True, 
                     label = label)
        #plt.plot(hist[1][1:],hist[0]/float(len(arr)),label=label)
#    kernel = kde.gaussian_kde(arr, bw_method = bin_width/np.std(arr))
#    pptykde = kernel(xaxis)
#    plt.plot(xaxis, pptykde, label = label)
    #, color = patches[0].get_facecolor())    
    plt.xlabel(xlabel)
    return (xaxis, pptykde)

#Distribution of Stellar Masses
plt.figure()
mstarxaxis, mstarkde = ppty_kde(eco.logmstar, xlabel = 'Stellar Mass',plttype = ptype, 
                                label = 'ECO SELs')
mstarxaxis, mstarkde = ppty_kde(res.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE SELs')
mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecodwarf], xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'ECO Gas Rich SEL Dwarfs ')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resdwarf], xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE Gas-Rich SEL Dwarfs')
mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecoflags.sftoagn], lw = 3,plttype = ptype,
                                xlabel = 'Stellar Mass', label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resflags.sftoagn], lw = 3,plttype = ptype,
                                xlabel = 'Stellar Mass', label = 'RESOLVE SFing-AGN')

plt.legend()  


plt.figure()
ecombary = np.log10(10**eco.logmstar + 10**eco.logmgas)
resmbary = np.log10(10**res.logmstar + 10**res.logmgas)
mbaryxaxis, mbarykde = ppty_kde(ecombary, xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'ECO SELs')
mbaryxaxis, mbarykde = ppty_kde(resmbary, label = 'RESOLVE SELs',plttype = ptype)
mbaryxaxis, mbarykde = ppty_kde(ecombary[ecoflags.sftoagn], xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'ECO SFing-AGN')
mbaryxaxis, mbarykde = ppty_kde(resmbary[resflags.sftoagn], xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'RESOLVE SFing-AGN')
plt.legend()


ra = res.radeg
flinsample = res.fl_insample
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
grpcz = res.grpcz
fallthreshold = 9.2
mgas = res.logmgas
mstars = res.logmstar
mbary = 10**mgas + 10**mstars
springinobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
            (((flinsample | (np.log10(mbary) > 9.2)) & inspring))

fallinobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                   ((flinsample | (np.log10(mbary) > fallthreshold)) & infall)) 

inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
        (((flinsample | (np.log10(mbary) > fallthreshold)) & infall) | \
         ((flinsample | (np.log10(mbary) > 9.2)) & inspring))

res = res[inobssample]
plt.figure()
ecombary = np.log10(10**eco.logmstar + 10**eco.logmgas)
resmbary = np.log10(10**res.logmstar + 10**res.logmgas)
mbaryxaxis, mbarykde = ppty_kde(ecombary, xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'ECO SELs')
mbaryxaxis, mbarykde = ppty_kde(resmbary, label = 'RESOLVE SELs')
mbaryxaxis, mbarykde = ppty_kde(ecombary[ecoflags.sftoagn], xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'ECO SFing-AGN')
mbaryxaxis, mbarykde = ppty_kde(resmbary[resflags.sftoagn], xlabel = 'Baryonic Mass',plttype = ptype,
                                label = 'RESOLVE SFing-AGN')
plt.legend()

plt.figure()
mstarxaxis, mstarkde = ppty_kde(eco.logmstar, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO SELs')
mstarxaxis, mstarkde = ppty_kde(res.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE SELs')
mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecodwarf], xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO Gas Rich SEL Dwarfs ')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resdwarf], xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE Gas-Rich SEL Dwarfs')
mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecoflags.sftoagn], plttype = ptype,
                                xlabel = 'Stellar Mass', label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resflags.sftoagn], plttype = ptype,
                                xlabel = 'Stellar Mass', label = 'RESOLVE SFing-AGN')

plt.legend()  

from scipy.io.idl import readsav
ecoint = readsav('eco_wresa_032918.dat')
resintphot = readsav('resolvecatalogphot.dat')
resint = readsav('resolvecatalog.dat')
ecofull = pd.read_csv('ECO_live22Oct2018.csv')

plt.figure()
endx = [x for x in range(len(ecoint.econames)) if ecoint.econames[x] in 
        list(eco.name)]
rndx = [x for x in range(len(resint.name)) if resint.name[x] in 
        list(res.name)]
mstarxaxis, mstarkde = ppty_kde(ecoint.mur90[endx], plttype = ptype,label = 'ECO SELs')
mstarxaxis, mstarkde = ppty_kde(resintphot.mur90[rndx], plttype = ptype,label = 'RESOLVE SELs')

ecosftoagn = list(ecoflags.galname[ecoflags.sftoagn])
ressftoagn = list(resflags.galname[resflags.sftoagn])

econdx = [x for x in range(len(ecoint.econames)) if ecoint.econames[x] in ecosftoagn]
resndx = [x for x in range(len(resint.name)) if resint.name[x] in ressftoagn]

ecodwarfndx = [x for x in range(len(ecoint.econames)) if ecoint.econames[x] \
                in list(eco.name[ecodwarf])]
resdwarfndx = [x for x in range(len(resint.name)) if resint.name[x] in \
               list(res.name[resdwarf])]

mstarxaxis, mstarkde = ppty_kde(ecoint.mur90[econdx], lw = 3,plttype = ptype,
                                label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(resintphot.mur90[resndx], lw = 3,plttype = ptype,
                                xlabel = 'mu_r90', label = 'RESOLVE SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(ecoint.mur90[ecodwarfndx], xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO Gas Rich SEL Dwarfs ')
mstarxaxis, mstarkde = ppty_kde(resintphot.mur90[resdwarfndx], plttype = ptype,
                                xlabel = 'mu_r90 (r-band surface brightness within 90% light radius)',
                                label = 'RESOLVE Gas-Rich SEL Dwarfs')

plt.legend()  

plt.figure()
mstarxaxis, mstarkde = ppty_kde(eco.grpn, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO SELs', bin_width = 1)
mstarxaxis, mstarkde = ppty_kde(res.grpn, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE SELs')
mstarxaxis, mstarkde = ppty_kde(eco.grpn[eco.grpn < 200][ecoflags.sftoagn], plttype = ptype,
                                xlabel = 'Stellar Mass', label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(res.grpn[resflags.sftoagn], plttype = ptype,
                                xlabel = 'Group N', label = 'RESOLVE SFing-AGN')

plt.legend()  

###############################################################################
#Plot ECO and RESOLVE distrubutions for every step
###############################################################################
resfull = pd.read_csv('RESOLVE_full_blend_dext_new.csv')
ptype = 'hist'
plt.figure()
#mstarxaxis, mstarkde = ppty_kde(ecofull.logmstar, xlabel = 'Stellar Mass', 
#                                label = 'ECO Parent Sample')
#mstarxaxis, mstarkde = ppty_kde(resfull.logmstar, xlabel = 'Stellar Mass', 
#                                label = 'RESOLVE Parent Sample')
#plt.legend()
plt.xlim(6.9,12.1)
plt.ylim(-0.0001,1.3)
ecoobs = pd.read_csv('ECO_inobssample.csv')
resobs = pd.read_csv('RESOLVE_inobssample.csv')

mstarxaxis, mstarkde = ppty_kde(ecoobs.logmstar, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO inobssample')
mstarxaxis, mstarkde = ppty_kde(resobs.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE inobssample')

plt.legend()

ecojhu = pd.read_csv('ECO_full_snr5.csv')
resjhu = pd.read_csv('RESOLVE_full_snr5.csv')

mstarxaxis, mstarkde = ppty_kde(ecojhu.logmstar, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO JHU SELs')
mstarxaxis, mstarkde = ppty_kde(resjhu.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE JHU SELs')
plt.legend()

ecoport = pd.read_csv('ECO_full_snr5_port.csv')
resport = pd.read_csv('RESOLVE_full_snr5_port.csv')
mstarxaxis, mstarkde = ppty_kde(ecoport.logmstar, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO Portsmouth SELs')
mstarxaxis, mstarkde = ppty_kde(resport.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE Portsmouth SELs')
plt.legend()


mstarxaxis, mstarkde = ppty_kde(eco.logmstar, xlabel = 'Stellar Mass', plttype = ptype,
                                label = 'ECO SELs')
mstarxaxis, mstarkde = ppty_kde(res.logmstar, xlabel = 'Stellar Mass',plttype = ptype,
                                label = 'RESOLVE SELs')
plt.legend()

mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecodwarf], xlabel = 'Stellar Mass', plttype = ptype,
                                lw = 2, label = 'ECO Gas Rich SEL Dwarfs ')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resdwarf], xlabel = 'Stellar Mass',plttype = ptype,
                                lw = 2, label = 'RESOLVE Gas-Rich SEL Dwarfs')
plt.legend()
mstarxaxis, mstarkde = ppty_kde(eco.logmstar[ecoflags.sftoagn], plttype = ptype,
                                lw = 4, xlabel = 'Stellar Mass', label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(res.logmstar[resflags.sftoagn], plttype = ptype,
                                lw = 4, xlabel = 'Stellar Mass', label = 'RESOLVE SFing-AGN')

plt.legend()  


##############################################################################
ptype = 'kde'
plt.figure()
mstarxaxis, mstarkde = ppty_kde(eco.logmh, xlabel = 'Halo Mass', plttype = ptype,
                                label = 'ECO SELs')
mstarxaxis, mstarkde = ppty_kde(res.logmh, xlabel = 'Halo Mass',plttype = ptype,
                                label = 'RESOLVE SELs')
plt.legend()

mstarxaxis, mstarkde = ppty_kde(eco.logmh[ecodwarf], xlabel = 'Halo Mass', plttype = ptype,
                                lw = 2, label = 'ECO SEL Dwarfs ')
mstarxaxis, mstarkde = ppty_kde(res.logmh[resdwarf], xlabel = 'Stellar Mass',plttype = ptype,
                                lw = 2, label = 'RESOLVE SEL Dwarfs')
plt.legend()
mstarxaxis, mstarkde = ppty_kde(eco.logmh[ecoflags.sftoagn], plttype = ptype,
                                lw = 4, xlabel = 'Halo Mass', label = 'ECO SFing-AGN')
mstarxaxis, mstarkde = ppty_kde(res.logmh[resflags.sftoagn], plttype = ptype,
                                lw = 4, xlabel = 'Halo Mass', label = 'RESOLVE SFing-AGN')

plt.legend()  


import matplotlib.colors as colors
import matplotlib.cm as cm
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
sf_colors_map = truncate_colormap(cm.gray_r)

fig, ax1 = plt.subplots()
definite = np.column_stack((eco.logmh, eco.logmstar))
xmin = 10; xmax = 15
nbins = 100
ymin = 8; ymax = 11.5
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((res.logmh, res.logmstar))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#               shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
res_plot = ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], colors='k',
            label = 'RESOLVE')
res_plot.collections[0].set_label('RESOLVE')
eco_plot = ax1.contourf(xgrid, ygrid, definite_z.reshape(xgrid.shape),
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
            label = 'ECO')
eco_plot.collections[0].set_label('ECO')
plt.legend()
ax1.set_xlabel('log(Halo Mass)')
ax1.set_ylabel('log(Stellar Mass)')

fig, ax1 = plt.subplots()
definite = np.column_stack((eco.logmh[ecodwarf], eco.logmstar[ecodwarf]))
xmin = 10; xmax = 15
nbins = 100
ymin = 8; ymax = 11.5
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((res.logmh[resdwarf], res.logmstar[resdwarf]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#               shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
res_plot = ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
            [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], colors='k',
            label = 'RESOLVE')
res_plot.collections[0].set_label('RESOLVE Dwarfs')
eco_plot = ax1.contourf(xgrid, ygrid, definite_z.reshape(xgrid.shape),
            [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
            label = 'ECO')
eco_plot.collections[0].set_label('ECO Dwarfs')
plt.legend()
ax1.set_xlabel('log(Halo Mass)')
ax1.set_ylabel('log(Stellar Mass)')

fig, ax1 = plt.subplots()
definite = np.column_stack((eco.logmh[ecoflags.sftoagn], 
                            eco.logmstar[ecoflags.sftoagn]))
xmin = 10; xmax = 15
nbins = 100
ymin = 8; ymax = 11.5
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((res.logmh[resflags.sftoagn], 
                               res.logmstar[resflags.sftoagn]))
#xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
#ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xmin_agn = xmin; ymin_agn = ymin
xmax_agn = xmax; ymax_agn = ymax
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#               shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
res_plot = ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
            3, colors='k',
            label = 'RESOLVE')
#ax1.clabel(res_plot, inline=True, fontsize=8)
res_plot.collections[0].set_label('RESOLVE SFing-AGN')
eco_plot = ax1.contourf(xgrid, ygrid, definite_z.reshape(xgrid.shape),
            3, label = 'ECO')
#ax1.clabel(eco_plot, inline=True, fontsize = 8)
eco_plot.collections[0].set_label('ECO SFing-AGN')
plt.legend()
ax1.set_xlabel('log(Halo Mass)')
ax1.set_ylabel('log(Stellar Mass)')


#fig,ax2 = plt.subplots()
#ax2.tricontour(definite[:,0], definite[:,1], definite_z.reshape(xgrid.shape), 
#               levels=14, linewidths=0.5, colors='k')
#cntr2 = ax2.tricontourf(definite, definite_z.reshape(xgrid.shape), 
#                        levels=14, cmap="RdBu_r")
#
#fig.colorbar(cntr2, ax=ax2)
#ax2.plot(definite, 'ko', ms=3)
##ax2.set(xlim=(-2, 2), ylim=(-2, 2))
##ax2.set_title('tricontour (%d points)' % npts)
#
#plt.subplots_adjust(hspace=0.5)
#plt.show()
#
#from astroML.plotting import scatter_contour as sc
#fig,ax2 = plt.subplots()
#sc(eco.logmh, eco.logmstar,  
#   levels = 10,threshold = 20)
#sc(res.logmh, res.logmstar, 
#   levels = 10,threshold = 20,
#   plot_args=dict(color='black'),
#                contour_args=dict(cmap=plt.cm.bone))
#
#ax2.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
#
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = plt.axes(projection='3d')
##ax.plot_surface(gridx1, gridy1, flux_idl[k], color = 'b')
#eco_surf = ax.plot_surface(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                cmap = 'coolwarm', edgecolor = 'none')#color = 'b')
#res_surf = ax.plot_surface(xgrid_agn, ygrid_agn, 
#                agn_contour_z.reshape(xgrid_agn.shape), 
#                cmap = 'coolwarm', edgecolor = 'none')#color = 'r')
##ax.scatter(test_input[:,0],test_input[:,1],test_output,c = 'r')
##ax.plot_surface(gridx, gridy, fluxarr[no], color = 'g')
#ax.set(xlabel="Log(Halo Mass)", ylabel="Log(Stellar Mass)", zlabel="Density")
#fig.colorbar(eco_surf, shrink=0.5, aspect=5)
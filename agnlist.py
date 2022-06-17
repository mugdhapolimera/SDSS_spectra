# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:10:59 2020

@author: mugdhapolimera

Compiling all AGN in RESOLVE
"""
import numpy as np
import pandas as pd
from scipy.io.idl import readsav
from astropy.coordinates import SkyCoord
import astropy.units as u

res = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_inobssample.csv')
res.index = res.name

ir = pd.read_csv("mid_ir/RESOLVE_WISE_AGN.csv")
opt = pd.read_csv("ECO/SEL/eco+resolve_emlineclass_dext_snr5_jhu.csv")
optport = pd.read_csv("eco+resolve_emlineclass_dext_snr5_port.csv")
optnsa = pd.read_csv("eco+resolve_emlineclass_dext_snr5_nsa.csv")
bpt = pd.read_csv("eco+resolve_emlineclass_full_bpt1snr5_jhu.csv")
#opt = pd.read_csv("ECO/SEL/eco_emlineclass_dext_snr5_jhu.csv")
#bpt = pd.read_csv("eco_emlineclass_full_bpt1snr5_jhu.csv")
xray = pd.read_csv("../xray/ECO+RESOLVE_inobssample_xray_chandra_new.csv")
xray_xmm = pd.read_csv("../xray/ECO+RESOLVE_inobssample_xray_xmm.csv")
xray = xray.append(xray_xmm)
iragn = np.array(ir.name)
bptagn = np.array(bpt.galname[bpt.defagn])
bptcomposite = np.array(bpt.galname[bpt.composite])
xrayagn = np.unique(xray.name[xray.xrayagn])
optagn = np.array(opt.galname[opt.defagn])
optcomposite = np.array(opt.galname[opt.composite])
optsfingagn = np.array(opt.galname[opt.sftoagn])
optagntosf = np.array(opt.galname[opt.agntosf])

portnotinjhu = np.setdiff1d(optport.galname[~optport.defstarform],opt.galname)
optport.index = optport.galname
optport = optport.loc[portnotinjhu]
optagnport = np.array(optport.galname[optport.defagn])
optcompositeport = np.array(optport.galname[optport.composite])
optsfingagnport = np.array(optport.galname[optport.sftoagn])
optagntosfport = np.array(optport.galname[optport.agntosf])

nsanotinjhu = np.setdiff1d(optnsa.galname[~optnsa.defstarform],opt.galname)
optnsa.index = optnsa.galname
optnsa = optnsa.loc[nsanotinjhu]
optagnnsa = np.array(optnsa.galname[optnsa.defagn])
optcompositensa = np.array(optnsa.galname[optnsa.composite])
optsfingagnnsa = np.array(optnsa.galname[optnsa.sftoagn])
optagntosfnsa = np.array(optnsa.galname[optnsa.agntosf])

optagn = np.unique(np.concatenate((optagn, optagnport, optagnnsa)))
optcomposite = np.unique(np.concatenate((optcomposite, optcompositeport, optcompositensa)))
optsfingagn = np.unique(np.concatenate((optsfingagn, optsfingagnport, optsfingagnnsa)))
optagntosf = np.unique(np.concatenate((optagntosf, optagntosfport, optagntosfnsa)))

bptagn = np.array(list(set(bptagn) | set(optagn)))
bptcomposite = np.array(list(set(bptcomposite) | set(optcomposite)))

allagn = np.unique(list(set(iragn) | set(bptagn) | set(bptcomposite) | set(xrayagn) | 
        set(optsfingagn) | set(optagntosf)))

irflag = [x in iragn for x in allagn]
xrayflag = [x in xrayagn for x in allagn]
bptagnflag = [x in bptagn for x in allagn]
bptcompositeflag = [x in bptcomposite for x in allagn]
sfingagnflag = [x in optsfingagn for x in allagn]
agntosfflag = [x in optagntosf for x in allagn]

agntype = pd.DataFrame({'name':allagn,
                    'midiragn':irflag,
                    'xrayagn':xrayflag,
                    'bptagn':bptagnflag,
                    'bptcomposite':bptcompositeflag,
                    'sfingagn':sfingagnflag,
                    'agntosf':agntosfflag,
                    })
agntype.index = agntype.name
#agn.to_csv("ECO+RESOLVE_AGN_list.csv")
allagn = np.intersect1d(allagn, list(res.name))
###############################################################################
cols = ['name','radeg','dedeg','vhel','logmstar','logmgas','logmh',
        'gmag','rmag','oiii_5007_flux','o3lum','nuvmag','rband_censb',
        'obstel','obsstat','broadlineagn','agntype','ifusource']
resagn = res.loc[allagn]
#agn = pd.DataFrame({'name':np.array(resagn.name),
#                    'ra':np.array(resagn.radeg),
#                    'dec':np.array(resagn.dedeg),
#                    'logmstar':np.array(resagn.logmstar),
#                    'logmgas':np.array(resagn.logmgas),
#                    'logmh':np.array(resagn.logmh),
#                    'gmag':np.array(resagn.logmgas),
#                    'rmag':np.array(resagn.logmgas),
#                    'oiii':np.array(resagn.logmgas),
#                    })

#get vhel for resagn
internal = readsav('resolvecatalog.dat')
internalphot = readsav('resolvecatalogphot.dat')
names, resndx, internalndx = np.intersect1d(resagn.name, internal.name, 
                                            return_indices = True)
resagn['vhel'] = np.zeros(len(resagn))
resagn['vhel'].loc[resagn.name[resndx]] = internal.vhel[internalndx]

ecointernal = readsav('eco_wresa_032918.dat')
names, resndxforeco, ecointernalndx = np.intersect1d(resagn.name, ecointernal.econames, 
                                            return_indices = True)
cf = np.pi/180.
econame = np.asarray(ecointernal.econames)
ecora = np.asarray(ecointernal.goodra)
ecodec = np.asarray(ecointernal.gooddec)
ecoeqcoords = SkyCoord(ra=ecora*u.degree, dec=ecodec*u.degree)
ecogalactic = ecoeqcoords.galactic
lgcorr = 300.*np.cos(ecogalactic.b.value*cf)*np.sin(ecogalactic.l.value*cf)
ecovhel = (3e5*ecointernal.goodz) - lgcorr
resagn['vhel'].loc[resagn.name[resndxforeco]] = ecovhel[ecointernalndx]

#o3lum and flux
jhu = pd.read_csv('ECO+RESOLVE_snr5_dext_jhu.csv')
jhu.index = jhu.name
port = pd.read_csv('ECO+RESOLVE_snr5_dext_port.csv')
port.index = port.name
nsa = pd.read_csv('ECO+RESOLVE_snr5_dext_nsa.csv')
nsa.index = nsa.name
resagn['oiii_5007_flux'] = [np.nan]*len(resagn)
resagn['oiii_5007_flux'] = jhu.loc[resagn.name]['oiii_5007_flux']
resagn['oiii_5007_flux_err'] = [np.nan]*len(resagn)
resagn['oiii_5007_flux_err'] = jhu.loc[resagn.name]['oiii_5007_flux_err']

resagn.loc[portnotinjhu]['oiii_5007_flux'] = port.loc[portnotinjhu]['oiii_5007_flux']
resagn.loc[portnotinjhu]['oiii_5007_flux_err'] = port.loc[portnotinjhu]['oiii_5007_flux_err']

resagn.loc[portnotinjhu]['oiii_5007_flux'] = nsa.loc[nsanotinjhu]['oiii_5007_flux']
resagn.loc[portnotinjhu]['oiii_5007_flux_err'] = nsa.loc[nsanotinjhu]['oiii_5007_flux_err']

o3flux = resagn['oiii_5007_flux'] * 1e-17 #ergs/s
vel = resagn['vhel']
H0 = 70 #km/s/Mpc
Mpc = 3.086e24 #cm
d = (vel/H0) * Mpc #cm
resagn['o3lum'] = 4*3.14*(d**2) * o3flux

#rband_censb
resagn['rband_censb'] = np.zeros(len(resagn))
resagn['rband_censb'].loc[resagn.name[resndx]] = internal.ifusb[internalndx]

#obstel
telarr = readsav('telarray.dat')['telarr']
resagn['obstel'] = ['N/A']*len(resagn)
resagn['obstel'].loc[resagn.name[resndx]] = telarr[internalndx]
for x in list(resagn.name): 
    if resagn.obstel.loc[x] == 'KINDONE':
        ndx = np.where(internal.name == x)
        if internal.blue[ndx]:
           resagn.obstel.loc[x] = 'SOARABS'     
        elif internal.red[ndx]:
           resagn.obstel.loc[x] = 'SOAREM'     
        elif internal.gemblue[ndx]:
           resagn.obstel.loc[x] = 'GEMABS'     
        elif internal.gemred[ndx]:
           resagn.obstel.loc[x] = 'GEMEM'     
        elif internal.saltlsblue[ndx]:
           resagn.obstel.loc[x] = 'SALTblue'     
        elif internal.saltls[ndx]:
           resagn.obstel.loc[x] = 'SALT'     

#obsstat
resagn['obsstat'] = ['notobs']*len(resagn)
blueflag = ((internal.blue) | (internal.gemblue) | \
        (internal.koalablue) | (internal.saltlsblue))
redflag = ((internal.red) | (internal.gemred) | \
        (internal.koalared) | (internal.saltls))
broadflag = internal.broad
for x in list(resagn.name): 
        ndx = np.where(internal.name == x)
        if blueflag[ndx] | redflag[ndx] | broadflag[ndx]:
            obsstr = ''
            if blueflag[ndx]:
               obsstr+='bluedone/'     
            if redflag[ndx]:
               obsstr+='reddone/'     
            if broadflag[ndx]:
               obsstr+='broaddone/'     
            resagn.obsstat.loc[x] = obsstr

#broadlineagn
broadline = pd.read_csv('ECO_BroadlineAGN_Liu19.csv')
broadline = broadline.append(pd.read_csv('RESOLVE_BroadlineAGN_Liu19.csv'))
broadline.index = broadline['name_1']
resagn['broadlineagn'] = [False]*len(resagn)
broadnames = np.intersect1d(resagn.name, broadline['name_1'])
resagn['broadlineagn'].loc[broadnames] = True

#agntype
resagn['agntype'] = np.zeros(len(resagn))
for key in agntype.keys():
    if key != 'name':
        resagn['agntype'].loc[agntype[key]] = key

#ifusource
manga = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\catalog_matching\optical\MANGA_ECO+RESOLVE_inobssample.csv')
manga.index = manga.name_1
sami = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\catalog_matching\optical\SAMI_ECO+RESOLVE_inobssample.csv')
sami.index = sami.name_1
resagn['ifusource'] = ['none']*len(resagn)
manganames= np.intersect1d(manga.name_1, resagn.name)
saminames= np.intersect1d(sami.name_1, resagn.name)
resagn['ifusource'].loc[manganames] = 'manga'
resagn['ifusource'].loc[saminames] = 'sami'

agnlist = resagn[cols]
agnlist.to_csv("ECO+RESOLVE_AGN_list.csv")


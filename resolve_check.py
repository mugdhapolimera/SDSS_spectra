# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:22:43 2020o

@author: mugdhapolimera

Counting numbers of RESOLVE-A and B objects to quantify the difference
between the two parts of the survey

Questions:
    1. Why does RESOLVE-A (spring) have more AGN than RESOLVE-B?
        a. Environment dependent?
        b. SDSS coverage dependent?
        c. Related to more low-sb galaxies in RESOLVE-B?
        d. Just chance?
    2. Does ECO follow the same pattern as RESOLVE-A, or is there further
       difference due to environment, completeness etc.?
       
"""
import numpy as np
import pandas as pd
from astropy.stats import binom_conf_interval
from IPython.display import Latex, Math
###############################################################################
# Full RESOLVE survey, including buffers
###############################################################################

res = pd.read_csv('RESOLVE_full_blend_dext_new.csv')
res.index = res.name

fallthreshold = 9.0
fullinobssample = 1
baryinobssample = 0
spring = [x for x in list(res.name) if x[1] == 's']
fall = [x for x in list(res.name) if x[1] == 'f']

totalres = len(res)
totalresspring = len(spring)
totalresfall = len(fall)
pcresspring = 100.0*totalresspring/totalres
pcresfall = 100.0*totalresfall/totalres
print('Full RESOLVE Survey (including buffers)')
print('Total: {} \nSpring: {} ({:.2f}%) \nFall: {} ({:.2f}%)'
      .format(totalres,totalresspring,pcresspring,totalresfall,pcresfall))

###############################################################################
# SDSS matches of full RESOLVE survey (including buffers)
###############################################################################

sdss = pd.read_csv('RESOLVE_full_cont.csv')
sdss.index = sdss.name

sdssspring = [x for x in list(sdss.name) if x[1] == 's']
sdssfall = [x for x in list(sdss.name) if x[1] == 'f']

totalsdss = len(sdss)
totalsdssspring = len(sdssspring)
totalsdssfall = len(sdssfall)
pcsdssspring = 100.0*totalsdssspring/totalsdss
pcsdssfall = 100.0*totalsdssfall/totalsdss

print('\n\nSDSS matches of full RESOLVE survey (including buffers)')
print('Total: {} \nSpring: {} ({:.2f}%) \nFall: {} ({:.2f}%)'
      .format(totalsdss,totalsdssspring,pcsdssspring,totalsdssfall,pcsdssfall))

###############################################################################
# RESOLVE + inobssample
###############################################################################
ra=res.radeg
dec=res.dedeg
flinsample = res.fl_insample
grpcz = res.grpcz
cz = res.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = res.logmgas
mstars = res.logmstar
mbary = 10**mgas + 10**mstars
if fullinobssample: 
    springinobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                (((flinsample | (np.log10(mbary) > 9.2)) & inspring))
    
    fallinobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                       ((flinsample | (np.log10(mbary) > fallthreshold)) & infall)) 
    
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
            (((flinsample | (np.log10(mbary) > fallthreshold)) & infall) | \
             ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
if baryinobssample : 
    springinobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                ((((np.log10(mbary) > 9.2)) & inspring))
    
    fallinobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                       (((np.log10(mbary) > fallthreshold )) & infall)) 
    
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
            ((((np.log10(mbary) > fallthreshold)) & infall) | \
             (((np.log10(mbary) > 9.2)) & inspring))

totalresinobs = np.sum(inobssample)
totalresspringinobs = np.sum(springinobssample)
totalresfallinobs = np.sum(fallinobssample)
pcresspringinobs = 100.0*totalresspringinobs/totalresinobs
pcresfallinobs = 100.0*totalresfallinobs/totalresinobs
print(np.sum(((fallinobssample) & (mstars < 9.5))))
print('\n\ninobssample in full RESOLVE survey')
print('Total: {} \nSpring: {} ({:.2f}%) \nFall: {} ({:.2f}%)'
      .format(totalresinobs, totalresspringinobs, pcresspringinobs,
              totalresfallinobs, pcresfallinobs))

###############################################################################
# RESOLVE in SDSS + inobssample
###############################################################################
ra=res.radeg.loc[sdss.name]
dec=res.dedeg.loc[sdss.name]
flinsample = res.fl_insample.loc[sdss.name]
grpcz = res.grpcz.loc[sdss.name]
cz = res.cz.loc[sdss.name]
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = res.logmgas.loc[sdss.name]
mstars = res.logmstar.loc[sdss.name]
mbary = 10**mgas + 10**mstars
if fullinobssample: 
    springinobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                (((flinsample | (np.log10(mbary) > 9.2)) & inspring))
    
    fallinobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                       ((flinsample | (np.log10(mbary) > fallthreshold)) & infall)) 
    
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
            (((flinsample | (np.log10(mbary) > fallthreshold)) & infall) | \
             ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
if baryinobssample : 
    springinobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                ((((np.log10(mbary) > 9.2)) & inspring))
    
    fallinobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
                       (((np.log10(mbary) > fallthreshold )) & infall)) 
    
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
            ((((np.log10(mbary) > fallthreshold)) & infall) | \
             (((np.log10(mbary) > 9.2)) & inspring))

totalsdssinobs = np.sum(inobssample)
totalsdssspringinobs = np.sum(springinobssample)
totalsdssfallinobs = np.sum(fallinobssample)

pcsdssspringinobs = 100.0*totalsdssspringinobs/totalsdssinobs
pcsdssfallinobs = 100.0*totalsdssfallinobs/totalsdssinobs

print('\n\ninobssample of SDSS matches of RESOLVE survey')
print('Total: {} \nSpring: {} ({:.2f}%) \nFall: {} ({:.2f}%)'
      .format(totalsdssinobs, totalsdssspringinobs,pcsdssspringinobs, 
              totalsdssfallinobs, pcsdssfallinobs))

###############################################################################
# RESOLVE in SDSS & inobssample & SNR > 5 & 10**-3 < flux < 1e5 in ALL 3 CATALOGS
###############################################################################
master = pd.read_csv('RESOLVE_snr5_master_new.csv')
master.index = master.name
masterspring = [x for x in list(master.name) if x in list(sdss[springinobssample].index)]
masterfall = [x for x in list(master.name) if x in list(sdss[fallinobssample].index)]

#masterspring = [x for x in list(master.name) if x[1] == 's']
#masterfall = [x for x in list(master.name) if x[1] == 'f']

totalmaster = len(master)
totalmasterspring = len(masterspring)
totalmasterfall = len(masterfall)
pcmasterspring = 100.0*totalmasterspring/totalmaster
pcmasterfall = 100.0*totalmasterfall/totalmaster
print('\n\nRESOLVE in SDSS & inobssample & SNR > 5 & 10**-3 < flux < 1e5 in ALL 3 CATALOGS')
print('Total: {} \nSpring: {} ({:.2f}%) \nFall: {} ({:.2f}%)'
      .format(totalmaster,totalmasterspring,pcmasterspring,
              totalmasterfall,pcmasterfall))

###############################################################################
# AGN in RESOLVE Master Catalog
###############################################################################
conf = pd.read_csv('RESOLVE_snr5_master_conf_new.csv')
conf.index = conf.name

agn = conf.confidence_level > 0.0
#confspring = [x for x in list(conf.name) if x[1] == 's']
#conffall = [x for x in list(conf.name) if x[1] == 'f']
confspring = [x for x in list(conf.name) if x in list(sdss[springinobssample].index)]
conffall = [x for x in list(conf.name) if x in list(sdss[fallinobssample].index)]

totalconf = len(conf)
confspringagn = len(conf.loc[confspring][agn])
conffallagn = len(conf.loc[conffall][agn])

pcconfspringagn = 100.0*confspringagn/totalconf
pcconffallagn  = 100.0*conffallagn/totalconf
print('\n\nRESOLVE Master Catalog All AGN from Confidence Levels')
print('Total: {} \nSpring AGN: {} ({:.2f}%) \nFall AGN: {} ({:.2f}%)'
      .format(totalconf,confspringagn,pcconfspringagn,
              conffallagn,pcconffallagn))
print('\nSpring AGN: {}/{} = {:.2f}% of spring sample'
      .format(confspringagn,len(confspring), 100.0*confspringagn/len(confspring)))
print('\nFall AGN: {}/{} = {:.2f}% of fall sample'
      .format(conffallagn,len(conffall), 100.0*conffallagn/len(conffall)))

###############################################################################
# Dwarf AGN in RESOLVE Master Catalog
###############################################################################
dwarf = master[master.logmstar < 9.5]
dwarfspring = [x for x in list(dwarf.name) if x[1] == 's']
dwarffall = [x for x in list(dwarf.name) if x[1] == 'f']

#dwarfspring = [x for x in list(dwarf.name) if x in list(sdss[springinobssample].index)]
#dwarffall = [x for x in list(dwarf.name) if x in list(sdss[fallinobssample].index)]

totaldwarf = len(dwarf)

pcdwarfspring = 100.0*len(dwarfspring)/totalmasterspring
pcdwarffall = 100.0*len(dwarffall)/totalmasterfall

dwarfagn = conf.loc[list(dwarf.name)].confidence_level > 0.0

dwarfspringagn = len(dwarf.loc[dwarfspring][dwarfagn])
dwarffallagn = len(dwarf.loc[dwarffall][dwarfagn])
totaldwarfagn = dwarfspringagn+dwarffallagn

pcdwarfspringagn = 100.0*dwarfspringagn/totaldwarfagn
pcdwarffallagn  = 100.0*dwarffallagn/totaldwarfagn

print('\n\nDwarfs from RESOLVE Master Catalog')
print('Number of Dwarfs: {} \nSpring Dwarfs: {} ({:.2f}% of spring sample) \
      \nFall Dwarfs: {} ({:.2f}% of fall sample)'
      .format(totaldwarf,len(dwarfspring),pcdwarfspring,
              len(dwarffall),pcdwarffall))
#print('Total: {} \nSpring Dwarf AGN: {} ({:.2f}%) \nFall Dwarf AGN: {} ({:.2f}%)'
#      .format(totaldwarfagn,dwarfspringagn,pcdwarfspringagn,
#              dwarffallagn,pcdwarffallagn))
pcspringdwarfagn = 100.0*dwarfspringagn/len(dwarfspring)
pcfalldwarfagn = 100.0*dwarffallagn/len(dwarffall)
pcresdwarfagn = 100.0*np.sum(dwarfagn)/totaldwarf

springlowlim, springuplim = 100*binom_conf_interval(dwarfspringagn,len(dwarfspring))
springup= springuplim -pcspringdwarfagn
springlow = springlowlim-pcspringdwarfagn
falllowlim, falluplim = 100*binom_conf_interval(dwarffallagn,len(dwarffall))
fallup= falluplim -pcfalldwarfagn
falllow = falllowlim-pcfalldwarfagn
reslowlim, resuplim = 100*binom_conf_interval(np.sum(dwarfagn),totaldwarf)
resup= resuplim -pcresdwarfagn
reslow = reslowlim-pcresdwarfagn

def pcprint(pc,up,low):
    pc = str(round(pc,2))+'^{+'+str(round(up,2))+'}_{'+str(round(low,2))+'}\%'
    display(Math(pc))

display(Math('Dwarf AGN'))
display(Math('Spring : '+str(dwarfspringagn)+'/'+str(len(dwarfspring))))
pcprint(pcspringdwarfagn,springup,springlow)
display(Math('Fall : '+str(dwarffallagn)+'/'+str(len(dwarffall))))
pcprint(pcfalldwarfagn,fallup,falllow)
display(Math('RESOLVE : '+str(np.sum(dwarfagn))+'/'+str(totaldwarf)))
pcprint(pcresdwarfagn,resup,reslow)
#print('\nSpring Dwarf AGN: {}/{} = {:.2f}% of spring dwarf sample'
#      .format(dwarfspringagn,len(dwarfspring), 100.0*dwarfspringagn/len(dwarfspring)))
#print('Fall Dwarf AGN: {}/{} = {:.2f}% of fall dwarf sample'
#      .format(dwarffallagn,len(dwarffall), 100.0*dwarffallagn/len(dwarffall)))

###############################################################################
# Full ECO survey, including buffers
###############################################################################

eco = pd.read_csv('ECO_live22Oct2018.csv')
eco.index = eco.name

fallthreshold = 9.2
fullinobssample = 1
baryinobssample = 0
dec05 = np.array(eco.name[(eco.dedeg > 0) & (eco.dedeg < 5)])
other = np.array(eco.name[(eco.dedeg < 0) | (eco.dedeg > 5)])

totaleco = len(eco)
totalecodec05 = len(dec05)
totalecoother = len(other)
pcecodec05 = 100.0*totalecodec05/totaleco
pcecoother = 100.0*totalecoother/totaleco
print('Full ECO Survey (including buffers)')
print('Total: {} \nDec 0-5: {} ({:.2f}%) \nOther: {} ({:.2f}%)'
      .format(totaleco,totalecodec05,pcecodec05,totalecoother,pcecoother))

###############################################################################
# SDSS matches of full ECO survey (including buffers)
###############################################################################
sdsseco = pd.read_csv('ECO_full_raw.csv')
sdsseco.index = sdsseco.name

sdssdec05 = [x for x in sdsseco.name if x in dec05]
#np.array(sdsseco.name[(sdsseco.dedeg > 0) & (sdsseco.dedeg < 5)])
sdssother = [x for x in sdsseco.name if x in other]
#np.array(sdsseco.name[(sdsseco.dedeg < 0) | (sdsseco.dedeg > 5)])

totalsdsseco = len(sdsseco)
totalsdssecodec05 = len(sdssdec05)
totalsdssecoother = len(sdssother)
pcsdssecodec05 = 100.0*totalsdssecodec05/totalsdsseco
pcsdssecoother = 100.0*totalsdssecoother/totalsdsseco
print('\n\nFull ECO Survey (including buffers)')
print('Total: {} \nDec 0-5: {} ({:.2f}%) \nOther: {} ({:.2f}%)'
      .format(totalsdsseco,totalsdssecodec05,pcsdssecodec05,
              totalsdssecoother,pcsdssecoother))


###############################################################################
# ECO + inobssample
###############################################################################
ra=eco.radeg
dec=eco.dedeg
grpcz = eco.grpcz
cz = eco.cz
inother = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = eco.logmgas
mstars = eco.logmstar
mbary = 10**mgas + 10**mstars
ineco = (130.05 < eco.radeg) & (eco.radeg < 237.45)
ecoinobssample = (((eco.grpcz >= 3000.) & (eco.grpcz <= 7000.)) & 
                   (eco.absrmag < -17.33) & ineco)
dec05inobssample = ecoinobssample & (dec > 0) & (dec < 5)
otherinobssample = ecoinobssample & ((dec < 0) | (dec > 5))
totalecoinobs = np.sum(ecoinobssample)
totalecodec05inobs = np.sum(dec05inobssample)
totalecootherinobs = np.sum(otherinobssample)
pcecodec05inobs = 100.0*totalecodec05inobs/totalecoinobs
pcecootherinobs = 100.0*totalecootherinobs/totalecoinobs
#print(np.sum(((otherinobssample) & (mstars < 9.5))))
print('\n\ninobssample in full ECO survey')
print('Total: {} \nDec 0-5: {} ({:.2f}%) \nOther: {} ({:.2f}%)'
      .format(totalecoinobs, totalecodec05inobs, pcecodec05inobs,
              totalecootherinobs, pcecootherinobs))

###############################################################################
# ECO in SDSS + inobssample
###############################################################################
ra=eco.radeg.loc[sdsseco.name]
dec=eco.dedeg.loc[sdsseco.name]
grpcz = eco.grpcz.loc[sdsseco.name]
cz = eco.cz.loc[sdsseco.name]
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = eco.logmgas.loc[sdsseco.name]
mstars = eco.logmstar.loc[sdsseco.name]
mbary = 10**mgas + 10**mstars
absrmag = eco.absrmag.loc[sdsseco.name]
ineco = (130.05 < ra) & (ra < 237.45)
sdssecoinobssample = (((grpcz >= 3000.) & (grpcz <= 7000.)) & 
                   (absrmag < -17.33) & ineco)
sdssdec05inobssample = sdssecoinobssample & (dec > 0) & (dec < 5)
sdssotherinobssample = sdssecoinobssample & ((dec < 0) | (dec > 5))
totalsdssecoinobs = np.sum(sdssecoinobssample)
totalsdssecodec05inobs = np.sum(sdssdec05inobssample)
totalsdssecootherinobs = np.sum(sdssotherinobssample)
pcsdssecodec05inobs = 100.0*totalsdssecodec05inobs/totalsdssecoinobs
pcsdssecootherinobs = 100.0*totalsdssecootherinobs/totalsdssecoinobs
#print(np.sum(((otherinobssample) & (mstars < 9.5))))
print('\n\ninobssample of SDSS matches of ECO survey')
print('Total: {} \nDec 0-5: {} ({:.2f}%) \nOther: {} ({:.2f}%)'
      .format(totalsdssecoinobs, totalsdssecodec05inobs, pcsdssecodec05inobs,
              totalsdssecootherinobs, pcsdssecootherinobs))

###############################################################################
# ECO in SDSS & inobssample & SNR > 5 & 10**-3 < flux < 1e5 in ALL 3 CATALOGS
###############################################################################
master = pd.read_csv('ECO_snr5_master_new.csv')
master.index = master.name
masterdec05 = np.array(master.name[(master.dedeg > 0) & (master.dedeg < 5)])
masterother = np.array(master.name[(master.dedeg < 0) | (master.dedeg > 5)])
#masterspring = [x for x in list(master.name) if x[1] == 's']
#masterfall = [x for x in list(master.name) if x[1] == 'f']

totalmaster = len(master)
totalmasterdec05 = len(masterdec05)
totalmasterother = len(masterother)
pcmasterdec05 = 100.0*totalmasterdec05/totalmaster
pcmasterother = 100.0*totalmasterother/totalmaster
print('\n\nECO in SDSS & inobssample & SNR > 5 & 10**-3 < flux < 1e5 in ALL 3 CATALOGS')
print('Total: {} \ndec05: {} ({:.2f}%) \nother: {} ({:.2f}%)'
      .format(totalmaster,totalmasterdec05,pcmasterdec05,
              totalmasterother,pcmasterother))

###############################################################################
# AGN in ECO Master Catalog
###############################################################################
conf = pd.read_csv('ECO_snr5_master_conf_new.csv')
conf.index = conf.name

agn = conf.confidence_level > 0.0
#confspring = [x for x in list(conf.name) if x[1] == 's']
#conffall = [x for x in list(conf.name) if x[1] == 'f']
confdec05 =  np.array(master.name[(master.loc[conf.name].dedeg > 0) \
                                  & (master.loc[conf.name].dedeg < 5)])
#[x for x in list(conf.name) if x in list(sdsseco[dec05inobssample].index)]
confother = np.array(master.name[(master.loc[conf.name].dedeg < 0) \
                                 | (master.loc[conf.name].dedeg > 5)])
#[x for x in list(conf.name) if x in list(sdsseco[otherinobssample].index)]

totalconf = len(conf)
confdec05agn = len(conf.loc[confdec05][agn])
confotheragn = len(conf.loc[confother][agn])

pcconfdec05agn = 100.0*confdec05agn/totalconf
pcconfotheragn  = 100.0*confotheragn/totalconf
print('\n\nECO Master Catalog All AGN from Confidence Levels')
print('Total: {} \ndec05 AGN: {} ({:.2f}%) \nother AGN: {} ({:.2f}%)'
      .format(totalconf,confdec05agn,pcconfdec05agn,
              confotheragn,pcconfotheragn))
print('\ndec05 AGN: {}/{} = {:.2f}% of dec05 sample'
      .format(confdec05agn,len(confdec05), 100.0*confdec05agn/len(confdec05)))
print('\nother AGN: {}/{} = {:.2f}% of other sample'
      .format(confotheragn,len(confother), 100.0*confotheragn/len(confother)))

###############################################################################
# Dwarf AGN in ECO Master Catalog
###############################################################################
dwarf = master[master.logmstar < 9.5]

dwarfdec05_condition = (master.loc[dwarf.name].dedeg > 0) \
                                  & (master.loc[dwarf.name].dedeg < 5)
dwarfother_condition = (master.loc[dwarf.name].dedeg < 0) \
                                 | (master.loc[dwarf.name].dedeg > 5)
dwarfdec05 = np.array(dwarfdec05_condition.index[dwarfdec05_condition])
dwarfother = np.array(dwarfother_condition.index[dwarfother_condition])
totaldwarf = len(dwarf)

pcdwarfdec05 = 100.0*len(dwarfdec05)/totalmasterdec05
pcdwarfother = 100.0*len(dwarfother)/totalmasterother

dwarfagn = conf.loc[list(dwarf.name)].confidence_level > 0.0

dwarfdec05agn = len(dwarf.loc[dwarfdec05][dwarfagn])
dwarfotheragn = len(dwarf.loc[dwarfother][dwarfagn])
totaldwarfagn = dwarfdec05agn+dwarfotheragn

pcdwarfdec05agn = 100.0*dwarfdec05agn/totaldwarfagn
pcdwarfotheragn  = 100.0*dwarfotheragn/totaldwarfagn

print('\n\nDwarfs from ECO Master Catalog')
print('Number of Dwarfs: {} \ndec05 Dwarfs: {} ({:.2f}% of dec05 sample) \
      \nother Dwarfs: {} ({:.2f}% of other sample)'
      .format(totaldwarf,len(dwarfdec05),pcdwarfdec05,
              len(dwarfother),pcdwarfother))
#print('Total: {} \ndec05 Dwarf AGN: {} ({:.2f}%) \nother Dwarf AGN: {} ({:.2f}%)'
#      .format(totaldwarfagn,dwarfdec05agn,pcdwarfdec05agn,
#              dwarfotheragn,pcdwarfotheragn))
pcdec05dwarfagn = 100.0*dwarfdec05agn/len(dwarfdec05)
pcotherdwarfagn = 100.0*dwarfotheragn/len(dwarfother)
pcecodwarfagn = 100.0*np.sum(dwarfagn)/totaldwarf

dec05lowlim, dec05uplim = 100*binom_conf_interval(dwarfdec05agn,len(dwarfdec05))
dec05up= dec05uplim -pcdec05dwarfagn
dec05low = dec05lowlim-pcdec05dwarfagn
otherlowlim, otheruplim = 100*binom_conf_interval(dwarfotheragn,len(dwarfother))
otherup= otheruplim -pcotherdwarfagn
otherlow = otherlowlim-pcotherdwarfagn
ecolowlim, ecouplim = 100*binom_conf_interval(np.sum(dwarfagn),totaldwarf)
ecoup= ecouplim -pcecodwarfagn
ecolow = ecolowlim-pcecodwarfagn

def pcprint(pc,up,low):
    pc = str(round(pc,2))+'^{+'+str(round(up,2))+'}_{'+str(round(low,2))+'}\%'
    display(Math(pc))

display(Math('Dwarf AGN'))
display(Math('dec05 : '+str(dwarfdec05agn)+'/'+str(len(dwarfdec05))))
pcprint(pcdec05dwarfagn,dec05up,dec05low)
display(Math('other : '+str(dwarfotheragn)+'/'+str(len(dwarfother))))
pcprint(pcotherdwarfagn,otherup,otherlow)
display(Math('ECO : '+str(np.sum(dwarfagn))+'/'+str(totaldwarf)))
pcprint(pcecodwarfagn,ecoup,ecolow)

#print('\ndec05 Dwarf AGN: {}/{} = {:.2f}% of dec05 dwarf sample'
#      .format(dwarfdec05agn,len(dwarfdec05), 100.0*dwarfdec05agn/len(dwarfdec05)))
#print('other Dwarf AGN: {}/{} = {:.2f}% of other dwarf sample'
#      .format(dwarfotheragn,len(dwarfother), 100.0*dwarfotheragn/len(dwarfother)))

###############################################################################
#Compare the distributions of ECO galaxies in Dec 0-5 and Dec > 5

import matplotlib.pyplot as plt
plt.figure('Stellar Mass')
n, bins, patches = plt.hist(dwarf.loc[dwarfdec05].logmstar, bins = 'fd', 
                            density = True, color = 'gray')
plt.hist(dwarf.loc[dwarfother].logmstar, bins = bins, histtype = 'step',
         density = True, color = 'orange', linewidth = 3)
med_dec05 = np.median(dwarf.loc[dwarfdec05].logmstar)
med_other = np.median(dwarf.loc[dwarfother].logmstar)
yaxis = np.arange(0,2.5)
plt.plot(med_dec05*np.ones(len(yaxis)),yaxis, color = 'blue')
plt.plot(med_other*np.ones(len(yaxis)),yaxis, color = 'orange')

plt.figure('Gas Mass')
n, bins, patches = plt.hist(dwarf.loc[dwarfdec05].logmgas, bins = 'fd', 
                            density = True, color = 'gray')
plt.hist(dwarf.loc[dwarfother].logmgas, bins = bins, histtype = 'step',
         density = True, color = 'orange', linewidth = 3)
med_dec05 = np.median(dwarf.loc[dwarfdec05].logmgas)
med_other = np.median(dwarf.loc[dwarfother].logmgas)
yaxis = np.arange(0,2.5)
plt.plot(med_dec05*np.ones(len(yaxis)),yaxis, color = 'blue')
plt.plot(med_other*np.ones(len(yaxis)),yaxis, color = 'orange')

plt.figure('Group Halo Mass')
n, bins, patches = plt.hist(dwarf.loc[dwarfdec05].logmh, bins = 'fd', 
                            density = True, color = 'gray')
plt.hist(dwarf.loc[dwarfother].logmh, bins = bins, histtype = 'step',
         density = True, color = 'orange', linewidth = 3)
med_dec05 = np.median(dwarf.loc[dwarfdec05].logmh)
med_other = np.median(dwarf.loc[dwarfother].logmh)
yaxis = np.arange(0,2.5)
plt.plot(med_dec05*np.ones(len(yaxis)),yaxis, color = 'blue')
plt.plot(med_other*np.ones(len(yaxis)),yaxis, color = 'orange')

plt.figure('r-band Luminosity')
n, bins, patches = plt.hist(dwarf.loc[dwarfdec05].absrmag, bins = 'fd', 
                            density = True, color = 'gray')
plt.hist(dwarf.loc[dwarfother].absrmag, bins = bins, histtype = 'step',
         density = True, color = 'orange', linewidth = 3)
med_dec05 = np.median(dwarf.loc[dwarfdec05].absrmag)
med_other = np.median(dwarf.loc[dwarfother].absrmag)
yaxis = np.arange(0,2.5)
plt.plot(med_dec05*np.ones(len(yaxis)),yaxis, color = 'blue')
plt.plot(med_other*np.ones(len(yaxis)),yaxis, color = 'orange')

plt.figure('Baryonic Mass')
dec05_bary = np.log10(10**dwarf.loc[dwarfdec05].logmgas + \
                      10**dwarf.loc[dwarfdec05].logmstar)
other_bary = np.log10(10**dwarf.loc[dwarfother].logmgas + \
                      10**dwarf.loc[dwarfother].logmstar)
n, bins, patches = plt.hist(dec05_bary, bins = 'fd', 
                            density = True, color = 'gray')
plt.hist(other_bary, bins = bins, histtype = 'step',
         density = True, color = 'orange', linewidth = 3)
med_dec05 = np.median(dec05_bary)
med_other = np.median(other_bary)
yaxis = np.arange(0,2.5)
plt.plot(med_dec05*np.ones(len(yaxis)),yaxis, color = 'blue')
plt.plot(med_other*np.ones(len(yaxis)),yaxis, color = 'orange')


#resdata = readsav('../SDSS_spectra/resolvecatalog.dat')
#resphot = readsav('../SDSS_spectra/resolvecatalogphot.dat')
#
#ndx1 = [x for x in range(len(resdata)) if resdata.name[x] in dwarfdec05]
#ndx2 = [x for x in range(len(resdata)) if resdata.name[x] in dwarfother]
#

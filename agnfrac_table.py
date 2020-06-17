# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:11:09 2020

@author: mugdhapolimera

Tables for RESOLVE and ECO Dwarf AGN %
"""
import pandas as pd
import numpy as np
from astropy.stats import binom_conf_interval

resfull = pd.read_csv('RESOLVE_snr5_master_bary.csv')
resflag = pd.read_csv('resolve_emlineclass_full_bary_master_new.csv')
resconf = pd.read_csv('RESOLVE_snr5_master_conf_bary.csv')
resfull['logmbary'] = np.log10(10**resfull.logmstar + 10**resfull.logmgas)
resfull.index = resfull.name
resflag.index = resflag.galname
resconf.index = resconf.name
resconf.confidence_level = resconf.jhu + resconf.port
resndx = (resfull.logmbary > 9.2) # | (resfull.absrmag < -17.33) 
res = resfull[resndx]
resflag = resflag[resndx]
resconf = resconf.loc[resndx]

ecofull = pd.read_csv('ECO_snr5_master_bary.csv')
ecofullflag = pd.read_csv('eco_emlineclass_full_bary_master_new.csv')
ecoconf = pd.read_csv('ECO_snr5_master_conf_bary.csv')
ecofull.index = ecofull.name
ecofullflag.index = ecofullflag.galname
ecofull['logmbary'] = np.log10(10**ecofull.logmstar + 10**ecofull.logmgas)
econdx = (((ecofull.cz < 4500) | (ecofull.dedeg < 0) | (ecofull.dedeg > 5)) \
          & ((ecofull.logmbary > 9.2))) #| \
#         (ecofull.resname == 'notinresolve') )
#eco = ecofull[econdx]
#ecoflag = ecofullflag[econdx]
ecoconf.index = ecoconf.name
ecoconf.confidence_level = ecoconf.jhu + ecoconf.port
eco = ecofull[econdx]# | (ecofull.absrmag < -17.33))] 
ecoflag = ecofullflag[econdx]# | (ecofull.absrmag < -17.33))]
ecoconf = ecoconf.loc[econdx]

table = pd.DataFrame(columns = ['resdwarf','resdwarfagn','resdwarfagnpc',\
                                   'ecodwarf','ecodwarfagn','ecodwarfagnpc'])

resdwarf = res[res.logmstar < 9.5]
ecodwarf = eco[eco.logmstar < 9.5]
res_agn = resconf.confidence_level >= 0
eco_agn = ecoconf.confidence_level >= 0
ecodwarfagn = ecodwarf[eco_agn]
resdwarfagn = resdwarf[res_agn]
resdwarfconf = resconf.loc[resdwarf.name]
ecodwarfconf = ecoconf.loc[ecodwarf.name]
indexes = ['JHU', 'Portsmouth', 'JHU OR Portsmouth', 'JHU AND Portsmouth']
colname = {'JHU': 'jhu', 'Portsmouth':'port',
           'JHU OR Portsmouth': 'jhu+port',
           'JHU AND Portsmouth': 'jhu&port'}
#print('\\begin{tabular}{|c | c | c | c| c| c| c| } \n \
#        \\toprule \n \
#        {} &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
#        &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \
#        \hline')
#print('\\begin{table*} \ \n\centering \
#    \n \t \\begin{tabularx}{\\textwidth} { \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X \
#  | >{\centering\\arraybackslash}X | }\
#    \n \t \\toprule \
#    \n \t Category &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
#    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \hline')
print('\\begin{table*} \ \n \centering \
      \n \t \\begin{tabular}{| C{2cm} | C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |} \
      \n \t \\toprule \ \n \t {} & \multicolumn{3}{|c|}{RESOLVE} & \multicolumn{3}{|c|}{ECO} \\\\ \hline \
         \n \t Category &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \n \t \hline')
         
for i in range(len(indexes)):
    index = indexes[i] 
    
    if ('+' in colname[index]): 
        rdwarf = np.sum((resdwarfconf.jhu != 0) | (resdwarfconf.port != 0))
        rdwarfagn = np.sum(resdwarfconf['confidence_level'] >= 0)        
        edwarf = np.sum((ecodwarfconf.jhu != 0) | (ecodwarfconf.port != 0))
        edwarfagn = np.sum(ecodwarfconf['confidence_level'] >= 0)        

    elif('&' in colname[index]):
        rdwarf = np.sum((resdwarfconf.jhu != 0) & (resdwarfconf.port != 0))
        rdwarfagn = np.sum(resdwarfconf['confidence_level'] == 2)        
        edwarf = np.sum((ecodwarfconf.jhu != 0) & (ecodwarfconf.port != 0))
        edwarfagn = np.sum(ecodwarfconf['confidence_level'] == 2)        

    else:
        rdwarf = np.sum(resdwarfconf[colname[index]] != 0)
        rdwarfagn = np.sum(resdwarfconf[colname[index]] > 0)        
        edwarf = np.sum(ecodwarfconf[colname[index]] != 0)
        edwarfagn = np.sum(ecodwarfconf[colname[index]] > 0)        
    rdwarfagnpc = round((100.0*rdwarfagn/rdwarf),2)
    r_edown, r_eup = 100.0*binom_conf_interval(rdwarfagn, rdwarf) - rdwarfagnpc
    r_edown = round(-r_edown,2)
    r_eup = round(r_eup,2)

    edwarfagnpc = round((100.0*edwarfagn/edwarf),2)
    e_edown, e_eup = 100.0*binom_conf_interval(edwarfagn, edwarf) - edwarfagnpc
    e_edown = round(-e_edown,2)
    e_eup = round(e_eup,2)

    print '\t'+index+' & '+ str(rdwarf)+' & ',str(rdwarfagn)+' & $'+str(rdwarfagnpc)+\
          '^{+'+str(r_eup)+'}'+'_{'+str(-r_edown)+'}$'+' & '\
          + str(edwarf)+' & ',str(edwarfagn)+' & $'+str(edwarfagnpc)+\
          '^{+'+str(e_eup)+'}'+'_{'+str(-e_edown)+'}$\\\\'
print('\t \hline \n \t \end{tabular} \n \label{table:2} \n \end{table*}')

print('\\begin{table*} \ \n \centering \
      \n \t \\begin{tabular}{| C{2cm} | C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |} \
      \n \t \\toprule \ \n \t {} & \multicolumn{3}{|c|}{RESOLVE} & \multicolumn{3}{|c|}{ECO} \\\\ \hline \
         \n \t Category &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \n \t \hline')

indexes = ['Confidence $>=$ 0', 'Confidence $>=$ 1', 'Confidence $=$ 2']
for i in [0,1,2]:
    index = indexes[i]
    if i == 2:
        rdwarf = np.sum((resdwarfconf.jhu != 0) & (resdwarfconf.port != 0))
        rdwarfagn = np.sum(resdwarfconf['confidence_level'] == 2)        
        edwarf = np.sum((ecodwarfconf.jhu != 0) & (ecodwarfconf.port != 0))
        edwarfagn = np.sum(ecodwarfconf['confidence_level'] == 2)        
    else:        
        rdwarf = len(resdwarfconf)
        rdwarfagn = np.sum(resdwarfconf['confidence_level'] >= i)        
        edwarf = len(ecodwarfconf)
        edwarfagn = np.sum(ecodwarfconf['confidence_level'] >= i)        

    rdwarfagnpc = round((100.0*rdwarfagn/rdwarf),2)
    r_edown, r_eup = 100.0*binom_conf_interval(rdwarfagn, rdwarf) - rdwarfagnpc
    r_edown = round(-r_edown,2)
    r_eup = round(r_eup,2)

    edwarfagnpc = round((100.0*edwarfagn/edwarf),2)
    e_edown, e_eup = 100.0*binom_conf_interval(edwarfagn, edwarf) - edwarfagnpc
    e_edown = round(-e_edown,2)
    e_eup = round(e_eup,2)

    print '\t'+index+' & '+ str(rdwarf)+' & ',str(rdwarfagn)+' & $'+str(rdwarfagnpc)+\
          '^{+'+str(r_eup)+'}'+'_{'+str(-r_edown)+'}$'+' & '\
          + str(edwarf)+' & ',str(edwarfagn)+' & $'+str(edwarfagnpc)+\
          '^{+'+str(e_eup)+'}'+'_{'+str(-e_edown)+'}$\\\\'
print('\t \hline \n \t \end{tabular} \n \label{table:2} \n \end{table*}')


#resjhuandport = (resconf.confidence_level == 2) | \
#            (resconf.confidence_level == -2) |(resconf.confidence_level == 0)
#resjhuandportname_agn = list(resconf.loc[resdwarfagn.name][resconf.confidence_level == 2].name)
#resjhuandportname = list(resconf.loc[resdwarf.name][resjhuandport].name)
#
#ecojhuandport = (ecoconf.confidence_level == 2) | \
#            (ecoconf.confidence_level == -2) |(ecoconf.confidence_level == 0)
#ecojhuandportname_agn = list(ecoconf.loc[ecodwarfagn.name][ecoconf.confidence_level == 2].name)
#ecojhuandportname = list(ecoconf.loc[ecodwarf.name][ecojhuandport].name)
#
#

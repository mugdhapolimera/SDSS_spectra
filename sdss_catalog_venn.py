# -*- coding: utf-8 -*-
"""
Created on Mon May 17 20:56:08 2021

@author: mugdhapolimera
"""
from matplotlib_venn import venn2, venn3
import pandas as pd
import matplotlib.pyplot as plt
#resolve = 0
#eco = 1
##Catalog crossmatch with full RESOLVE/ECO
#if resolve:
#    jhu = pd.read_csv("RESOLVE_raw_jhu.csv")
#    port = pd.read_csv("RESOLVE_raw_port.csv")
#    nsa = pd.read_csv("RESOLVE_raw_nsa.csv")
#if eco:
#    jhu = pd.read_csv("ECO_raw_jhu.csv")
#    port = pd.read_csv("ECO_raw_port.csv")
#    nsa = pd.read_csv("ECO_raw_nsa.csv")
#set1 = set(jhu.name)
#set2 = set(nsa.name)
#set3 = set(port.name)
#plt.figure()
#venn3([set1, set2, set3], ('JHU', 'NSA', 'Port'), set_colors = ('r', 'b', 'y'))
#
##EL samples
#if resolve:
#    jhu = pd.read_csv("RESOLVE_full_hasnr5_dext_jhu.csv")
#    port = pd.read_csv("RESOLVE_full_hasnr5_dext_port.csv")
#    nsa = pd.read_csv("RESOLVE_full_hasnr5_dext_nsa.csv")
#if eco:
#    jhu = pd.read_csv("ECO_full_hasnr5_dext_jhu.csv")
#    port = pd.read_csv("ECO_full_hasnr5_dext_port.csv")
#    nsa = pd.read_csv("ECO_full_hasnr5_dext_nsa.csv")
#set1 = set(jhu.name)
#set2 = set(nsa.name)
#set3 = set(port.name)
#print(len(set1), len(set3), len(set2))
#print(len(set1 & set3), len(set2 & set3), len(set1 & set2))
#plt.figure()
#venn3([set1, set2, set3], ('JHU', 'NSA', 'Port'), set_colors = ('r', 'b', 'y'))
#
##SEL samples
#if resolve:
#    jhu = pd.read_csv("RESOLVE_full_snr5_dext_jhu.csv")
#    port = pd.read_csv("RESOLVE_full_snr5_dext_port.csv")
#    nsa = pd.read_csv("RESOLVE_full_snr5_dext_nsa.csv")
#if eco:
#    jhu = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_jhu.csv")
#    port = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_port.csv")
#    nsa = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_nsa.csv")
#set1 = set(jhu.name)
#set2 = set(nsa.name)
#set3 = set(port.name)
#print(len(set1), len(set3), len(set2))
#print(len(set1 & set3), len(set2 & set3), len(set1 & set2))
#plt.figure()
#venn3([set1, set2, set3], ('JHU', 'NSA', 'Port'), set_colors = ('r', 'b', 'y'))
#


# set1 = set(jhu[jhuflagagn & (jhu.logmstar < 9.5)].name)
# set2 = set(nsa[nsaflagagn & (nsa.logmstar < 9.5)].name)
# set3 = set(port[portflagagn & (port.logmstar < 9.5)].name)

#S06 and SEL samples - RESOLVE + ECO
jhus06 = pd.read_csv("ECO+RESOLVE_bpt1snr5_dext_jhu.csv")
ports06 = pd.read_csv("ECO+RESOLVE_bpt1snr5_dext_port.csv")
nsas06 = pd.read_csv("ECO+RESOLVE_bpt1snr5_dext_nsa.csv")

jhusel = pd.read_csv("ECO+RESOLVE_snr5_dext_jhu.csv")
portsel = pd.read_csv("ECO+RESOLVE_snr5_dext_port.csv")
nsasel = pd.read_csv("ECO+RESOLVE_snr5_dext_nsa.csv")

#set1 = set(jhus06.name)
#set2 = set(jhusel.name)
#plt.figure()
#venn2([set1, set2], ('JHU S06', 'JHU SEL'), set_colors = ('r', 'b'))
#
#set1 = set(ports06.name)
#set2 = set(portsel.name)
#plt.figure()
#venn2([set1, set2], ('Port S06', 'Port SEL'), set_colors = ('r', 'b'))
#
#set1 = set(nsas06.name)
#set2 = set(nsasel.name)
#plt.figure()
#venn2([set1, set2], ('NSA S06', 'NSA SEL'), set_colors = ('r', 'b'))

#IR, S06, SEL  samples - RESOLVE
res_ir = pd.read_csv('mid_ir/RESOLVE_WISE_good.csv')
jhus06 = pd.read_csv("RESOLVE_full_bpt1snr5_dext_jhu.csv")
ports06 = pd.read_csv("RESOLVE_full_bpt1snr5_dext_port.csv")
nsas06 = pd.read_csv("RESOLVE_full_bpt1snr5_dext_nsa.csv")

jhusel = pd.read_csv("RESOLVE_full_snr5_dext_jhu.csv")
portsel = pd.read_csv("RESOLVE_full_snr5_dext_port.csv")
nsasel = pd.read_csv("RESOLVE_full_snr5_dext_nsa.csv")

fig,(ax1,ax2,ax3) = plt.subplots(1,3)
set1 = set(jhus06.name)
set2 = set(res_ir.name)
set3 = set(jhusel.name)
ax1 = venn3([set1, set2, set3], ('JHU S06', 'Mid IR', 'JHU SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax1)

set1 = set(ports06.name)
set2 = set(res_ir.name)
set3 = set(portsel.name)
ax2 = venn3([set1, set2, set3], ('Port S06', 'Mid IR', 'Port SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax2)

set1 = set(nsas06.name)
set2 = set(res_ir.name)
set3 = set(nsasel.name)
ax3 = venn3([set1, set2, set3], ('NSA S06', 'Mid IR', 'NSA SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax3)

#ECO
res_ir = pd.read_csv('mid_ir/ECO_WISE_good.csv')
jhus06 = pd.read_csv("ECO_full_bpt1snr5_dext_jhu.csv")
ports06 = pd.read_csv("ECO_full_bpt1snr5_dext_port.csv")
nsas06 = pd.read_csv("ECO_full_bpt1snr5_dext_nsa.csv")

jhusel = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_jhu.csv")
portsel = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_port.csv")
nsasel = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_nsa.csv")

fig,(ax1,ax2,ax3) = plt.subplots(1,3)
set1 = set(jhus06.name)
set2 = set(res_ir.name)
set3 = set(jhusel.name)
ax1 = venn3([set1, set2, set3], ('JHU S06', 'Mid IR', 'JHU SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax1)

set1 = set(ports06.name)
set2 = set(res_ir.name)
set3 = set(portsel.name)
ax2 = venn3([set1, set2, set3], ('Port S06', 'Mid IR', 'Port SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax2)

set1 = set(nsas06.name)
set2 = set(res_ir.name)
set3 = set(nsasel.name)
ax3 = venn3([set1, set2, set3], ('NSA S06', 'Mid IR', 'NSA SEL'), 
      set_colors = ('r', 'b', 'y'), ax = ax3)


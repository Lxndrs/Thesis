
"""Plot data from table of pressures and fluxes from proplyd arc fits
   Original: William Henney.
   Modified Version: Jorge A. Tarango-Yong
"""

from astropy.table import Table
from matplotlib import pyplot as plt
import json
import numpy as np
import seaborn as sns
import argparse

def split_list(tab):
    """
    Some elements in the input table include uncertainties in the form 'x +/- dx'. 
    So is neccesary to separate them into two arrays.
    Input: Element from astropy table.
    Output:
        x: Measurements
        dx: Uncertainties
    """
    x = []
    dx = []
    for t in tab:
        x.append(float(t.split("+/-")[0]))
        try:
            dx.append(float(t.split("+/-")[1]))
        except IndexError:
            dx.append(0.0)
    return x, dx

parser = argparse.ArgumentParser(
    description=""" Choose a Winds Parameters table""")

parser.add_argument("--table", type=str,
                    default="wind-fits.tab",
                    help=" Choose a Winds Parameters table ")

cmd_args = parser.parse_args()
table = cmd_args.table

# Useful Physical and Astronomical constants.
AU = 1.49597870691e13
PC = 3.085677582e18
d_Orion = 414.0 # Use Menten et al. (2007) Measurement
k_Boltzmann = 1.3806503e-16
cos80 = 0.173648177667

# Determine which model we are using
# Model A: Uses inner wind density from HA1998 but inclinations obtained in this work
# Model B: Derive inner wind density from our determinations of $\beta$, depends of outer wind pressure.
# Model B.1: uses outer wind parameters from GAH:2002
# Model B.2: uses cool model from Gagne:2005
# Model HA: Uses inner wind density and inclinations from HA1998
# the input table assumes certain model.
# In sake of consistency, is neccesary to distinguish model B from the others

model = "A" # assumes model A by default

if "beta" in table:
    model = "B.1"
    if "Gagne" in table:
        model = "B.2"
if "HA98" in table:
    model = "HA"

# Read parameters table
tab = Table.read('../'+table, format='ascii.tab')
sources = sorted(set(tab['Fuente']))
n = len(sources)
#colors = sns.color_palette('Set1', n)
sns.set_style("whitegrid")
# Physical separations in parsec
D, Delta_D = split_list(tab["D(as)"])
D_pc = np.array(D)*d_Orion*AU/PC
Delta_D_pc = D_pc*(np.array(Delta_D)/np.array(D))

#collection of hex colors
dark_blue = "#1e25b6"
pearl_turquoise ="#32c6a6"
mexican_pink = "#e4007c"
crimson = "#dc143c"
leaf_green = "#15ae26"
brown = "#b6451e"
gray = "#515952"
guinda = "#aa1c47"
gold = "#FFD700"
orange = "#E08000"
#Create a dictionary with hex colors for the objects
colordict = {"LV2":dark_blue, "LV2b":pearl_turquoise, "LV3":mexican_pink, "LV4":crimson, "LV5":brown, "168-328":leaf_green, "169-338":gray, "177-341":guinda, "180-331":orange}


#split F(star), P(wind) and F(ph)/F(*)+  columns
#Fs, dFs = split_list(tab["F(star)"])
#Fs, dFs = np.array(Fs), np.array(dFs)
#Pw,
# Measuring stellar wind RAM pressure
# Units in cgs system
Mdot = 2.206e19
Vw = 1.2e8

# Mass loss rate and terminal velocity of \thC{} taken from Gagne et al. 2005, table 3.
# There are two models: cool model (for an O7 star) and hot model (for a O5.5 star) (hot model removed)
Mdot_cool = 3.47e19
#Mdot_hot = 8.83e19
Vw_cool = 2.76e8
#Vw_hot = 2.98e8
D_arr = np.logspace(-6, 0)
Pw = Mdot*Vw/(4*np.pi*(D_arr*PC)**2)
Pw_cool = Mdot_cool*Vw_cool/(4*np.pi*(D_arr*PC)**2)
#Pw_hot = Mdot_hot*Vw_hot/(4*np.pi*(D_arr*PC)**2)

# Measuring stellar flux
Qh = 1e49 # Number of ionizing photons per second
fd = 0.5 # Asortion by dust factor
Fs = Qh*(1-fd)/(4*np.pi*(D_arr*PC)**2)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.loglog(D_arr, Fs, c='k', alpha=0.1, lw=10, label='')
if model == "A": # in model A can display both outer winds: GAH:2002 and cool model
    ax2.loglog(D_arr, Pw/k_Boltzmann, c='k', alpha=0.1, lw=10, label='')
    ax2.loglog(D_arr, Pw_cool/k_Boltzmann, c='b', alpha=0.1, lw=10, label='')
elif (model == "B.1") or (model=="HA"): # in model B.1 or HA only the GAH:2002 must be displayed
    ax2.loglog(D_arr, Pw/k_Boltzmann, c='k', alpha=0.1, lw=10, label='')
else: # in the B.2 model, only the cool model must be displayed
    ax2.loglog(D_arr, Pw_cool/k_Boltzmann, c='b', alpha=0.1, lw=10, label='')  
#ax2.loglog(D_arr, Pw_hot/k_Boltzmann, c='r', alpha=0.1, lw=10, label='')
mm = tab["*"] == "**" # Subsamples where photoionization balance is accomplished
nn = tab["*"] == "*"  # Subsamples where photoionization balance is at least weakly accomplished.
used = tab["*"] == "x" # Subsamples where photoionization balance is at least poorly accomplished
#out_colnames = ['Source' ,'D prime', 'R0/D prime', 'Rc/R0 prime full',
#                'Rc/R0 prime select', 'beta', 'xi', 'inc', 'D', 'R0/D']

#out_formats = {
#    'D prime': '%.1f',
#    'R0/D prime':         '%.4f',
#    'Rc/R0 prime full':   '%.4f',
#    'Rc/R0 prime select': '%.4f',
#    'beta':               '%.4f',
#    'xi':                 '%.4f',
#    'inc':                '%.4f',
#    'D':                  '%.4f',
#    'R0/D':               '%.4f',
#
#
#out_rows = []

#def var_range(x, dx):
#    return x, dx


for source in sources:
    m = tab['Fuente'] == source
    Dprime = tab["D'(as)"][m][0] * d_Orion*AU/PC
    try:
        Fp, dFp = split_list(tab['F(photo)'])
        Fp, dFp = np.array(Fp), np.array(dFp)
        Pi, dPi = split_list(tab["P(in)"])
        Pi, dPi = np.array(Pi), np.array(dPi)
    except KeyError:
        Fp, dFp_plus = split_list(tab['F(photo)+'])
        Fp, dFp_minus = split_list(tab['F(photo)-'])
        Fp = np.array(Fp)
        dFp = np.array([dFp_minus, dFp_plus])
        Pi, dPi_plus = split_list(tab["P(in)+"])
        Pi, dPi_minus = split_list(tab["P(in)-"])
        Pi = np.array(Pi)
        dPi = np.array([dPi_minus, dPi_plus])
    # Get data from variational fits
#    combined_file = '../saves/LV-bowshocks-xyfancy-variations-{}.save'.format(source)
#    vardata = json.load(open(combined_file))
    # Combine with data from org table to fill in a new output table (I don't need this anymore)
#    Rcp_select = tab['\Pi\''][m & mm]
#    Rcp_full = np.array(vardata['Rc'])/np.array(vardata['R0'])
#    beta = tab[r'\beta'][m & mm]
#    xi = tab['xi'][m & mm]
#    inc = tab['i (deg)'][m & mm]
#    DD = tab['D(as)'][m & mm] * d_Orion*AU/PC
#    R0_D = tab['q'][m & mm]

# Actually we don't need the table work since i already did it
#    out_rows.append({
#        'Source': source,
#        'D prime': tab['D\'(as)'][m][0],
#        'R0/D prime':         var_range(np.mean(vardata['R0']), np.std(vardata['R0'])),
#        'Rc/R0 prime full':    var_range(np.mean(Rcp_full), np.std(Rcp_full)), 
#        'Rc/R0 prime select':  var_range(np.mean(Rcp_select), np.std(Rcp_select)),
#        'beta':               var_range(np.mean(beta), np.std(beta)), 
#        'xi':                 var_range(np.min(xi), np.max(xi)),
#        'inc':                var_range(np.mean(inc), np.std(inc)),
#        'D':                   var_range(np.mean(DD), np.std(DD)),
#        'R0/D':               var_range(np.mean(R0_D), np.std(R0_D)),
#    })
                   
    ax1.plot([Dprime, Dprime/cos80], [Fp[m], Fp[m]], ':', c=colordict[source], alpha=0.4, label='')
    ax1.plot(D_pc[m & used], Fp[m & used], linestyle="None", marker=".", c=colordict[source], alpha=0.4, label='')
    ax1.plot(D_pc[m & mm], Fp[m & mm],
               linestyle="None", marker="o", lw=3, c=colordict[source], label=source)
    ax2.plot(D_pc[m & used], Pi[m & used]/k_Boltzmann,
               linestyle="None", marker=".", c=colordict[source], alpha=0.4, label='')
    ax2.plot(D_pc[m & mm], Pi[m & mm]/k_Boltzmann,
               linestyle="None", marker="o", c=colordict[source], lw=3, label=source)
    ax1.plot(D_pc[m & nn], Fp[m & nn],
               linestyle="None", marker="o", lw=0.5, c=colordict[source], label='', alpha=0.4)
    ax2.plot(D_pc[m & nn], Pi[m & nn]/k_Boltzmann,
               linestyle="None", marker="o", c=colordict[source], lw=0.5, label='', alpha=0.4)
    try:
        ax1.errorbar(D_pc[m & used], Fp[m & used], xerr=Delta_D_pc[m & used],yerr=dFp[m & used], c=colordict[source], alpha=0.4)
        ax1.errorbar(D_pc[m & mm], Fp[m & mm], xerr=Delta_D_pc[m & mm], yerr=dFp[m & mm], c=colordict[source])
        ax2.errorbar(D_pc[m & used], Pi[m & used]/k_Boltzmann, xerr=Delta_D_pc[m & used], yerr=dPi[m & used]/k_Boltzmann, c=colordict[source], alpha=0.4)
        ax2.errorbar(D_pc[m & mm], Pi[m & mm]/k_Boltzmann, xerr=Delta_D_pc[m & mm], yerr=dPi[m & mm]/k_Boltzmann, c=colordict[source])
        ax1.errorbar(D_pc[m & nn], Fp[m & nn], xerr = Delta_D_pc[m & (nn & ~mm)], yerr=dFp[m & nn], c=colordict[source], alpha=0.4)
        ax2.errorbar(D_pc[m & nn], Pi[m & nn]/k_Boltzmann, xerr=Delta_D_pc[m & nn], yerr=dPi[m & nn]/k_Boltzmann, c=colordict[source], alpha=0.4)
    except:
        ax1.errorbar(D_pc[m & used], Fp[m & used], xerr=Delta_D_pc[m & used],yerr=dFp[:, m & used], c=colordict[source], alpha=0.4)
        ax1.errorbar(D_pc[m & mm], Fp[m & mm], xerr=Delta_D_pc[m & mm], yerr=dFp[:, m & mm], c=colordict[source])
        ax2.errorbar(D_pc[m & used], Pi[m & used]/k_Boltzmann, xerr=Delta_D_pc[m & used], yerr=dPi[:, m & used]/k_Boltzmann, c=colordict[source], alpha=0.4)
        ax2.errorbar(D_pc[m & mm], Pi[m & mm]/k_Boltzmann, xerr=Delta_D_pc[m & mm], yerr=dPi[:, m & mm]/k_Boltzmann, c=colordict[source])
        ax1.errorbar(D_pc[m & nn], Fp[m & nn], xerr = Delta_D_pc[m & (nn & ~mm)], yerr=dFp[:, m & nn], c=colordict[source], alpha=0.4)
        ax2.errorbar(D_pc[m & nn], Pi[m & nn]/k_Boltzmann, xerr=Delta_D_pc[m & nn], yerr=dPi[:, m & nn]/k_Boltzmann, c=colordict[source], alpha=0.4)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.legend(ncol=2, loc='lower left')
ax2.set_xlim(0.008, 0.3)
ax2.set_ylim(1e7, 4e9)
ax1.set_ylim(2e12, 3e14)
ax1.set_ylabel(r'Ionizing Flux, $\mathrm{cm^{-2}\ s^{-1}}$')
ax2.set_ylabel(r'Stagnation Pressure: $P/k$, $\mathrm{cm^{-3}\ K}$')
ax2.set_xlabel('Distance, parsec')
fig.set_size_inches(5, 8)
fig.tight_layout()
output = table.replace(".tab", ".pdf")
fig.savefig('../Figures/plot-' + output)

#out_tab = Table(names=out_colnames, rows=out_rows)

#print(out_tab.pprint(max_width=-1, max_lines=-1))

#out_tab.write('../arc-fit-table-for-paper.tab',
#              format='ascii.fixed_width', formats=out_formats)


"""Original Version:
from astropy.table import Table
from matplotlib import pyplot as plt
import json
import numpy as np
import seaborn as sns

AU = 1.49597870691e13
PC = 3.085677582e18
d_Orion = 440.0
k_Boltzmann = 1.3806503e-16
cos80 = 0.173648177667

tab = Table.read('../doc/wind-fits.tab', format='ascii.tab')

sources = sorted(set(tab['Source']))
n = len(sources)
colors = sns.color_palette('Set1', n)

# Physical separations in parsec
D = tab['D'] * d_Orion*AU/PC
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.loglog(D, tab['F(star)'], c='k', alpha=0.1, lw=10, label='')
ax2.loglog(D, tab['P(wind)']/k_Boltzmann, c='k', alpha=0.1, lw=10, label='')
mm = (tab['F(ph)/F(*)'] > 0.707) & (tab['F(ph)/F(*)'] < 1.414)
out_colnames = ['Source' ,'D prime', 'R0/D prime', 'Rc/R0 prime full',
                'Rc/R0 prime select', 'beta', 'xi', 'inc', 'D', 'R0/D']

out_formats = {
    'D prime': '%.1f',
    'R0/D prime':         '%.4f',
    'Rc/R0 prime full':   '%.4f',
    'Rc/R0 prime select': '%.4f',
    'beta':               '%.4f',
    'xi':                 '%.4f',
    'inc':                '%.4f',
    'D':                  '%.4f',
    'R0/D':               '%.4f',
}

out_rows = []

def var_range(x, dx):
    return x, dx


for source, color in zip(sources, colors):
    m = tab['Source'] == source
    Dprime = tab['D\''][m][0] * d_Orion*AU/PC
    F = tab['F(photo)'][m][0]

    # Get data from variational fits
    combined_file = '../read-shapes/LV-bowshocks-xyfancy-variations-{}.save'.format(source)
    vardata = json.load(open(combined_file))
    # Combine with data from org table to fill in a new output table 
    Rcp_select = tab['Rc\'/R0\''][m & mm]
    Rcp_full = np.array(vardata['Rc'])/np.array(vardata['R0'])
    beta = tab[r'\beta'][m & mm]
    xi = tab['xi'][m & mm]
    inc = tab['i'][m & mm]
    DD = tab['D'][m & mm] * d_Orion*AU/PC
    R0_D = tab['R0/D'][m & mm]

    out_rows.append({
        'Source': source,
        'D prime': tab['D\''][m][0],
        'R0/D prime':         var_range(np.mean(vardata['R0']), np.std(vardata['R0'])),
        'Rc/R0 prime full':    var_range(np.mean(Rcp_full), np.std(Rcp_full)), 
        'Rc/R0 prime select':  var_range(np.mean(Rcp_select), np.std(Rcp_select)),
        'beta':               var_range(np.mean(beta), np.std(beta)), 
        'xi':                 var_range(np.min(xi), np.max(xi)),
        'inc':                var_range(np.mean(inc), np.std(inc)),
        'D':                  var_range(np.mean(DD), np.std(DD)),
        'R0/D':               var_range(np.mean(R0_D), np.std(R0_D)),
    })
                   
    ax1.loglog([Dprime, Dprime/cos80], [F, F], ':', c=color, alpha=0.4, label='')
    ax1.loglog(D[m], tab['F(photo)'][m],
               '-', c=color, alpha=0.4, label='')
    ax1.loglog(D[m & mm], tab['F(photo)'][m & mm],
               'o-', lw=3, c=color, label=source)
    ax2.loglog(D[m], tab['P(in)'][m]/k_Boltzmann,
               '-', c=color, alpha=0.4, label='')
    ax2.loglog(D[m & mm], tab['P(in)'][m & mm]/k_Boltzmann,
               'o-', c=color, lw=3, label=source)
ax2.legend(ncol=2, loc='lower left')
ax2.set_xlim(0.008, 0.3)
ax2.set_ylim(1e7, 4e9)
ax1.set_ylim(2e12, 3e14)
ax1.set_ylabel(r'Ionizing Flux, $\mathrm{cm^{-2}\ s^{-1}}$')
ax2.set_ylabel(r'Stagnation Pressure: $P/k$, $\mathrm{cm^{-3}\ K}$')
ax2.set_xlabel('Distance, parsec')
fig.set_size_inches(5, 8)
fig.tight_layout()
fig.savefig('plot-wind-fits.pdf')

out_tab = Table(names=out_colnames, rows=out_rows)

print(out_tab.pprint(max_width=-1, max_lines=-1))

out_tab.write('arc-fit-table-for-paper.tab',
              format='ascii.fixed_width', formats=out_formats)


"""

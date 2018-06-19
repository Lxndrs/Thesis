"""Plot data from table of pressures and fluxes from proplyd arc fits
   Original: William Henney.
   Modified Version: Jorge A. Tarango-Yong
"""
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

from astropy.table import Table
from matplotlib import pyplot as plt
import json
import numpy as np
import seaborn as sns


# Useful Physical and Astronomical constants.
AU = 1.49597870691e13
PC = 3.085677582e18
d_Orion = 414.0 # Use Menten et al. (2007) Measurement
k_Boltzmann = 1.3806503e-16
cos80 = 0.173648177667

# Read parameters table
tab = Table.read('../wind-fits.tab', format='ascii.tab')
sources = sorted(set(tab['Fuente']))
n = len(sources)
colors = sns.color_palette('Set1', n)

# Physical separations in parsec
D = tab['D(as)'] * d_Orion*AU/PC
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.loglog(D, tab['F(star)'], c='k', alpha=0.1, lw=10, label='')
ax2.loglog(D, tab['P(wind)']/k_Boltzmann, c='k', alpha=0.1, lw=10, label='')
mm = (tab['F(ph)/F(*)'] > 0.707) & (tab['F(ph)/F(*)'] < 1.414) # Subsamples where photoionization balance is accomplished
nn = (tab['F(ph)/F(*)'] > 0.33) & (tab['F(ph)/F(*)'] < 3.) # Subsamples where photoionization balance is  weakly accomplished.
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
    m = tab['Fuente'] == source
    Dprime = tab['D\'(as)'][m][0] * d_Orion*AU/PC
    F = tab['F(photo)'][m][0]

    # Get data from variational fits
    combined_file = '../saves/LV-bowshocks-xyfancy-variations-{}.save'.format(source)
    vardata = json.load(open(combined_file))
    # Combine with data from org table to fill in a new output table 
    Rcp_select = tab['\Pi\''][m & mm]
    Rcp_full = np.array(vardata['Rc'])/np.array(vardata['R0'])
    beta = tab[r'\beta'][m & mm]
    xi = tab['xi'][m & mm]
    inc = tab['i (deg)'][m & mm]
    DD = tab['D(as)'][m & mm] * d_Orion*AU/PC
    R0_D = tab['q'][m & mm]

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
#        'D':                  var_range(np.mean(DD), np.std(DD)),
#        'R0/D':               var_range(np.mean(R0_D), np.std(R0_D)),
#    })
                   
    ax1.loglog([Dprime, Dprime/cos80], [F, F], ':', c=color, alpha=0.4, label='')
    ax1.loglog(D[m], tab['F(photo)'][m],
               '-', c=color, alpha=0.4, label='')
    ax1.loglog(D[m & mm], tab['F(photo)'][m & mm],
               'o-', lw=3, c=color, label=source)
    ax2.loglog(D[m], tab['P(in)'][m]/k_Boltzmann,
               '-', c=color, alpha=0.4, label='')
    ax2.loglog(D[m & mm], tab['P(in)'][m & mm]/k_Boltzmann,
               'o-', c=color, lw=3, label=source)
    # Now plot the subsamples which accomplish weakly photoionization balance.
    ax1.loglog(D[m], tab['F(photo)'][m],
               '-', c=color, alpha=0.4, label='')
    ax1.loglog(D[m & (nn & ~mm)], tab['F(photo)'][m & (nn & ~mm)],
               'o-', lw=0.5, c=color, label='', alpha=0.4)
    ax2.loglog(D[m], tab['P(in)'][m]/k_Boltzmann,
               '-', c=color, alpha=0.4, label='')
    ax2.loglog(D[m & (nn & ~mm)], tab['P(in)'][m & (nn & ~mm)]/k_Boltzmann,
               'o-', c=color, lw=0.5, label='', alpha=0.4)
ax2.legend(ncol=2, loc='lower left')
ax2.set_xlim(0.008, 0.3)
ax2.set_ylim(1e7, 4e9)
ax1.set_ylim(2e12, 3e14)
ax1.set_ylabel(r'Ionizing Flux, $\mathrm{cm^{-2}\ s^{-1}}$')
ax2.set_ylabel(r'Stagnation Pressure: $P/k$, $\mathrm{cm^{-3}\ K}$')
ax2.set_xlabel('Distance, parsec')
fig.set_size_inches(5, 8)
fig.tight_layout()
fig.savefig('../Figures/plot-wind-fits.pdf')

#out_tab = Table(names=out_colnames, rows=out_rows)

#print(out_tab.pprint(max_width=-1, max_lines=-1))

#out_tab.write('../arc-fit-table-for-paper.tab',
#              format='ascii.fixed_width', formats=out_formats)

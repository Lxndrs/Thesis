
import sys
sys.path.insert(0,"../bowshock-shape/Dust-wave/")
sys.path.insert(0,"../bowshock-shape/")
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import bow_projection as bp
import ancantoid_shape
import bow_diagnostic

#####################################
# Modify program to plot            #
# \Pi' vs inclination               #
# instead of \Lambda' vs \Pi'       #
# Document everyhing as I           #
# understand what each command does #
#####################################

# Maybe I won't need this
#try: 
#    xiset = sys.argv[1] # Additional argument in command line to enter the anisotropy parameter 'xi'
#    plotfile = sys.argv[0].replace('.py', f'-{xiset}.pdf') # The output pdf file will be the name of the program itself with extension
#    assert xiset in 'ab' #Test searching potential errors  # '.pdf' instead of '.py' 
#    istart = -2 if xiset == 'a' else -1 # Honestly I don't know what this means
#except:
#    sys.exit(f"Usage: {sys.argv[0]} a|b") # Exit in case of failure

#sns.set_style('ticks') #Set plot axis style
# Adapt the style to the other graphs I have so far
sns.set_style("white") 
#fig, ax = plt.subplots(figsize=(4, 4)) # set subplot size
f = plt.figure()
ax1 = f.add_subplot(1, 3, 1, adjustable="box") # wilkinoid + cantoid plot
ax2 = f.add_subplot(1, 3, 2, adjustable="box") # Ancantoid xi=0.8 plot
ax3 = f.add_subplot(1, 3, 3, adjustable="box") # Ancantoid xi=0.4 plot

bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03 
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01     #Stuff from bow_projection classes

#left_annotate_pars = dict(xytext=(-5, 5), ha='right', va='bottom')
#right_annotate_pars = dict(xytext=(5, -5), ha='left', va='top') # set location oftext inside plot


#Rc_grid = np.linspace(0.0, 10.0, 2000)
#R90_T0_grid = np.sqrt(2*Rc_grid)
#R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
#R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 #set grids for shaded regions

#ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
#ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1) # shade the different regions in diagram
#ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5) # Plot the parabolic interface
#ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1) #plot horizontal line
#ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1) #plot vertical line
#ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1) #Plot diagonal identity line

XI_LIST = [None, 0.8, 0.4]
BETA_LIST = [0.005, 0.01, 0.05, 0.08, 0.5]
nxi, nbeta = len(XI_LIST), len(BETA_LIST) # set shells parameters (xi=None for cantoid shell)
cols = sns.color_palette('magma', n_colors=nbeta+1) # color palette of curves
# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
#ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)
# And plot the projected wilkinoids 
shape = bp.wilkinoid_R_theta
th_inf = bp.theta_infinity(shape)
inc = np.linspace(0.0, th_inf - np.pi/2, 50)
tab = bow_diagnostic.parameter_table(inc, shape)
Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
ax1.plot(np.degrees(inc), Rc, '-', c=cols[0], label="wilkinoid", lw=2.0, alpha=1.0)
#sini = np.linspace(0.0, 1.0, 20)
#inc_e = np.arcsin(sini)
#tab_e = bow_diagnostic.parameter_table(inc_e, shape)
#Rc_e, R90_e = tab_e['tilde R_c prime'], tab_e['tilde R_90 prime']
#ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
#           linewidths=0.1, edgecolors='none',
#           c='w', alpha=0.5, label="_nolabel_")

#annot_pars_list = [right_annotate_pars]*2 + [left_annotate_pars]*2 
#for beta in BETA_LIST[::-1]:
#    for xi, col, annot_pars in list(zip(XI_LIST, cols, annot_pars_list))[istart::-2]: #start loops in beta and xi
for xi in XI_LIST:
    k = None if xi is None else 2/xi - 2
    for beta, col in zip(BETA_LIST, cols[1:]):    
#        if beta == BETA_LIST[0]:
#            label = "Cantoid" if k is None else fr"Ancantoid $k = {k:.1f}$" # set label into plot
#        else:
#            label = "_nolabel_"
#
        if xi is None: #cantoid case
            shape = bp.Spline_R_theta_from_function(
                ngrid=1000,
                shape_func=bp.cantoid_R_theta,
                shape_func_pars=(beta,))
        else: #ancantoid case
            shape = ancantoid_shape.Ancantoid(xi=xi, beta=beta, n=301)

        th_inf = bp.theta_infinity(shape)
        inc = np.linspace(0.0, th_inf - np.pi/2, 200)
        tab = bow_diagnostic.parameter_table(inc, shape)
        Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
        label = r"$\beta={}$".format(beta)
        if xi is None:
            ax1.plot(np.degrees(inc), Rc, '-', c=col, label=label, lw=1.0, alpha=1.0) #Plot Rc vs i
        elif xi==0.8:
            ax2.plot(np.degrees(inc), Rc, '-', c=col, label=label, lw=1.0, alpha=1.0) #Plot Rc vs i
        else:
            ax3.plot(np.degrees(inc), Rc, '-', c=col, label=label, lw=1.0, alpha=1.0) #Plot Rc vs i
        # Get points evenly spaced in sin i
#        sini = np.linspace(0.0, 1.0, 20)
#        inc_e = np.arcsin(sini)
#        inc_e = inc_e[inc_e < th_inf - np.pi/2]
#        tab_e = bow_diagnostic.parameter_table(inc_e, shape)
#        Rc_e, R90_e = tab_e['tilde R_c prime'], tab_e['tilde R_90 prime']
#        ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
#                   linewidths=0.1, edgecolors='none',
#                   c=col, alpha=0.5, label="_nolabel_")

        # Put a dot at the i=0 case
#        ax.plot(Rc[0:1], R90[0:1], 'o', mec='none', c=col, label="_nolabel_", alpha=0.7)
        # Label the dot with the cross-over inclination
#        beta_label = rf'$\beta = \mathrm{{{beta:g}}}$'
#        if beta_label.endswith('1}$'):
            # But only for some of them
#            ax.annotate(beta_label, xy=(Rc[0], R90[0]),
#                        textcoords='offset points',
#                        fontsize='x-small', color=col, **annot_pars)


ax3.legend(ncol=1, fontsize='small', frameon=True, title=r"Ancantoid $k=3.0$") # legends board
ax3.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.8, 6.0],
#    ylim=[-3.0, 1.1],
#    xlabel=r"inclination (deg)",
#    ylabel=r"Projected planitude: $\Pi'$", #Plot settings
)        

ax1.legend(ncol=1, fontsize='small', frameon=True, title="Isotropic inner wind") # legends board
ax1.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.8, 6.0],
#    ylim=[-3.0, 1.1],
    xlabel=r"inclination (deg)",
    ylabel=r"Projected planitude: $\Pi'$", #Plot settings
)

ax2.legend(ncol=1, fontsize='small', frameon=True, title=r"Ancantoid $k=0.5$") # legends board
ax2.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.8, 6.0],
#    ylim=[-3.0, 1.1],
#    xlabel=r"inclination (deg)",
#    ylabel=r"Projected planitude: $\Pi'$", #Plot settings
)
#sns.despine()
ax1.text(5, 5.8, "( a )")
ax2.text(5, 5.8, "( b )")
ax3.text(5, 5.8, "( c )")
f.tight_layout()
f.set_size_inches(17, 10)
f.savefig("./Figures/Pi-vs-i.pdf")
#print(plotfile, end='')
# The End

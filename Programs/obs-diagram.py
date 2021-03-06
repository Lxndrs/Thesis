
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import seaborn as sns
import sys
sys.path.insert(0,"../bowshock-shape/Dust-wave/")
sys.path.insert(0,"../bowshock-shape/")
import json
import glob
import bow_projection as bp
import ancantoid_shape
import bow_diagnostic
import matplotlib.ticker as mpl

# Set graph style
f = plt.figure()

sns.set_style("ticks")


# Set theoretical curves

bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03 
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01     #Stuff from bow_projection classes

XI_LIST = [None, 0.8, 0.4, 0.2, 0.1]
BETA_LIST = [5e-4, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.06, 0.1]
nxi, nbeta = len(XI_LIST), len(BETA_LIST) # set shells parameters (xi=None for cantoid shell)
cols = sns.color_palette('magma', n_colors=nbeta) # color palette of curves

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

m_savefiles = glob.glob("../saves/LV-bowshocks-xyfancy-positionswill-*.save")
dict_xtext = {"LV2":0.18, "LV2b":0.18, "LV3":0.18, "LV4":0.18, "LV5":0.25, "168-328":0.25, "169-338":0.25, "177-341":0.25, "180-331":0.32}
dict_ytext = {"LV2":0.9, "LV2b":0.85, "LV3":0.8, "LV4":0.75, "LV5":0.9, "168-328":0.85, "169-338":0.8, "177-341":0.75, "180-331":0.9}

for xi in XI_LIST:
    k = None if xi is None else 2/xi - 2
    ax = f.add_subplot(1, 1, 1, adjustable="box") 
    for beta, col in zip(BETA_LIST, cols):    
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
        Rc, R0pR0 = tab['tilde R_c prime'], tab['R_0 prime']
        R0D = np.sqrt(beta)/(1+np.sqrt(beta))
        DDp = 1./np.cos(inc)
        R0 = R0pR0*R0D*DDp
        label = r"$\beta={}$".format(beta)
        ax.plot(R0, Rc, '-', c=col, label=label, lw=1.0, alpha=1.0)
        # Get points evenly spaced every 15 degrees (and minor marks every 5 degrees)
        inc_e = np.radians(np.array([15, 30, 45, 60, 75, 90]))
        inc_e2 = np.radians(np.array([5, 10, 20, 25, 35, 40, 50, 55, 65, 70, 80, 85]))
        inc_e = inc_e[inc_e < th_inf - np.pi/2]
        inc_e2 = inc_e2[inc_e2 < th_inf - np.pi/2]
        tab_e = bow_diagnostic.parameter_table(inc_e, shape)
        tab_e2 = bow_diagnostic.parameter_table(inc_e2, shape)
        Rc_e, R0pR0_e = tab_e['tilde R_c prime'], tab_e['R_0 prime']
        Rc_e2, R0pR0_e2 = tab_e2['tilde R_c prime'], tab_e2['R_0 prime']
        DDp_e = 1./np.cos(inc_e)
        R0_e = R0pR0_e*R0D*DDp_e
        DDp_e2 = 1./np.cos(inc_e2)
        R0_e2 = R0pR0_e2*R0D*DDp_e2
        ax.scatter(R0_e, Rc_e, marker='o', s=3**2,
                   linewidths=0.1, edgecolors='none',
                   c=col, alpha=0.8, label="_nolabel_")
        ax.scatter(R0_e2, Rc_e2, marker='|', s=3**2,
                   linewidths=0.08, edgecolors='none',
                   c=col, alpha=0.5, label="_nolabel_")

        # Put a dot at the i=0 case
        ax.plot(R0[0:1], Rc[0:1], 'o', mec='none', c=col, label="_nolabel_", alpha=0.7)



    #Add the observational points
    for savefile in m_savefiles:
        data = json.load(open(savefile))
        combined_file = savefile.replace('positionswill', 'variations')
        vardata = json.load(open(combined_file))
        ax.plot(data["R0"], data["Rc"]/data["R0"],
               # color=colordict[data["proplyd"]],
               color='k',
               marker="o")
        ax.annotate(data["proplyd"], xy=(data["R0"], data["Rc"]/data["R0"]),
                   xytext=(dict_xtext[data["proplyd"]], dict_ytext[data["proplyd"]]),
                   textcoords="figure fraction", fontsize="xx-small",
                   bbox=dict(boxstyle='round, pad=0.5',
                             fc=colordict[data["proplyd"]],
                             alpha=0.5))
        # Plot the variations of the fits with points removed
        R0_d = data["R0"]
        A = data["Rc"]/data["R0"]
        var_R0 = vardata["R0"]
        var_A = np.array(vardata["Rc"])/np.array(vardata["R0"])
        for vR0, vA in zip(var_R0, var_A):
#        # Scale gives fractional deviation from typical value
            scale = np.hypot((vR0 - R0_d)/0.25, (vA - A)/1.5)
            alpha = 1./(1 + 20.0*scale)
            ax.plot([R0_d, vR0], [A, vA], '-',
                    lw=2, alpha=alpha, color=colordict[data["proplyd"]])
    ktitle = "Cantoid" if k is None else r"$k={}$".format(k)
    filesuffix = "Cantoid" if k is None else "k{:02.0f}".format(10*k)
    ax.legend(loc="upper right", title=ktitle, fontsize="x-small", ncol=2, frameon=True)
    ax.set_xlabel(r"Projected apex radius: $R'_0/D'$")
    ax.set_ylabel(r"Projected Planitude: $\Pi'$")
    ax.get_xaxis().set_minor_locator(mpl.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.AutoMinorLocator())
    ax.grid(b=True, which='major', linewidth=1.0)
    ax.grid(b=True, which='minor', linewidth=0.5)
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 4.0)
    f.set_size_inches(6, 6)
    f.tight_layout()
    f.savefig("../Figures/obs-diagnostic-Pi-R0-{}.pdf".format(filesuffix))
    f.clf()

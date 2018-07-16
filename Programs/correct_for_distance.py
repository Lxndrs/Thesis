
import numpy as np
import json

"""
Using columns 2, 5 and 6 from Henney & Arthur 1998 I deduced that they used
a distance to Theta 1 C Ori of 460 pc, but I'm using a distance of 414 +/- 6.8 (Menten et al. 2007)
for that purpose. So this would make some inconsistencies in the measurements of the derived stagnation pressure
and flux. The purpose of this program is to scale the ionization front radius and the surface brightness for the
most recent distance measurement.
"""

# set observational data from thesis and HA1998

sources = ["168-328", "169-338", "177-341", "180-331", "LV2", "LV2b", "LV3", "LV4", "LV5"]
d_prime = {"168-328":6.8, "169-338":16.4, "177-341":25.6, "180-331":25.1, "LV2":7.8, "LV2b":7.2, "LV3":6.9, "LV4":6.2, "LV5":9.6}
r14_HA = {"168-328":2.8, "169-338":2.8, "177-341":20.4, "180-331":12.2, "LV2":7.9, "LV2b":2.5, "LV3":5.0, "LV4":3.5, "LV5":6.3}
dr14_HA = {"168-328":0.3, "169-338":0.3, "177-341":1.6, "180-331":1.2, "LV2":0.3, "LV2b":0.6, "LV3":0.6, "LV4":0.3, "LV5":0.6}
N6_HA = {"168-328":4.00, "169-338":1.40, "177-341":0.41, "180-331":0.48, "LV2":2.53, "LV2b":4.13, "LV3":3.11, "LV4":4.13, "LV5":2.33}
dN6_HA_P = {"168-328":2.11, "169-338":0.91, "177-341":0.11, "180-331":0.25, "LV2":1.06, "LV2b":1.85, "LV3":1.48, "LV4":1.97, "LV5":1.13}
dN6_HA_m = {"168-328":0.02, "169-338":0.05, "177-341":0.03, "180-331":0.05, "LV2":0.17, "LV2b":0.25, "LV3":0.44, "LV4":0.34, "LV5":0.33}
d_ori_HA = 460
d_ori_new = 414
Delta_d_ori_new = 6.8

# Compute new data
savedata = {"r14":"", "dr14":"", "N6":"", "dN6+":"", "dN6-":""}
r14_new = {}
dr14_new = {}
N6_new = {}
dN6_P_new = {}
dN6_m_new = {}
for source in sources:
    r14 = r14_HA[source]*d_ori_new/d_ori_HA
    dr14 = r14*(Delta_d_ori_new/d_ori_new + dr14_HA[source]/r14_HA[source])
    N6 = N6_HA[source]*np.sqrt(r14_HA[source]/r14)
    N6_new[source] = "{:.2f}".format(N6)
    dN6_P = N6*(dN6_HA_P[source]/N6_HA[source] + 0.5*dr14/r14 + 0.5*dr14_HA[source]/r14_HA[source])
    dN6_m = N6*(dN6_HA_m[source]/N6_HA[source] + 0.5*dr14/r14 + 0.5*dr14_HA[source]/r14_HA[source])
    r14_new[source] = "{:.1f}".format(r14)
    dr14_new[source] = "{:.1f}".format(dr14)
    dN6_P_new[source] = "{:.2f}".format(dN6_P)
    dN6_m_new[source] = "{:.2f}".format(dN6_m)

savedata["r14"] = r14_new
savedata["dr14"] = dr14_new
savedata["N6"] = N6_new
savedata["dN6+"] = dN6_P_new
savedata["dN6-"] = dN6_m_new

savefile = "distances_corrected.json"
with open(savefile, 'w') as f:
    json.dump(savedata, f, indent=2)

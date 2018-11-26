
import numpy as np
import json

"""
Using columns 2, 5 and 6 from Henney & Arthur 1998 I deduced that they used
a distance to Theta 1 C Ori of 460 pc, but I'm using a distance of 414 +/- 6.8 (Menten et al. 2007)
for that purpose. So this would make some inconsistencies in the measurements of the derived stagnation pressure
and flux. The purpose of this program is to scale the ionization front radius and the surface brightness for the
most recent distance measurement.
Update: Use García-Arredondo et al. 2001 for more updated measurements with smaller errors and they mention explicitly
the distance assumed to \thC{} of 430 pc. Most measurements of r14 and n6 remains almost unchanged, which suggests me that this
was the correct answer all the time.
"""

# set observational data from thesis and HA1998

sources = ["168-328", "169-338", "177-341", "180-331", "LV2", "LV2b", "LV3", "LV4", "LV5"]
d_prime = {"168-328":6.8, "169-338":16.4, "177-341":25.6, "180-331":25.1, "LV2":7.8, "LV2b":7.2, "LV3":6.9, "LV4":6.2, "LV5":9.6}
r14_GAH = {"168-328":2.8, "169-338":2.8, "177-341":20.4, "180-331":12.2, "LV2":7.9, "LV2b":2.5, "LV3":5.0, "LV4":3.5, "LV5":6.3}
dr14_GAH = {"168-328":0.3, "169-338":0.3, "177-341":1.6, "180-331":1.2, "LV2":0.3, "LV2b":0.6, "LV3":0.6, "LV4":0.3, "LV5":0.6}
N6_GAH = {"168-328":4.00, "169-338":1.40, "177-341":0.41, "180-331":0.48, "LV2":2.60, "LV2b":4.13, "LV3":3.13, "LV4":4.13, "LV5":2.33}
dN6_GAH = {"168-328":0.01, "169-338":0.91, "177-341":0.02, "180-331":0.03, "LV2":0.11, "LV2b":0.16, "LV3":0.3, "LV4":0.23, "LV5":0.22}
d_ori_GAH = 430
d_ori_new = 414
Delta_d_ori_new = 6.8

# Compute new data
savedata = {"r14":"", "dr14":"", "N6":"", "dN6":""}
r14_new = {}
dr14_new = {}
N6_new = {}
dN6_new = {}

for source in sources:
    r14 = r14_GAH[source]*d_ori_new/d_ori_GAH
    dr14 = r14*(Delta_d_ori_new/d_ori_new + dr14_GAH[source]/r14_GAH[source])
    N6 = N6_GAH[source]*np.sqrt(r14_GAH[source]/r14)
    N6_new[source] = "{:.2f}".format(N6)
    dN6 = N6*(dN6_GAH[source]/N6_GAH[source] + 0.5*dr14/r14 + 0.5*dr14_GAH[source]/r14_GAH[source])
    r14_new[source] = "{:.1f}".format(r14)
    dr14_new[source] = "{:.1f}".format(dr14)
    dN6_new[source] = "{:.2f}".format(dN6)

savedata["r14"] = r14_new
savedata["dr14"] = dr14_new
savedata["N6"] = N6_new
savedata["dN6"] = dN6_new


savefile = "distances_corrected.json"
with open(savefile, 'w') as f:
    json.dump(savedata, f, indent=2)

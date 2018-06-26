
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
d_ori_HA = 460
d_ori_new = 414

# Compute new data
r14_new = {}
for source in sources:
    r14_new[source] = "{:.1f}".format(r14_HA[source]*d_ori_new/d_ori_HA)

savefile = "distances_corrected.json"
with open(savefile, 'w') as f:
    json.dump(r14_new, f, indent=2)

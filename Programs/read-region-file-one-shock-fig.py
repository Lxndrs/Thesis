import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
from scipy import optimize
import seaborn as sns
parser = argparse.ArgumentParser(
    description=""" Choose a region file to work and the angle
    to measure radius""")

parser.add_argument("--region", type=str,
                    default="LV-positions-new-will.reg",
                    help=" Choose a region file to work ")

parser.add_argument('--proplyd', type=str, default='LV3',
                    help='choose a proplyd to work')
parser.add_argument('--on-axis', action="store_true",
                    help='Force circle center to lie on star-proplyd axis')
parser.add_argument('--tfit', type=float, default=45, help='upper theta value for fitting data')

cmd_args = parser.parse_args()
regionfile = cmd_args.region
name, ext = regionfile.split('.')
regfl_chsn = name.split('-')[1] + name.split('-')[-1]
tfit = cmd_args.tfit
sns.set_style("whitegrid")
def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)


#@countcalls
def f_2(c):
    """calculate the algebraic distance between the 2D points and the
    mean circle centered at c=(xc, yc)"""
    Ri = calc_R(*c)
    return Ri - Ri.mean()


def extract_data(line):
    """Find the coordinates and description text from the line of a ds9
    regions file

    """
    coordstring, paramstring = line.split("# ")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0]
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0]
    else:
        text = "NO TEXT"
    return ra, dec, text

Shapes = {}
Centers = {}

with open(regionfile) as f:
    lines = f.readlines()
    for line in lines:
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line
        if skipthisline:
            continue
        ra, dec, text = extract_data(line)
        # name of source (LV knot or star)
        source = text.split()[0]
        if "OW" in text or text == "NO TEXT":
            # We do not use the OW positions
            continue

        shr, smin, ssec = ra.split(':')
        hr, mn, sec = float(shr), float(smin), float(ssec)
        ra_arcsec = 15*(3600.*hr+60.*mn+sec)
        sdeg, samin, sasec = dec.split(':')
        deg, amin, asec = float(sdeg), float(samin), float(sasec)
        dec_arcsec = 3600.*deg + np.sign(deg)*(60*amin + asec)

        if "shock" in text:
            if source not in Shapes:
                Shapes[source] = []
            Shapes[source].append(np.array([ra_arcsec, dec_arcsec]))
        else:
            Centers[source] = np.array([ra_arcsec, dec_arcsec])


# print Centers
# print
# print Shapes

proplyds = Shapes.keys()
proplyd = cmd_args.proplyd


# vector star->proplyd
try:
    vecD = Centers['th1C'] - Centers[proplyd]
except KeyError:
    # no center for this proplyd, so skip it
    print('No center')
# scalar distance star->proplyd
D = np.hypot(*vecD)
# unit vector in direction star->proplyd
xunit = vecD/D
# unit vector perpendicular to star->proplyd
yunit = np.array([-xunit[1], xunit[0]])

print("Proplyd ", proplyd)
print("theta cutoff",tfit)
print("On axis",cmd_args.on_axis)
   # print "D = ", D
   # print "xunit = ", xunit
   # print "yunit = ", yunit

x = []
y = []
R = []
theta = []
for shockpos in Shapes[proplyd]:
    # vector radius proplyd->shock normalized by separation D
    vecR = (shockpos - Centers[proplyd])/D
    R.append(np.hypot(*vecR))
    x.append(np.dot(vecR, xunit))
    y.append(np.dot(vecR, yunit))
    theta.append(np.arctan2(np.dot(vecR, yunit), np.dot(vecR, xunit)))

x = np.array(x)
y = np.array(y)
R = np.array(R)
theta = np.array(theta)

#choosing the desired data to fit
th_fit = theta[np.abs(theta)<=np.radians(tfit)]
R_fit = R[np.abs(theta)<=np.radians(tfit)]
#print len(th_fit),len(R_fit)
#Here is where the fit start
method_2 = "leastsq"
x_m = 0.0
y_m = 0.0
pixel_size = 0.05
# Assume 1-pixel uncertainty in positions
sigma = np.ones(len(R_fit))*pixel_size


def circle_function(theta, xc, yc, rc):
    """Calculate distance from the origin as a function of theta for a
    circle that is centered at xc, yc, with radius rc

    """
    fac1 = xc*np.cos(theta) + yc*np.sin(theta)
    fac2 = xc*np.sin(theta) - yc*np.cos(theta)
    return fac1 + np.sqrt(rc**2 - fac2**2)


def circle_on_axis(theta, xc, rc):
    """Same as circle_function but with yc=0.
    """
    return circle_function(theta, xc, 0.0, rc)


def data_minus_model(thdata, rdata, rmodel_func, model_pars):
    return rdata - rmodel_func(thdata, *model_pars)

R_m = calc_R(x_m, y_m).mean()
if cmd_args.on_axis:
    # Two-parameter version
    p0 = x_m, R_m
    f = circle_on_axis
else:
    # Three-parameter version
    p0 = x_m, y_m, R_m
    f = circle_function

print("Initial parameter:")
print("Circle x, y, r = ", x_m, y_m, R_m)
popt, pcov = optimize.curve_fit(f, th_fit, R_fit, p0, sigma)

if cmd_args.on_axis:
    xc_2, R_2 = popt
    yc_2 = 0.0
else:
    xc_2, yc_2, R_2 = popt

print("Fit results:")
print("Circle x, y, r = ", xc_2, yc_2, R_2)
print("Covariance matrix: ", pcov)

# Deviation betwen model and data
error = data_minus_model(theta, R, f, popt)
print("Deviations between model and data:")
print("Max, rms: ", abs(error).max(), np.sqrt(np.mean(error**2)))

theta_fit = np.linspace(-np.pi, np.pi, 180)
x_fit2 = xc_2 + R_2*np.cos(theta_fit)
y_fit2 = yc_2 + R_2*np.sin(theta_fit)


#***************************************************************************
#We need two points to draw the R0 and Rc. Better choose R0 always on axis
xrcline, yrcline = xc_2 + R_2/np.sqrt(2.), yc_2 + R_2/np.sqrt(2.)
delta = np.arctan2(yc_2, xc_2)
xr0line, yr0line = xc_2 + R_2, 0
#***************************************************************************
#Plotting data (normalized with R0)
R0 = np.sqrt(xr0line**2+yr0line**2)
plt.plot(x/R0, y/R0, "o", label=r"{}: $D = {:.2f}''$".format(proplyd, D))
#***************************************************************************
#Plotting the circle fit
plt.plot(x_fit2/R0, y_fit2/R0, 'k--', label="Circle fit to data", lw=2)
plt.plot([xc_2/R0], [yc_2/R0], 'bx')#Drawing R0 and Rc:
r0x,r0y = np.array([0, xr0line]), np.array([0, yr0line])
rcx,rcy = np.array([xc_2, xrcline]), np.array([yc_2, yrcline])
plt.plot(r0x/R0, r0y/R0, 'g-', label=r"$R'_0/D'={:.3f}$".format(R0))
plt.plot(rcx/R0, rcy/R0, 'm-', label=r"$\Pi'={:.3f}$".format(R_2/R0))
#***************************************************************************


plt.annotate(r"$R'_0$", xy=(0.5*xr0line/R0, 0.5*yr0line/R0), xytext=(20, -20), fontsize='x-small',
    alpha=1.0, textcoords='offset points', ha ='right', va='bottom',
    bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.5),
    arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
plt.annotate(r"$\Pi'$", xy=(0.5*(xc_2+xrcline)/R0, 0.5*(yc_2+yrcline)/R0), xytext=(-20, 20),
    fontsize='x-small', alpha=1.0,textcoords='offset points', ha='left', va='top',
     bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.5),
     arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
# Proplyd position (at the origin in this frame)
plt.plot(0.0, 0.0, "rx", label=None)

# xmin, xmax = plt.xlim()
# ymin, ymax = plt.ylim()

# vmin = min(xmin, ymin)
# vmax = max(xmax, ymax)

# plt.xlim(xmin=vmin, xmax=vmax)
# plt.ylim(ymin=vmin, ymax=vmax)


plt.legend(loc="lower left", frameon=True)
plt.xlabel(r"$z'/R'_0$", fontsize = "large")
plt.ylabel(r"$r'/R'_0$", fontsize = "large")

#plt.grid()
#plt.title("{} fit circle".format(proplyd))

plt.xlim(-3.0, 3.0)
plt.ylim(-2.0, 2.0)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

if cmd_args.on_axis:
    plotfile="LV-bowshocks-xyfancy-onaxis-{}-{}.pdf".format(regfl_chsn,proplyd)
else:
    plotfile="LV-bowshocks-xyfancy-{}-{}.pdf".format(regfl_chsn, proplyd)
plt.savefig(plotfile)



# Write data to a save file
savefile = plotfile.replace('.pdf', '.save')
savedict = {
    'proplyd': proplyd,
    'D': D,
    'R0': R0,
    'Rc': R_2
}

with open(savefile, 'w') as f:
    json.dump(savedict, f)


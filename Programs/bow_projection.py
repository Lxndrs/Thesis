import sys
import numpy as np
from scipy.optimize import brentq
from scipy.misc import derivative


# * Module parameters
#
# The delta theta that is used in the central difference approximation
# to the derivative of the R(theta) function.  For optimum balance
# between round-off and discretization error, this should be of order
# ~sqrt(eps)~, where ~eps~ is the machine precision
DX_FOR_NUMERICAL_DERIVATIVE = 3.0*np.sqrt(np.finfo(1.0).resolution)

# If True, then print out some diagnostic info
DEBUG = False

# * Functions to find plane-of-sky shape
#
# All these functions should have argument lists of the form:
#
# :    theta, [inc], func_R, *args_for_func_R
#
# where ~func_R~ has signature ~func_R(theta, *args_for_func_R)~ and
# ~inc~ is the inclination (for those functions that depend on that).
#
# They should also be written as element-wise functions of a vector
# ~theta~, so no ~if~ statements are allowed, but ~inc~ must be a
# scalar, as must all of the extra args for ~func_R~.
#
def omega(theta, func_R, *args_for_func_R):
    """Find omega = R^{-1} d R/d theta 

Note that theta may be an array. Any extra arguments are passed to
`func_R` after `theta`

    """
    def log_R(theta, *args):
        return np.log(func_R(theta, *args))

    return derivative(log_R, theta,
                      dx=DX_FOR_NUMERICAL_DERIVATIVE, args=args_for_func_R)


def alpha(theta, func_R, *args_for_func_R):
    """Find alpha = tan^{-1} |dy/dx|, the slope angle

Note that theta may be an array. Any extra arguments are passed on to omega

    """
    om = omega(theta, func_R, *args_for_func_R)
    tan_theta = np.tan(theta)
    return np.arctan((1 + om*tan_theta)/(tan_theta - om))


def sin_phi_t(theta, inc, func_R, *args_for_func_R):
    """Returns sin(phi_t), where phi_t is azimuth along tangent line"""
    if np.tan(inc) == 0.0:
        # Avoid NaNs in the zero inclination case
        return np.zeros_like(theta)

    om = omega(theta, func_R, *args_for_func_R)
    tan_theta = np.tan(theta)
    return np.tan(inc)*(1.0 + om*tan_theta)/(om - tan_theta) 


def xyprime_t(theta, inc, func_R, *args_for_func_R):
    """Returns observer-frame (x', y') coordinates of tangent line"""
    R = func_R(theta, *args_for_func_R)
    sphi_t = sin_phi_t(theta, inc, func_R, *args_for_func_R)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    xx = cos_theta*np.cos(inc) - sin_theta*sphi_t*np.sin(inc)
    with np.errstate(all='ignore'):
        yy = sin_theta*np.sqrt(1.0 - sphi_t**2)
    R[R < 0.0] = np.nan
    return R*xx, R*yy


def radius_of_curvature(theta, func_R, *args_for_func_R):
    """Returns R_c = (R^2 + R'^2)^{3/2} / |R^2 + 2 R'^2 - R R''| 

Uses numerical differentiation.  NOT RECOMMENDED SINCE NOT ACCURATE ON
THE AXIS.  Use `axis_Rc` instead.

    """
    R = func_R(theta, *args_for_func_R)
    dR = derivative(func_R, theta,
                    dx=DX_FOR_NUMERICAL_DERIVATIVE, args=args_for_func_R)
    d2R = derivative(func_R, theta, n=2,
                     dx=DX_FOR_NUMERICAL_DERIVATIVE, args=args_for_func_R)
    return (R**2 + dR**2)**1.5 / np.abs(R**2 + 2*dR**2 - R*d2R)


# * Projected R_c and R_{90}
#

# How close we try to get to the asymptotic theta
TOLERANCE_THETA_INFINITY = 1.e-6

def theta_infinity(func_R, *args_for_func_R):
    """Return maximum theta where R its derivative are still finite"""
    th0, dth = 0.0, np.pi
    thinf_m = np.pi - TOLERANCE_THETA_INFINITY
    with np.errstate(all='ignore'):
        # Keep repeating on a finer and finer grid until we get to within
        # the required tolerance
        while dth > TOLERANCE_THETA_INFINITY:
            # This will divide dth by 50 on each iteration
            th, dth = np.linspace(th0, min(thinf_m, th0 + dth), retstep=True)
            # It is more stringent to insist that omega must be
            # finite, since that needs to be an extra distance (=
            # DX_FOR_NUMERICAL_DERIVATIVE) away from the asymptote in
            # order to calculate the numerical derivative
            om = omega(th, func_R, *args_for_func_R)
            # The largest th for which omega is finite must be within at most
            # (dth + DX_FOR_NUMERICAL_DERIVATIVE) of the true asymptote
            if np.isfinite(om).sum() > 0:
                th0 = th[np.isfinite(om)].max()

    return th0



def theta_0_90(inc, func_R, *args_for_func_R):
    """Return (theta_0, theta_90), corresponding projected x and y axes
    """

    # Easy case first
    if inc == 0.0:
        return 0.0, np.pi/2

    # Second, check tangent line existence
    th_inf = theta_infinity(func_R, *args_for_func_R)
    if np.abs(inc) > th_inf - np.pi/2:
        return np.nan, np.nan

    # Otherwise, use root finding

    tani = np.tan(inc)
    sinsq2i = np.sin(2*inc)**2
    cossqi = np.cos(inc)**2

    def _f0(theta):
        """Function that is zero at theta = theta_0"""
        om = omega(theta, func_R, *args_for_func_R)
        return np.sin(theta)*(1.0 - om*tani) - np.cos(theta)*(om + tani)

    def _f90(theta):
        """Function that is zero at theta = theta_90"""
        om = omega(theta, func_R, *args_for_func_R)
        return (1.0/np.tan(theta)
                - (1.0 - np.sqrt(1.0 + om**2 * sinsq2i))/(2*om*cossqi))

    th_inf = theta_infinity(func_R, *args_for_func_R)
    # If the inclination is too high, then there may be no solution
    if np.abs(inc) > th_inf - np.pi/2:
        return np.nan, np.nan

    # The theta_0 value could be anywhere in range 0 -> th_inf, but we
    # set the lower limit higher than that to avoid some rare errors.
    # This should be alright unless R_c/R_0 < 0.1, which is not true
    # for any of the models I am interested in
    th1, th2 = 0.001*inc, th_inf
    # Make sure we do indeed bracket the root
    if _f0(th1)*_f0(th2) <= 0.0:
        # And use Brent's method to find the root
        th0 = brentq(_f0, th1, th2)
    else:
        th0 = np.nan

    # Repeat for the theta_90 value, which must be higher than theta_0
    th1, th2 = 1.001*th0, th_inf
    if _f90(th1)*_f90(th2) <= 0.0:
        th90 = brentq(_f90, th1, th2)
    else:
        th90 = np.nan

    return th0, th90

# Maximum angle used in fitting Rc'
DELTA_THETA_RC = np.radians(30.0)
# Degree of theta'^2 polynomial used in fitting R'(theta') @ theta' = 0
DEGREE_POLY_NEIGHBORHOOD = 1

# Maximum angle used in fitting R90'
DELTA_THETA_R90 = np.radians(5.0)
# Same for around theta' = 90.  These need to be different since we
# are fitting a polynomial in just theta' instead of theta'^2
DEGREE_POLY_NEIGHBORHOOD_90 = 2

def characteristic_radii_projected_new(inc, func_R, *args_for_func_R):
    """Return all the characteristic radii for projected bow shock

    Returns dict of 'R_0 prime', 'tilde R_c prime', 'theta_0', 'tilde
    R_90 prime', 'theta_90'.  This is a new implementation, since I am
    getting fed up with the previous implementation (see
    'characteristic_radii_projected').  Tha main difference is
    that we define the neighborhood in terms of th-prime instead of
    th.
    """
    # Zeroth step, check that we do have a tangent line
    th_inf = theta_infinity(func_R, *args_for_func_R)

    # What to return when there is no solution 
    no_solution = {'R_0 prime': np.nan, 'theta_inf': th_inf,
                   'tilde R_c prime': np.nan, 'theta_0': np.nan,
                   'tilde R_90 prime': np.nan, 'theta_90': np.nan}

    if np.abs(inc) > th_inf - np.pi/2:
        # No tangent line, so return all NaNs
        return no_solution

    th0, th90 = theta_0_90(inc, func_R, *args_for_func_R)
    th = np.linspace(th0, th_inf, 1001)
    xp, yp = xyprime_t(th, inc, func_R, *args_for_func_R)
    m = np.isfinite(xp) & np.isfinite(yp)
    # Reflect to give full projected bow
    xxp = np.concatenate((xp[m][::-1], xp[m]))
    yyp = np.concatenate((-yp[m][::-1], yp[m]))

    # Polar coordinates on plane of sky
    Rp = np.hypot(xxp, yyp)
    thp = np.arctan2(yyp, xxp)

    # Now fit R'(th') around th' = 0 to determine R_0' and R_c'
    m = np.isfinite(Rp*thp) & (np.abs(thp) < DELTA_THETA_RC)
    if m.sum() <= 3*DEGREE_POLY_NEIGHBORHOOD:
        # If not enough good points, then give up
        return no_solution

    # Fit R' with a polynomial in (theta')^2, and use the constant and
    # linear coefficients to find the projected R_0 and R_c
    coeffs = np.polyfit(thp[m]**2, Rp[m],
                        deg=DEGREE_POLY_NEIGHBORHOOD)
    R0_prime = coeffs[-1]
    gamma = coeffs[-2]/coeffs[-1]
    Rc_prime = 1./(1. - 2*gamma)

    # Next, R90'
    m = np.isfinite(Rp*thp) & (np.abs(thp - 0.5*np.pi) < DELTA_THETA_R90)

    # **************** 29 Nov 2017: TO FINISH ****************
    pass

# Number of neighborhood points to use when fitting around projected
# axes (theta' = 0 and theta' = 90)
N_NEIGHBORHOOD = 50
# Theta scale of neighborhood around theta' = 0 in units of (th90 - th0)
SCALE_NEIGHBORHOOD = 0.2
SCALE_NEIGHBORHOOD_90 = 0.03

def characteristic_radii_projected(inc, func_R, *args_for_func_R):
    """Return all the characteristic radii for projected bow shock

Returns dict of 'R_0 prime', 'tilde R_c prime', 'theta_0', 'tilde R_90
prime', 'theta_90'

    """

    # Zeroth step, check that we do have a tangent line
    th_inf = theta_infinity(func_R, *args_for_func_R)

    # What to return when there is no solution 
    no_solution = {'R_0 prime': np.nan, 'theta_inf': th_inf,
                   'tilde R_c prime': np.nan, 'theta_0': np.nan,
                   'tilde R_90 prime': np.nan, 'theta_90': np.nan}

    if np.abs(inc) > th_inf - np.pi/2:
        # No tangent line, so return all NaNs
        return no_solution

    # First, the quantities at th0, which is theta on the projected
    # symmetry axis (y' = 0) for this inclination
    th0, th90 = theta_0_90(inc, func_R, *args_for_func_R)

    # Make a grid of theta in the neighborhood of th0
    #
    # Earlier, I was using (pi - th0), which I am now (29 Nov 2017)
    # changing to (th90 - th0), so I multiply by 2 to maintain the
    # same scale in the case of i=0.  So if we want dtheta interval,
    # we need to put SCALE_NEIGHBORHOOD = dtheta/180, so 0.33 is 60
    # degrees
    dth = 2*SCALE_NEIGHBORHOOD*(th90 - th0)
    th = np.linspace(th0, th0 + dth, N_NEIGHBORHOOD)
    if DEBUG:
        print("theta", th, file=sys.stderr)
        print("R", func_R(th, *args_for_func_R), file=sys.stderr)
        print("sin(phi_t)", sin_phi_t(th, inc, func_R, *args_for_func_R),
              file=sys.stderr)

    # Now find the tangent line and convert back to polar coordinates
    xprime, yprime = xyprime_t(th, inc, func_R, *args_for_func_R)
    Rprime = np.hypot(xprime, yprime)
    thprime = np.arctan2(yprime, xprime)
    if DEBUG:
        print("x'", xprime, file=sys.stderr)
        print("y'", yprime, file=sys.stderr)
        print("R'", Rprime, file=sys.stderr)
        print("theta'", thprime, file=sys.stderr)
    # Filter out any NaNs in the projected coordinates
    m = np.isfinite(Rprime*thprime)
    if m.sum() <= 3*DEGREE_POLY_NEIGHBORHOOD:
        # If not enough good points, then give up
        return no_solution

    # Fit R' with a cubic in (theta')^2, and use the constant and
    # linear coefficients to find the projected R_0 and R_c
    #
    # It seems to be enough to use deg=2 on 8 points
    coeffs = np.polyfit(thprime[m]**2, Rprime[m],
                        deg=DEGREE_POLY_NEIGHBORHOOD)
    R0_prime = coeffs[-1]
    gamma = coeffs[-2]/coeffs[-1]
    Rc_prime = 1./(1. - 2*gamma)
    if DEBUG:
        print("Polynomial coefficients", coeffs/coeffs[-1], file=sys.stderr)


    # Second, the quantities at th90, which is the theta on the projected
    # perpendicular axis (x' = 0)
    dth = SCALE_NEIGHBORHOOD_90*np.pi/2
    th = np.linspace(th90 - dth/2, th90 + dth/2, N_NEIGHBORHOOD)
    xprime, yprime = xyprime_t(th, inc, func_R, *args_for_func_R)
    Rprime = np.hypot(xprime, yprime)
    thprime = np.arctan2(yprime, xprime)
    if DEBUG:
        print("90 x'", xprime, file=sys.stderr)
        print("90 y'", yprime, file=sys.stderr)
        print("90 R'", Rprime, file=sys.stderr)
        print("90 theta'", thprime, file=sys.stderr)
    m = np.isfinite(Rprime*thprime)
    if m.sum() <= 3*DEGREE_POLY_NEIGHBORHOOD_90:
        # If not enough good points, then give up
        return no_solution
    # Fit a polynomial to thprime, Rprime in the vecinity of pi/2
    p = np.poly1d(np.polyfit(thprime[m], Rprime[m],
                             deg=DEGREE_POLY_NEIGHBORHOOD_90))
    # Evaluate polynomial at pi/2 to find R90_prime
    R90_prime = p(np.pi/2)/R0_prime
    if DEBUG:
        print("90 Polynomial coefficients", coeffs/coeffs[-1], file=sys.stderr)

    return {'R_0 prime': R0_prime, 'theta_inf': th_inf,
            'tilde R_c prime': Rc_prime, 'theta_0': th0,
            'tilde R_90 prime': R90_prime, 'theta_90': th90}



# * Example analytic shape functions
#

def wilkinoid_R_theta(theta):
    """Wilkin solution for wind-stream interaction"""
    # Convert to array if scalar
    theta = np.atleast_1d(theta)
    # Equation (9) of Wilkin (1996)
    R = np.sqrt(3*(1.0 - theta/np.tan(theta)))/np.sin(np.abs(theta))
    # Equation is not stable for very small theta, so we use a Taylor
    # expansion instead
    small_angle = np.abs(theta) < 1e-5
    R[small_angle] = 1.0 + 0.2*theta[small_angle]**2
    return R

def cantoid_R_theta(theta, beta):
    """Cantoid solution from CRW for wind-wind interaction

Returns R(theta), normalized to the stagnation radius. Extra parameter
`beta` is the momentum ratio of the two winds.  Note that this will
not be accurate if beta is too close to zero, but it seems to be OK
with beta >= 1.e-6.  For lower values than this, the results will be
almost identical to the Wilkinoid, so `wilkinoid_R_theta` should be
used instead.

    """

    theta = np.atleast_1d(theta)

    # Approximate solution for theta_1, the polar angle measured from
    # the "other" star
    theta1 = np.sqrt(7.5*(-1.0 + np.sqrt(
        1.0 + 0.8*beta*(1.0 - theta/np.tan(theta)))))
    # Make sure theta1 and theta have the same sign
    theta1 *= np.sign(theta)


    # On-axis (theta = 0) radius to stagnation point, in units of
    # star-star separation D
    R0 = np.sqrt(beta)/(1.0+np.sqrt(beta))

    R = np.where(np.abs(theta + theta1) > np.pi,
                 # theta > theta_inf => no solution
                 np.nan,
                 # theta <= theta_inf => return radius in units of R0
                 np.sin(theta1) / np.sin(theta + theta1) / R0)

    # Replace with Taylor expansion close to axis
    C = (1.0 - beta)/30.0
    gamma = C/(1.0 + np.sqrt(beta)) + (1.0 + 2*np.sqrt(beta))/6
    small_angle = np.abs(theta) < 1.e-5
    R[small_angle] = 1.0 + gamma*theta[small_angle]**2

    return R


def paraboloid_R_theta(theta):
    """This is the special parabola with Rc = 2"""
    return 2.0 / (1.0 + np.cos(theta))


def paraboloid_omega_true(theta):
    """Analytic omega for special parabola"""
    return np.sin(theta)  / (1.0 + np.cos(theta))


# * Non-analytic shape functions
#
# These are bow shock shapes for which it is "non-trivial" to
# calculate each R(theta).  E.g., requiring numerical root finding, so
# hard to write naturally in an element-wise vector form
#
# For efficiency, we therefore calculate R(theta) once on a grid, and
# then use a spline interpolation for fast subsequent evaluation
# of R(theta) and its derivative

import scipy.interpolate

class _Spline_R_theta(object):
    """Base class for non-analytic shapes

The R(theta) shape is initialized once on a grid when the class is
instantiated, and fitted by a B-spline.  The object can then be called
as a function of theta, which will be very fast since it just
evaluates the B-spline.

    """

    thgrid = None
    Rgrid = None
    def _init_grid(self, ngrid, **shape_kwds):
        raise NotImplementedError("Override this method in a sub-class")

    def _init_spline(self, kspline, Rmax, smooth):
        """Fit a smoothing spline to the R(theta) grid. 

We fit B-splines to the parametric [x(theta), y(theta)] representation
of the bow shock. `kspline` is the order of the splines (default:
cubic). `Rmax` is the maximum radius to be included in the spline fit.
`smooth` is the spline smoothing condition (see docs for
`scipy.interpolate.splprep`).

        """
        mgood = np.isfinite(self.Rgrid) & (self.Rgrid <= Rmax)
        x = self.Rgrid[mgood]*np.cos(self.thgrid[mgood])
        y = self.Rgrid[mgood]*np.sin(self.thgrid[mgood])
        self.spline_tck, u = scipy.interpolate.splprep(
            [x, y], u=self.thgrid[mgood], k=kspline, s=smooth)

    def __call__(self, theta):
        """Evaluate R(theta) from spline interpolation"""
        theta = np.atleast_1d(theta)
        x, y = scipy.interpolate.splev(theta, self.spline_tck, ext=1)
        # The ext=1 option to splev return 0 for points outside range of theta
        R = np.hypot(x, y)
        # Then we set those out-of-bound points to NaN
        R[(x == 0.0) & (y == 0.0)] = np.nan
        return R

    def __init__(self, ngrid=100, kspline=3, Rmax=100, smooth=0, **shape_kwds):
        """"""
        # Set up grids of theta and R
        self._init_grid(ngrid, **shape_kwds)
        # Set up spline interpolation
        self._init_spline(kspline, Rmax, smooth)


class Spline_R_theta_from_function(_Spline_R_theta):
    """Spline-interpolated bow shock shape from explicit function

Extra parameters for initialization: `shape_func` and
`shape_func_pars`. THIS IS FOR TESTING ONLY!!! It checks that the
interpolation machinery works for simple shapes. Outside of such
tests, there is really no need to use the spline interpolation for
these cases.

    """

    def _init_grid(self, ngrid,
                   shape_func=paraboloid_R_theta,
                   shape_func_pars=()):
        # Include the negative branch so the spline will have the
        # right gradient on the axis
        self.th_inf = theta_infinity(shape_func, *shape_func_pars)
        self.thgrid = np.linspace(0.0, self.th_inf, ngrid)
        self.Rgrid = shape_func(self.thgrid, *shape_func_pars)


class Spline_R_theta_from_grid(_Spline_R_theta):
    """Spline-interpolated bow shock shape from user-specified arrays

Extra parameters for initialization: `theta_grid` and `R_grid`.  This
is the main way that the spline fits will be used.

    """
    def _init_grid(self, ngrid, theta_grid=None, R_grid=None):
        # Note that ngrid is ignored in this implementation
        if theta_grid is not None and R_grid is not None:
            self.thgrid = theta_grid
            self.Rgrid = R_grid
        else:
            raise ValueError("Both theta_grid and R_grid must be specified")


# * Basic tests of functionality
#

if __name__ == "__main__":
    import sys
    from matplotlib import pyplot as plt
    import seaborn as sns

    lib_name = sys.argv[0].replace('.py', '')

    sns.set_style('ticks')
    fig, ax = plt.subplots()

    th = np.linspace(0.0, np.pi, 501)
    th_dg = np.degrees(th)
    ax.plot(th_dg, omega(th, paraboloid_R_theta),
            label="paraboloid")
    ax.plot(th_dg, omega(th, wilkinoid_R_theta),
            label="wilkinoid")
    for beta in 0.001, 0.01, 0.1:
        ax.plot(th_dg, omega(th, cantoid_R_theta, beta),
                label=fr"cantoid $\beta = {beta:.3f}$")
    ax.legend(title=r"Analytic $R(\theta)$ functions")
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='white', lw=4, zorder=100)
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='k', lw=1, ls=':', zorder=101)
    ax.axhspan(0.0, 1.0, color='k', alpha=0.05, ec='none')
    ax.set_yscale('symlog', linthreshy=1.0, linscaley=0.5)
    ax.set(
        xlabel=r"Polar angle: $\theta$, degrees",
        ylabel=r"$\omega \equiv R^{-1} d R / d \theta$",
        xlim=[0, 180],
        ylim=[0.0, 2e2],
        xticks=[0, 30, 60, 90, 120, 150, 180],
    )
    sns.despine()
    fig.tight_layout()
    figfile = f"test_{lib_name}_omega.pdf"
    fig.savefig(figfile)
    print(figfile, end='')

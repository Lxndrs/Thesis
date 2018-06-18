"""Module for solving Equation 6 of Canto, Raga, & Wilkin 1996 (CRW96)

This finds the radius of a stellar wind interaction bowshock in terms
of the momentum and angular momentum injected by the two winds.

CRW96 concentrated on the case of isotropic winds, but this module
will work generally with any cylindrically symmetric wind

All positions and distances are in units of D, the interstar
separation

The philosophy of this module is to be as straightforward and
transparent an implementation as possible of the equations in the
paper.  There is no attempt to be efficient.

"""
from __future__ import print_function
import numpy as np
import scipy.integrate
import scipy.optimize
from scipy.special import gamma as gamma_func

DEBUG_LEVEL = 0

###
### Public functions 
###
def isotropic_momentum(theta):
    """Momentum as a function of angle for an isotropic wind"""
    return 1.0

DIFFUSE_BETA = 0.0              # Parameter giving relative strength of diffuse field
def proplyd_momentum(theta): 
    """Momentum as a function of angle for a proplyd wind

    Proportional to sqrt(cos(theta)) in the head (theta < pi/2), and
    then a small constant value (zero by default) in the tail (theta >
    pi/2).  The tail value is set via the module-level variable
    DIFFUSE_BETA.

    """
    return DIFFUSE_BETA + (1.0 - DIFFUSE_BETA)*np.sqrt(max(0.0,np.cos(theta)))
    
global MOMENTUM_K
MOMENTUM_K = 1.0
def anisotropic_momentum(theta):
    """Momentum as a power law in cos(theta) on the forward hemisphere only"""
    if theta <= np.pi/2:
        return np.cos(theta)**MOMENTUM_K
    else:
        return 0.0

###
### Public classes
###

class Wind(object):
    """Class to represent a stellar wind (or proplyd, etc)

    axial_momentum_flux is the momentum flux at theta=0 (or mu=1) in
    arbitrary units

    momentum_law is a function that describes how the momentum flux
    per unit solid angle varies with theta

    origin is a flag that is True if the center of this wind is at the
    co-ordinate origin.  If origin is False, then the center of the
    wind is at unit distance from the origin along the z axis.

    """
    def __init__(self, axial_momentum_flux=1.0, momentum_law=isotropic_momentum, origin=True):
        self.axial_momentum_flux = axial_momentum_flux
        self.momentum_law = momentum_law
        self.origin = origin



    def Jdot(self, theta):
        """Angular momentum injection rate about the origin,
        integrated between axis and theta

        """
        if self.origin:
            return 0.0
        else:
            if self.momentum_law == isotropic_momentum:
                return 0.25*self.axial_momentum_flux*(theta - np.sin(theta)*np.cos(theta))
            else:
                # I haven't implemented the numerical integration yet
                # in this case, but hopefully we will not need it
                raise NotImplementedError


    def Pidot_z(self, theta):
        """Linear z-momentum injection rate, integrated between axis
        and theta

        """
        if self.momentum_law == isotropic_momentum:
            # Analytic solution for isotropic case
            Pdz = 0.25*np.sin(theta)**2
        else:
            # Numerical integration for the general case
            Pdz, err = scipy.integrate.quad(self._integrand_Pdz, 0.0, theta)
        if self.origin:
            return Pdz*self.axial_momentum_flux
        else:
            # The second star has oppositely directed axial momentum
            return -Pdz*self.axial_momentum_flux

    def Pidot_r(self, theta):
        """Linear r-momentum injection rate, integrated between axis
        and theta

        """
        if self.momentum_law == isotropic_momentum:
            # Analytic solution for isotropic case
            Pdr = 0.25*(theta - np.sin(theta)*np.cos(theta))
        else:
            # Numerical integration for the general case
            Pdr, err = scipy.integrate.quad(self._integrand_Pdr, 0.0, theta)
        return Pdr*self.axial_momentum_flux

    def _integrand_Pdz(self, t):
        return 0.5*np.cos(t)*self.momentum_law(t)*np.sin(t)
    def _integrand_Pdr(self, t):
        return 0.5*np.sin(t)*self.momentum_law(t)*np.sin(t)

    
                            
class BaseShell(object):
    """Class to represent a two-wind interaction shell"""

    def __init__(self, w, w1):
        """The arguments w and w1 should be instances of the class
        Wind()

        The inner wind, w, should have origin=True, while the outer
        wind, w1, should have origin=False

        See the Shell() class for an easier to use wrapper around this
        class

        """
        self.w = w              # "inner" wind
        self.w1 = w1            # "outer" wind

        # We save the values of theta and theta1, so we can use them
        # to find an initial estimate of theta1 for the next angle
        # theta
        self.th1_save = None
        self.th_save = None

        # Pre-calculate the on-axis radius of the shell
        self.beta = self.w.axial_momentum_flux / self.w1.axial_momentum_flux
        self.R0 = np.sqrt(self.beta)/(1.0 + np.sqrt(self.beta))

    def radius(self, theta, method='brent', full=False):
        """Find the spherical radius of the shell as a function of angle

        Should work with scalar or vector argument `theta`. Returns
        `radius`, but if positional argument `full` is `True`, then
        return tuple: `radius`, `theta1`

        """
        def _radius(theta):
            """Helper function to find the shell radius for a single angle, theta"""
            if theta == 0.0:
                # special treatment for the axis
                return self.R0
            elif theta >= self.th_infty:
                # Indicate that we have gone too far
                return -1.0
            else:
                if method == 'fsolve':
                    if self.th1_save is None:
                        # For the first off-axis angle, we use the fact
                        # that R0 tan(theta) ~= (1 - R0) tan(theta1) for
                        # small theta
                        th1_guess = theta*self.R0 / (1.0 - self.R0)
                    else:
                        # For subsequent angles, we do geometric extrapolation
                        th1_guess = self.th1_save*theta/self.th_save 
                    # The tricky bit here is getting th1_guess to be close
                    # enough to the true solution.  If it is not, then the
                    # solver will fail
                    theta1 = _solve_for_th1(self.w, self.w1, theta,
                                            th1_guess, method=method)
                else:
                    # Assume other methods require root to be bracketed
                    # Must be between 0 and th1_infty
                    if self.th1_save is None:
                        a, b = 1e-10, self.th1_infty
                    else:
                        a, b = self.th1_save, self.th1_infty

                    theta1 = _solve_for_th1(self.w, self.w1, theta,
                                            bounds=[a, b], method=method)
                    if DEBUG_LEVEL > 0:
                        print('+++', self.th_infty - theta, self.th1_infty - theta1)
                self.th_save = theta
                self.th1_save = theta1
                return _radius_eq23(theta, theta1)

        try:
            # case where theta is iterable
            rslt = np.empty_like(theta)
            th1_rslt = np.empty_like(theta)
            for i, t in enumerate(theta):
                r = _radius(t)
                if r > 0.0:
                    rslt[i] = r
                    th1_rslt[i] = self.th1_save
                else:
                    # assume we have got to th_max
                    # so fill the remainder with NaNs
                    rslt[i:] = np.nan
                    th1_rslt[i:] = np.nan
                    break
            if full:
                return rslt, th1_rslt
            else:
                return rslt
        except TypeError:
            # fall-over case where theta is scalar
            if full:
                return _radius(theta), self.th1_save
            else:
                return _radius(theta)



class Shell(BaseShell):
    """Easy-to-use class to represent a two-wind interaction shell"""
    def __init__(self, beta=1.0, innertype="isotropic", outertype="isotropic", xi=None):
        """Parameters: 
        beta: axial momentum flux ratio (inner/outer)
        innertype: either 'proplyd' or 'isotropic'
        outertype: must be 'isotropic'
        """
        global MOMENTUM_K 
        if innertype == "anisotropic":
            mlaw = anisotropic_momentum
            # xi = 2 / (k + 2) => k = (2/xi) - 2 
            MOMENTUM_K = 2.0/xi - 2.0
        elif innertype == "proplyd":
            mlaw = proplyd_momentum
        elif innertype == "isotropic":
            mlaw = isotropic_momentum
        else:
            raise NotImplementedError("Inner wind must be isotropic or proplyd")

        if not outertype == "isotropic":
            raise NotImplementedError("Outer wind must be isotropic for now")

        w = Wind(axial_momentum_flux=beta, momentum_law=mlaw, origin=True)
        w1 = Wind(origin=False)
        # Initialise the base class
        BaseShell.__init__(self, w, w1)

        # Asymptotic angles for far wing/tail
        a, b = 1.0001*np.pi/2.0, np.pi
        if innertype == 'anisotropic':
            self.th_infty = scipy.optimize.brentq(_finf, a, b, args=(beta, xi))
        elif innertype == 'isotropic':
            self.th_infty = scipy.optimize.brentq(_finfCRW, a, b, args=(beta,))
        elif innertype == 'proplyd':
            self.th_infty = scipy.optimize.brentq(_finf, a, b, args=(beta, 0.8))

        self.th1_infty = np.pi - self.th_infty
            
###
### Private functions - these are implementation details that should
### not be accesible from outside
###

def _finf(th, beta, xi):
    """Function that gives f(theta) = 0 when theta = theta_infty
    
    Version for hemispheric flow with anisotropy xi
    """
    k = 2./xi-2
    C = (k+2*(1-beta))/(k+2)
    I = np.sqrt(np.pi)*gamma_func(0.5*(k+1))/(4*gamma_func(0.5*k+2))
    D = np.pi + 2*beta*I
    return th - C*np.tan(th) - D


def _finfCRW(th, beta):
    """Function that gives f(theta) = 0 when theta = theta_infty

    Version for spherically symmetric flow, as in CRW
    """
    return th - np.tan(th) - np.pi/(1.0 -beta)


def _radius_eq6(w, w1, th, th1):
    """Literal implementation of CRW96 Eq 6 for two winds w, w1

    Returns the radius for a given pair of angles th and th1 in terms
    of the momentum rates injected by the two winds

    """
    numerator = w.Jdot(th) + w1.Jdot(th1)
    denominator = (w.Pidot_r(th) + w1.Pidot_r(th1))*np.cos(th) \
                  - (w.Pidot_z(th) + w1.Pidot_z(th1))*np.sin(th)
    return numerator/denominator


def _radius_eq23(th, th1):
    """
    Literal implementation of CRW Eq 23 

    Gives the radius in terms of the two angles th and th1
    """
    return np.sin(th1)/np.sin(th+th1)

def _solve_for_th1(w, w1, th, th1_estimate=None, bounds=None, method='brent'):
    """For two winds (w and w1) and an angle (th) wrt the origin of w,
    find the angle th1 wrt the origin of w1

    It is necessary to give an initial estimate (th1_estimate) for
    th1. 

    Note that the internal function call is very expensive, since it
    potentially includes numerical integrations.  But who cares,
    right?!
    """

    def _f(th1, w, w1, th):
        """This should be zero when we have the correct th1"""
        return _radius_eq6(w, w1, th, th1) - _radius_eq23(th, th1)

    if method == 'fsolve':
        # This was the original method
        assert th1_estimate is not None, 'fsolve method needs guess for th1'
        th1, = scipy.optimize.fsolve(_f, th1_estimate, args=(w, w1, th))
    elif method == 'brent':
        assert bounds is not None and len(bounds) >= 2, 'Bounds must be 2-sequence'
        fa = _f(bounds[0], w, w1, th)
        fb = _f(bounds[1], w, w1, th)
        if DEBUG_LEVEL > 0:
            print(bounds, [fa, fb])
        if fa*fb < 0.0:
            # Hurray, we have bracketed the root
            a, b = bounds
        else:
            # Look for a change of sign in the middle somewhere
            xgrid = np.linspace(bounds[0], bounds[1], 2000)
            # Make sure it has the opposite sign from upper bound
            fgrid = -fb*np.array([_f(x, w, w1, th) for x in xgrid])
            if not np.any(fgrid > 0.0):
                # No roots at all - this is bad... bail out
                return np.nan
            i0 = np.argmax(fgrid)
            a, b = xgrid[i0], bounds[1]
            if DEBUG_LEVEL > 0:
                print(xgrid[i0], fgrid[i0])
        th1 = scipy.optimize.brentq(_f, a=a, b=b, args=(w, w1, th))
    else:
        raise NotImplementedError
    return th1


###
### What happens if we run the module as a script
###
if __name__ == "__main__":

    # Define a shell between two equal and isotropic winds
    shell = Shell(beta=1.0)

    # Define an array of angles
    theta = np.linspace(0.0, np.pi)

    # Calculate the shell radius for each angle
    R = shell.radius(theta)

    # Print the z coordinate of the shell: R cos(theta)
    print(R*np.cos(theta))       # These should all be 0.5

    
    # Now do the same for the proplyd case
    shell = Shell(beta=1.0, innertype="proplyd")
    R = shell.radius(theta)
    print(R*np.cos(theta))
    

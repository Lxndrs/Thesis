import sys
import json
import numpy as np
from bow_projection import Spline_R_theta_from_grid
sys.path.append("../CRW-shapes")
import equation6

DEBUG = False


class Ancantoid(object):
    def __init__(self, xi, beta, n=101):
        if DEBUG:
            print("Initialising Ancantoid(xi={xi:.2g}, beta={beta:.2g}, n={n})",
                  file=sys.stderr)
        self.xi = xi
        self.beta = beta
        self.n = n
        self.thgrid = np.linspace(0.0, np.pi, n)
        try:
            # First, look for cached file
            self._load_Rgrid_from_cache()
        except:
            # Failing that, use equation6 to find grid of R, theta
            self.shell = equation6.Shell(innertype='anisotropic', beta=beta, xi=xi)
            self.Rgrid = self.shell.radius(self.thgrid) / self.shell.R0
            self._save_Rgrid_to_cache()

        if DEBUG:
            print("thgrid =", self.thgrid, file=sys.stderr)
            print("Rgrid = ", self.Rgrid, file=sys.stderr)
        # Then set up the spline fit to the grid points
        self.splinefit = Spline_R_theta_from_grid(
            theta_grid=self.thgrid, R_grid=self.Rgrid)

    def __call__(self, theta):
        # When called as a function, give the spline fitted result
         return self.splinefit(theta)

    def _load_Rgrid_from_cache(self):
        with open(self._cache_filename()) as f:
            data = json.load(f)
        self.thgrid = np.array(data['theta'])
        self.Rgrid = np.array(data['R'])

    def _save_Rgrid_to_cache(self):
        data = {'theta': list(self.thgrid), 'R': list(self.Rgrid)}
        with open(self._cache_filename(), 'w') as f:
            json.dump(data, f, indent=4)

    def _cache_filename(self, suffix=".json"):
        fn = "ancantoid"
        fn += f"-xi{int(100*self.xi):03d}"
        fn += f"-beta{int(100000*self.beta):06d}"
        fn += f"-n{self.n:05d}"
        fn += suffix
        return fn


if __name__ == "__main__":

    from matplotlib import pyplot as plt
    import seaborn as sns

    lib_name = sys.argv[0].replace('.py', '')
    figfile = f"test_{lib_name}_radius.pdf"


    sns.set_style('ticks')
    fig, ax = plt.subplots()

    th = np.linspace(-np.pi, np.pi, 1001)
    th_dg = np.degrees(th)

    for xi, beta in [[0.8, 0.001],
                     [0.8, 0.01],
                     [0.8, 0.1],
                     [0.4, 0.001],
                     [0.4, 0.01],
                     [0.4, 0.1],]:
        label = fr"$\beta = {beta:.3f}$, $\xi = {xi:.1f}$"
        shape = Ancantoid(xi=xi, beta=beta)
        ax.plot(np.degrees(shape.thgrid), shape.Rgrid,
                color='b', alpha=0.2, lw=2, label='_nolabel_')
        ax.plot(th_dg, shape(th), lw=0.8, label=label)

    ax.legend(title=r"Ancantoid shapes")
    ax.set(
        xlabel=r"Polar angle: $\theta$, degrees",
        ylabel=r"$R$",
        xlim=[0, 180],
        yscale='log',
        ylim=[0.9, 200.0],
        xticks=[0, 30, 60, 90, 120, 150, 180],
    )
    sns.despine()
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end='')

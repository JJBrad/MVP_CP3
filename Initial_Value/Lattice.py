import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from PIL import Image


def laplacian(lat):
    """
    Returns the (discretised) laplacian at each point in a 2d lattice with periodic BCs (assuming dx == 1)
    """
    # np.roll moves all rows/cols (3rd arg) along by N spots (2nd arg), periodic conds.
    lap = np.roll(lat, 1, 0) + np.roll(lat, -1, 0) + np.roll(lat, 1, 1) + np.roll(lat, -1, 1) - 4. * lat
    return(lap)

def nablaLatSquared(lat):
    # TODO: Check this is right
    nab2 = np.square(np.roll(lat, 1, 0) - np.roll(lat, -1, 0)) + np.square(np.roll(lat, 1, 1) - np.roll(lat, -1, 1))
    return(nab2)

class lattice(object):
    """
    Lattice object for the SIRS model with built in dynamics and periodic boundary conditions. Each site has 
    one of 4 states; 0 is susceptible, 1 is infected, -1 is recovered, and 2 is immune.
    """
    def __init__(self, xDim=100, yDim=0, M=0.1, a=0.1, b=None, kappa=0.1, phi0=0., noise=0.01, dx=1., dt=1., measure=False, outDir="Data", label="Run", initCond="Uniform", status=True):
        """
        Constructor for the lattice object. Defaults to square lattice
        :param xDim: The x dimension of the lattice. Defaults to 50
        :param yDim: The y dimension of the lattice (optional), defaults to square lattice
        :param a, b, kappa: Parameters for the chemical potential, mu = -a phi + b phi^3 -kappa grad^2(phi)
        :param M: The motility of the suspension
        :param phi0: The initial value of the order parameter
        :param noise: The noise in the initial order paramter (gaussian standard deviation)
        :param dx: The step size in space
        :param dt: The step size in time
        :param measure: Whether to record measurements
        :param animate: Whether to animate the system
        :param outDir: 
        """
        self.M = M
        self.a = a
        if b is None:
            self.b = a
        else:
            self.b = b
        self.kappa = kappa
        self.dx = dx
        self.dt = dt
        
        self.outDir = outDir
        self.label = label
        self.path = self.outDir + "/" + self.label
        
        # Setup measurement 
        self.measure = measure
        
        if self.measure:
            self.nList = []
            self.FList = []        
            if os.path.exists(self.path):
                # Don't overwrite previous data
                print("Error. Output path {} already exists and measure is on. Exiting to avoid overwrite.".format(self.path))
                exit()
            elif os.path.exists(self.outDir): # Create run directory
                os.mkdir(self.path)
            else: # Create output and run directories
                os.mkdir(self.outDir)
                os.mkdir(self.path)
        
        # Initialise lattice
        self.n = 0 # Number of sweeps performed.
        self.xDim = xDim
        if yDim > 0:
            self.yDim = yDim
        else:
            self.yDim = xDim
        self.size = self.xDim*self.yDim
                
        if initCond == "binary":
            # Add gaussian noise to initial sample
            self.phi0 = phi0 + np.random.normal(0., noise)
            if self.phi0 > 1.: self.phi0 = 1.
            if self.phi0 < -1.: self.phi0 = -1.
            # Number of sites in state 1
            N = float(self.size)*(1. + self.phi0)/2.
            # Round to an integer, with probability accounting for rounding error.
            res = N - float(np.floor(N)) # Fraction of an extra site
            N = np.floor(N) # Round down
            if np.random.rand() > res: N = N + 1 # Add one depending on extra site
            
            # Define lattice
            sites = np.array([1]*N + [-1]*(self.size - N))
            np.random.shuffle(sites)
            self.phi = sites.reshape(self.xDim, self.yDim)
        else:
            self.phi0 = phi0
            self.phi = np.random.normal(self.phi0, noise, size=(self.xDim, self.yDim))
            #self.phi = np.random.uniform(-1.*noise, 1.*noise, size=(self.xDim, self.yDim)) + self.phi0
        
        if status:
            # Print state
            print(self)
            
    def freeEnergyDensity(self):
        f = -0.5*self.a*np.power(self.phi, 2.) + 0.25*self.b*np.power(self.phi, 4.) + 0.5*self.kappa*(1./(4.*self.dx**2.))*nablaLatSquared(self.phi)
        return(f)
            
    def __str__(self):
        """
        Returns string for printing key details of object.
        """
        return("Array has shape {}. Average phi is {:.3e}.".format(self.phi.shape, np.average(self.phi)))
    
    def chemPot(self):
        """
        Returns chemical potential at each discrete point in the lattice.
        """
        # np.roll moves all rows/cols (3rd arg) along by N spots (2nd arg), periodic conds.
        mu = -1.*self.a*self.phi + self.b*np.power(self.phi, 3.) - (self.kappa/(self.dx**2.))*laplacian(self.phi)
        return(mu)
    
    def next(self):
        """
        Perform one sweep.
        """
        self.phi = self.phi + (self.M*self.dt/(self.dx**2.)) * laplacian(self.chemPot())
        self.n += 1
        if self.measure:
            self.nList.append(self.n)
            self.FList.append(np.sum(self.freeEnergyDensity()))
            
    def analyse(self, showPlot = True):
        with open(self.path + "/Data.csv", "w") as outFile:
            outFile.write("n, F\n")
            for r in range(0, len(self.nList)):
                outFile.write("{},{}".format(self.nList[r], self.FList[r]))
        pyplot.plot(self.nList, self.FList, "k-")
        pyplot.xlabel("Steps (n)")
        pyplot.ylabel("Free Energy")
        pyplot.savefig(self.path + "/FreeEnergy.png", dpi=150)
        if showPlot:
            pyplot.show()
        else:
            pyplot.close()
        return()

    def animate(self, f, tMax, rate):
        for _ in range(0, rate):
            if self.n <= tMax:
                self.next()
        if f % 100 == 0:
            print("{} steps completed".format(f * rate))
            #print("Max: {}, min: {}".format(np.max(self.phi), np.min(self.phi)))
        if self.n == tMax:
            print("Animation complete.")
            self.n += rate
        im = axis.imshow(self.phi, cmap="bwr", interpolation="nearest", vmin=-1, vmax=1) #"seismic"
        pyplot.clim(-1, 1)
        return([im])
        
    def display(self, tMax=20000, rate=10):
        global axis, fig, cax
        fig, axis = pyplot.subplots()
        for s in axis.spines: axis.spines[s].set_color("k")
        axis.spines['left'].set_position(('outward', 1))
        axis.spines['top'].set_position(('outward', 0.5))
        # Position colourbar
        im = pyplot.imshow(self.phi, cmap="bwr", interpolation="nearest", vmin=-1, vmax=1) #"bwr"
        pyplot.clim(-1, 1)
        cbar = pyplot.colorbar(im)
        anim = FuncAnimation(fig, self.animate, tMax/rate + 1, interval=0, blit=True, repeat=False, fargs=(tMax,rate))
        pyplot.show()
        print("Animation finished. Please close the animation window.")

    def run(self, tMax=1000):
        for i in range(0, tMax):
            self.next()
        

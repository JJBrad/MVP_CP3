import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from PIL import Image
from scipy.optimize import curve_fit

def plotLat(lat, xCross, yCross, z, interp, norm, fileName):
        im = pyplot.imshow(lat, cmap="plasma", interpolation=interp, norm=norm) #"bwr"
        pyplot.axhline(y=yCross, color="k", lw=0.5, ls="--")
        pyplot.axvline(x=xCross, color="k", lw=0.5, ls="--")
        pyplot.title("i = {}".format(z))
        pyplot.xlabel("j")
        pyplot.ylabel("k")
        pyplot.colorbar()
        pyplot.savefig(fileName, dpi=180)
        pyplot.close()

# Define function for fitting (linear in this case)
def fitFunc(x, m, c):
    return(m*x + c)

class lattice(object):
    """
    Lattice object for the SIRS model with built in dynamics and periodic boundary conditions. Each site has 
    one of 4 states; 0 is susceptible, 1 is infected, -1 is recovered, and 2 is immune.
    """
    def __init__(self, xDim=50, yDim=0, zDim=0, dx=1., dt=-1., q=1., tolerance=1E-3, omega=1., outDir="Data", label="Run", initCond="Uniform", method="Jacobi", status=True):
        """
        Constructor for the lattice object. Defaults to square lattice
        :param xDim: The x dimension of the lattice. Defaults to 50
        :param yDim: The y dimension of the lattice (optional), defaults to square lattice
        :param a, b, kappa: Parameters for the chemical potential, mu = -a phi + b phi^3 -kappa grad^2(phi)
        :param M: The motility of the suspension
        :param phi0: The initial value of the order parameter
        :param noise: The noise in the initial order paramter (gaussian standard deviation)
        :param omega: Over-relaxation parameter
        :param dx: The step size in space
        :param dt: The step size in time
        :param measure: Whether to record measurements
        :param animate: Whether to animate the system
        :param outDir: 
        """
        self.dx = dx
        if dt > 0:
            self.dt = dt
        else:
            self.dt = 1.
        
        self.omega = omega
        self.tolerance = tolerance
        if method not in ["Jacobi", "GS"]:
            print("Method {} not recognised. Using Jacobi.")
            method = "Jacobi"
        self.method = method
        self.outDir = outDir
        self.label = label
        self.path = self.outDir + "/" + self.label

        if os.path.exists(self.path):
            # Don't overwrite previous data
            print("Error. Output path {} already exists and measure is on. Exiting to avoid overwrite.".format(self.path))
            # exit() #TODO
        elif os.path.exists(self.outDir): # Create run directory
            os.mkdir(self.path)
        else: # Create output and run directories
            os.mkdir(self.outDir)
            os.mkdir(self.path)

        # Initialise lattice
        self.n = 0 # Number of sweeps performed.
        self.xDim = xDim
        if yDim > 2:   #Outer layer fixed at zero, so dimensions must be larger than 2 in each direction
            self.yDim = yDim
        else:
            self.yDim = xDim
        if zDim > 2:
            self.zDim = zDim
        else:
            self.zDim = xDim
        self.size = self.xDim*self.yDim*self.zDim
        
        # Prepare lattice        
        self.phi = np.zeros((self.xDim, self.yDim, self.zDim))
        self.rho = np.zeros((self.xDim, self.yDim, self.zDim))
        self.approxMid = (self.xDim/2, self.yDim/2, self.zDim/2)
        self.rho[self.approxMid] = q
        
        if status:
            # Print state
            print(self)

    def __str__(self):
        """
        Returns string for printing key details of object.
        """
        return("Array has shape {}. Total charge is {:.3e}.".format(self.rho.shape, np.sum(self.rho)))
    
    def next(self):
        if self.method == "Jacobi":
            return(self.jacobiUpdate())
        elif self.method == "GS":
            return(self.GSUpdate())
        
    def jacobiUpdate(self):
        """
        Perform one sweep using the Jacobi algorithm. Return total absolute change in phi.
        """
        # Calculate terms in new new lattice
        update = np.roll(self.phi, 1, 0) + np.roll(self.phi, -1, 0)\
               + np.roll(self.phi, 1, 1) + np.roll(self.phi, -1, 1)\
               + np.roll(self.phi, 1, 2) + np.roll(self.phi, -1, 2)\
               + self.rho
        # Create new lattice of zeros
        newPhi = np.zeros((self.xDim, self.yDim, self.zDim))
        # Update all other than boundaries
        newPhi[1:-1, 1:-1, 1:-1] = update[1:-1, 1:-1, 1:-1]/6.
        
        # Calculate difference in potentials
        dPhi = np.sum(np.abs(newPhi - self.phi))
        # Update lattice
        self.phi = newPhi
        self.n += 1
        return(dPhi)
        
    def GSUpdate(self):
        """
        Perform one sweep using the Gauss-Siedel method. Return total absolute change in phi.
        """
        # Keep copy of old lattice so a threshold can be defined.
        oldPhi = np.copy(self.phi)
        
        # Update all entries other than boundaries
        for i in range(1, self.xDim-1):
            for j in range(1, self.yDim-1):
                for k in range(1, self.zDim-1):
                    self.phi[i, j, k] = (1./6.) * (self.phi[i+1, j, k] + self.phi[i-1, j, k] + \
                                                   self.phi[i, j+1, k] + self.phi[i, j-1, k] + \
                                                   self.phi[i, j, k+1] + self.phi[i, j, k-1] + \
                                                   self.rho[i, j, k])
                    # Apply over-relaxation:
                    self.phi[i,j,k] = self.omega*self.phi[i,j,k] + (1. - self.omega)*oldPhi[i,j,k]
        
        # Calculate total absolute difference in potentials
        dPhi = np.sum(np.abs(self.phi - oldPhi))
        # Update time
        self.n += 1
        return(dPhi)
            
    def analyse(self, showPlot = True):
        """
        Function to perform analysis of result for Poisson
        """
        # Show 2d heatmap of phi
        self.showCut()
        # Obtain lists of phi as a function of radius
        rVals = []
        phiVals = []
        # Open file to write results
        outFile = open(self.path + "/vals.csv", "w")
        outFile.write("i,j,k,r,phi\n")
        # Loop through array
        for i in range(0, self.xDim):
            for j in range(0, self.yDim):
                for k in range(0, self.zDim):
                    # Calculate distance
                    r = self.dx* np.sqrt(float(i - self.approxMid[0])**2. +\
                                         float(j - self.approxMid[1])**2. +\
                                         float(k - self.approxMid[2])**2.)
                    # Keep results only for r > 0 (to avoid undefined logarithms)
                    if r > 0.:
                        phiVals.append(self.phi[i,j,k])
                        rVals.append(r)
                    # Write to file
                    outFile.write("{},{},{},{},{}\n".format(i, j, k, r,self.phi[i,j,k]))
        outFile.close()
        # Convert to logarithms
        rVals = np.log10(rVals)
        phiVals = np.log10(phiVals)
        # Fit line in logarithmic space
        # Initial guess for fitting
        p0 = [-1, 0.1]
        
        # Perform fit up to log(R) = 0.75
        params, var = curve_fit(fitFunc, rVals[rVals < 0.75], phiVals[rVals < 0.75], p0=p0)
        err = np.sqrt(np.diag(var))
        
        xFit = np.linspace(min(rVals), max(rVals), 500)
        yFit = fitFunc(xFit, *params)
        pyplot.plot(xFit, yFit, "r--")
        fitStr = "Fit: $\log(\phi) = ({:.4f} \pm {:.4f})\log(R) + ({:.5f} \pm {:.5f})$".format(params[0], err[0], params[1], err[1])
        pyplot.figtext(0.89, 0.84, fitStr, ha="right")
        # Plot phi as a function of radius
        pyplot.plot(rVals, phiVals, "b.", ms=0.5)
        pyplot.savefig(self.path + "\RvsPhi.png", dpi=150)
        pyplot.xlabel(r"Log$_{10}$ R")
        pyplot.ylabel(r"Log$_{10}$ $\phi$")
        if showPlot: pyplot.show()
        else: pyplot.close()
        
        # Select only xy plane with charge
        phiSlice = self.phi[:,:,self.approxMid[2]]
        # Electric field is -del(phi). Get x and y components (assume dx = dy = 1)
        yField = -(np.roll(phiSlice, -1, 0) - np.roll(phiSlice, 1, 0))/(2. * self.dx**2.) # numpy arrays indexed like matrices (Row, column)
        xField = -(np.roll(phiSlice, -1, 1) - np.roll(phiSlice, 1, 1))/(2. * self.dx**2.)
        norm = np.sqrt(np.square(xField) + np.square(yField))[1:-1,1:-1]
        xField = xField[1:-1,1:-1]/norm
        yField = yField[1:-1,1:-1]/norm
        xCoords = range(1, self.xDim-1)
        yCoords = range(1, self.yDim-1)
        quiv = pyplot.quiver(xCoords, yCoords, xField, yField)
        pyplot.savefig(self.path + "\E_Field.png")
        if showPlot: pyplot.show()
        else: pyplot.close()
        
    def showCut(self):
        plotSlice = self.phi[self.approxMid[0], :, :]
        plotLat(plotSlice, self.approxMid[2], self.approxMid[1], self.approxMid[0], "nearest", None, self.path+"/Cut_1.png")
        plotLat(plotSlice, self.approxMid[2], self.approxMid[1], self.approxMid[0], "bilinear", None, self.path+"/Cut_2.png")
        plotLat(plotSlice, self.approxMid[2], self.approxMid[1], self.approxMid[0], "nearest", colors.LogNorm(), self.path+"/Cut_3.png")
        plotLat(plotSlice, self.approxMid[2], self.approxMid[1], self.approxMid[0], "bilinear", colors.LogNorm(), self.path+"/Cut_4.png")

    def run(self, tMax=100000, showPlot=True):
        for i in range(0, tMax):
            dPhi = self.next()
            if self.n % 200 == 0:
                print("Step {} complete. dPhi is {:.3e}".format(self.n, dPhi))
            if dPhi < self.tolerance:
                self.complete = True
                print("Converged after step {}. dPhi is {:.3e}".format(self.n, dPhi))
                break
        if not self.complete:
            print("Maximum iterations ({}) reached. Last dPhi was {:.3e} and tolerance not matched. Results may not be valid.".format(tMax, dPhi))
        self.analyse(showPlot=showPlot)
        return(self.n)
        
"""
TODO:
In magnetic field plot, consider only Az, so this replaces phi. Also, the field is no longer E = -del phi, it becomes B = del x A.
Magnetic field thing needs only 2 dimensional. 
Also do overcompensating correction thingy. Steps of 0.01 between 1 and 2
Plot E field as a function of R too.
"""

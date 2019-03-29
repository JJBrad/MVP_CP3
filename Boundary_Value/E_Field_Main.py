import numpy as np
import Lattice as lat
import interactive as interact
import matplotlib.pyplot as pyplot
import sys

# Default values:
params = {"X Dimension":50,
          "Y Dimension":-1,
          "phi0" : -0.5,
          "a" : 0.1,
          "b" : None,
          "kappa" : 0.1,
          "M" : 0.1,
          "noise" : 0.01,
          "tolerance" : 1.0E-3,
          "dx" : 1.,
          "dt" : 1.,
          "Dynamics" : "Jacobi",
          "Seed" : None,
          "tMax" : 100000,
          "omega" : 1.,
          "Rate" : 100, # Number of sweeps per update frame
          "RunLabel" : "Run",
          "outDir" : "Data"
          }

# Get input from command line
args = sys.argv[1:]
interact.readArgs(args, params)
np.random.seed(params["Seed"]) #None is default, changes each run.

lattice = lat.lattice(label=params["RunLabel"],
                      method=params["Dynamics"],
                      tolerance=params["tolerance"],
                      xDim=params["X Dimension"],
                      omega=params["omega"],
                      outDir=params["outDir"])
                      
lattice.run()

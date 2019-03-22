import numpy as np
import Lattice as lat
import interactive as interact
import matplotlib.pyplot as pyplot
import sys
from scipy.optimize import curve_fit

# Default values:
params = {"X Dimension":50,
          "Y Dimension":-1,
          "phi0" : -0.5,
          "a" : 0.1,
          "b" : None,
          "kappa" : 0.1,
          "M" : 0.1,
          "noise" : 0.01,
          "dx" : 1.,
          "dt" : 1.,
          "Seed" : None,
          "tMax" : 100000,
          "Rate" : 100, # Number of sweeps per update frame
          "Measure" : True,
          "Animate" : True,
          "RunLabel" : "Run",
          "outDir" : "Data"
          }

# Get input from command line
args = sys.argv[1:]
interact.readArgs(args, params)
np.random.seed(params["Seed"]) #None is default, changes each run.

if not (params["Measure"] or params["Animate"]):
    print("Error, system is set to neither animate nor measure. Exiting...")
    exit()

lattice = lat.lattice(params["X Dimension"], 
                      params["Y Dimension"], 
                      params["M"], 
                      params["a"], 
                      params["b"], 
                      params["kappa"], 
                      params["phi0"], 
                      params["noise"],
                      params["dx"],
                      params["dt"],
                      measure=params["Measure"],
                      label=params["RunLabel"],
                      outDir=params["outDir"])

if params["Animate"]:
    lattice.display(tMax=params["tMax"], rate=params["Rate"])
else:
    lattice.run(tMax=params["tMax"])
    
if params["Measure"]:
    lattice.analyse()

print("#"*40 + "\nNote: If animation was exited manually then an error may appear above.\nDisregard this error.\n" + "#"*40)

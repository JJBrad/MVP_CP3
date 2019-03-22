###### MVP Checkpoint 2 Part 1 ######
Python Cahn-Hilliard simulation.
Chemical potential is mu = -a*phi + b*phi - k * del**2(phi)
Defaults:
Lattice dimensions: 50 x 50
Initial state: Random
Number of Sweeps: 10000
Animation: On
Random seed: varies
Equilibration time: 150 sweeps
Autocorrelation time: 20 sweeps


Usage: python Ising.py <Tags>
Optional tags:
-x <value>        Lattice x dimension.
-y <value>        Lattice y dimension, defaults to match x.
-a <value>        Value of a for chemical potential
-b <value>        Value of b for chemical potential, defaults to match a.
-k <value>        Value of k for chemical potential
-M <value>        Motility
-dx <value>       Step size (space)
-dt <value>       Step size (time)
-i <values>       Initial value of the order parameter.
-rs <value>       Set random seed.
-N <values>       Number of sweeps to perform.
-D <Y/N>          Display animation? Yes (Y) or No (N)
-P <Y/N>          Turn plotting (and measurements) on or off
-u <value>        Update rate of the animation (sweeps between frames)
-r <name>         Run name
-o <dir>          Output directory name
-H                Print this dialogue and exit.

Example with 40 x 25 lattice with 1000 sweeps, and p1=p2=p3=0.5:
python Main.py -x 40 -y 25 -N 1000 -p [0.5,0.5,0.5]

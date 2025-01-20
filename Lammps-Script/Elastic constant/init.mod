# NOTE: This script was modified for calculate elastic constants at finite temperatures, 
# units, etc. See in.elastic for more info.
#

#------------------------ Random number generator -----------------------------#
  # Initalizes random number generator
  variable          rnd equal round(random(0,999999,1234532))
#------------------------------------------------------------------------------#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
include lmp.in
replicate 6 6 3
variable up equal 1.0e-2
 
# metal units, elastic constants in GPa
#units		metal
variable cfac equal 1.01325e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal 0.0 
variable ftol equal 1.0e-10
variable maxiter equal 1000
variable maxeval equal 1000
variable dmax equal 1.0e-2
variable N equal 20000

# Define dynamic parameters
variable dt equal 1
variable T equal 300.0
variable press equal 0.0
variable tnpt equal 16020
variable tequi equal 10000.0
variable tequihalf equal ${tequi}/2
#variable tdumpeq equal ${tequi}/10

variable erate equal 1e-4   ## this gives strain rate: 1e-4/1e-12 = 1e8/s.
variable trun equal 2*${up}/(${erate}*${dt})
#variable tdump equal ${trun}/10

# set up structure
#read_data       ../models/${T}K.dat
#read_restart ../equi_${T}K/restart.npt_${T}K
read_dump dump.npt_${T}K ${N} x y z vx vy vz



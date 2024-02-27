import numpy as np
from matplotlib import pyplot as plt
from CheapTools import *

# Represents the [O/Fe] and [Si/Fe] for the DTDs used in our work
TypeIa_SNe_ratio = 0.54/100.*1E9;# +/-0.12 events/cent
Area = np.pi*(20.**2-3.**2)*1E6;# In pc**2
today = 13.8; # Gyr
Solar_values = {"FeH":-2.752, "OFe":0.646,"SiFe":-0.291}# Fe/H, O/Fe, Si/Fe. Related to the ratio of solar Fe, O, and Si densities.

# Integration time
t_gyr = np.arange(0.01, today, 0.00125)# Gyr
t_gyr = t_gyr[t_gyr<today]

# Model parameters
chemdict = dict()
chemdict["omega"] = 0.4
chemdict["R"] = 0.285
chemdict["nuL"] = 2.

# Infall parameters:
chemdict["tauj"] = np.array([7.])
chemdict["tj"] = [0.]
chemdict["sigma_gas_0"] = 1E-8 # Sigma_gas_0 should be very close to zero but not zero
chemdict["Aj"] =  [9.98032680842189] # Makes 54 Msun/pc**2 today (Vincenzo et al. 2017)

# Now add the parameters associated with the iron element, with zero initial density.
#    We can it both by:
chemdict["yx"] = 5.6E-4# Use your preferred yield here
chemdict["mx1a"] = 6.26E-01# Use your preferred yield here
chemdict["sigmaX_0"] = 0.0# Let's assume no initial content of iron
#    ... or using the add_element() function for the values considered in Palicio et al. (submitted).
chemdict = add_element(chemdict, "Fe", 0.0)# For the moment, this function only works for Fe, O and Si,

# DTD parameters:
chemdict = Load_MR01_dict( chemdict )# For example, let's use the MR01 DTD
# Using "chemdict" as input, the output will have all the key-values of the input

# Now we have to provide a value for CIa, but instead of setting CIa directly we make
# use of the present-day Type Ia ratio:
chemdict["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, chemdict, present_day_time=today)# This computes CIa

# Solve the Chemical Evolution Model equation:
Sigma_MR01_Fe = SolveChemEvolModel( t_gyr, chemdict)# The output is the density of iron as a function of time

# It is more intuitive to work with [Fe/H] rather than sigma_Fe:
FeH_MR01 = FromSigmaToAbundance(t_gyr, Sigma_MR01_Fe, chemdict)# From sigma (surface density) to abundance:

# We have to implement the correction due to the solar values for the iron:
FeH_MR01 = FeH_MR01 -Solar_values["FeH"] + 0.125 # The factor 0.125 comes from log10(0.75), since 3/4 of the gas is made by Hydrogen
#-------------------------------------------------------------------



# Plotting the results
# -------------------------------------------------------------
# Figure
fig = plt.figure("Iron")

# In the upper panel, plot the density of iron
plt.subplot(211)
plt.plot(t_gyr, Sigma_MR01_Fe, label=r"$\sigma_{Fe}(\rm t)$", color="b")
plt.xlabel( "Time (Gyr)", fontsize=14)
plt.ylabel( r"$\sigma_{Fe} (M_{\odot}/pc^2)$", fontsize=14)

# In the lower panel, plot [Fe/H]
plt.subplot(212)
plt.plot(t_gyr, FeH_MR01, label="[Fe/H]", color="b")
plt.xlabel( "Time (Gyr)", fontsize=14)
plt.ylabel( "[Fe/H] (dex)", fontsize=14)

plt.show()
fig.savefig("QuickTestFigure.png")
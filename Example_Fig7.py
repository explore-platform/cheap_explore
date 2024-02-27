import numpy as np
from matplotlib import pyplot as plt
from CheapTools import *

# Represents the [O/Fe] and [Si/Fe] for the DTDs used in our work
TypeIa_SNe_ratio = 0.54/100.*1E9;# +/-0.12 events/cent
Area = np.pi*(20.**2-3.**2)*1E6;# In pc**2
today = 13.8; # Gyr
Ninfall = 1;# Number of infalls

fontsize = 26

Solar_values = {"FeH":-2.752, "OFe":0.646,"SiFe":-0.291}# Fe/H, O/Fe, Si/Fe

# Integration time
t_gyr = np.arange(0.01, today, 0.00125)# Gyr
t_gyr = t_gyr[t_gyr<today]

Model_and_infall_parameters = {"omega" : 0.4, "R" : 0.285, "nuL" : 2., "tauj" : 1*np.array(Ninfall*[7.]), "tj" : Ninfall*[0.], "sigma_gas_0" : 1E-8}# Sigma_gas_0 should be very close to zero but not zero
Model_and_infall_parameters["Aj"] =  [9.98032680842189] # Makes 54 Msun/pc**2 today (Vincenzo et al. 2017)

# ----------------------------------------------------------------------------------------------------




# Fit of the MR01 DTD
#------------------------------------------------------------
# Load the fit of the MR01 DTD
dict_MR01 = Load_MR01_dict();
for key, value in Model_and_infall_parameters.items(): dict_MR01[key] = value; # Add the parameters of the model and the infall

# Compute CIa
dict_MR01["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_MR01, present_day_time=today)

# Now add the values of the elements, with zero initial density:
dict_MR01_Fe = add_element(dict_MR01, "Fe", 0.0);
dict_MR01_O = add_element(dict_MR01, "O", 0.0);
dict_MR01_Si = add_element(dict_MR01, "Si", 0.0);
del(dict_MR01)

# Solve the models for all the elements simultaneously:
Sol_MR01_Fe = SolveChemEvolModel( t_gyr, dict_MR01_Fe)
Sol_MR01_O = SolveChemEvolModel( t_gyr, dict_MR01_O)
Sol_MR01_Si = SolveChemEvolModel( t_gyr, dict_MR01_Si)

# From sigma (surface density) to abundance:
Abund_MR01_Fe = FromSigmaToAbundance(t_gyr, Sol_MR01_Fe, dict_MR01_Fe)
Abund_MR01_O = FromSigmaToAbundance(t_gyr, Sol_MR01_O, dict_MR01_O)
Abund_MR01_Si = FromSigmaToAbundance(t_gyr, Sol_MR01_Si, dict_MR01_Si)
#-------------------------------------------------------------------



# Mannucci+06 DTD
#------------------------------------------------------------
dict_M06 = { "AG" : 19.95, "taup" : 0.05, "sigma_p" : 0.01, "AE" : 0.17, "tauD" : 3.,  "tau1" : 0.03, "tau2" : 10.05}
for key, value in Model_and_infall_parameters.items(): dict_M06[key] = value; # Add the parameters of the model and the infall
dict_M06 = prepare_chemdict(dict_M06) # Check the format
# Compute CIa
dict_M06["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_M06, present_day_time=today)
# Now add the values of the elements, with zero initial density:
dict_M06_Fe = add_element(dict_M06, "Fe", 0.0);
dict_M06_O = add_element(dict_M06, "O", 0.0);
dict_M06_Si = add_element(dict_M06, "Si", 0.0);
del(dict_M06)

# Solve the models for all the elements simultaneously:
Sol_M06_Fe = SolveChemEvolModel( t_gyr, dict_M06_Fe)
Sol_M06_O = SolveChemEvolModel( t_gyr, dict_M06_O)
Sol_M06_Si = SolveChemEvolModel( t_gyr, dict_M06_Si)

# From sigma (surface density) to abundance:
Abund_M06_Fe = FromSigmaToAbundance(t_gyr, Sol_M06_Fe, dict_M06_Fe)
Abund_M06_O = FromSigmaToAbundance(t_gyr, Sol_M06_O, dict_M06_O)
Abund_M06_Si = FromSigmaToAbundance(t_gyr, Sol_M06_Si, dict_M06_Si)
#-------------------------------------------------------------------




# Strolger+05 DTD
#------------------------------------------------------------
dict_S05 = { "AG" : 1., "taup" : 3.4, "sigma_p" : 0.68, "tau1" : 0.25, "tau2" : today+0*4.4}# tau1 is in Matteucci 09, tau2 in page 217 in Strolger+05
for key, value in Model_and_infall_parameters.items(): dict_S05[key] = value; # Add the parameters of the model and the infall
dict_S05 = prepare_chemdict(dict_S05) # Check the format
# Compute CIa
dict_S05["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_S05, present_day_time=today)
# Now add the values of the elements, with zero initial density:
dict_S05_Fe = add_element(dict_S05, "Fe", 0.0);
dict_S05_O = add_element(dict_S05, "O", 0.0);
dict_S05_Si = add_element(dict_S05, "Si", 0.0);
del(dict_S05)

# Solve the models for all the elements simultaneously:
Sol_S05_Fe = SolveChemEvolModel( t_gyr, dict_S05_Fe)
Sol_S05_O = SolveChemEvolModel( t_gyr, dict_S05_O)
Sol_S05_Si = SolveChemEvolModel( t_gyr, dict_S05_Si)

# From sigma (surface density) to abundance:
Abund_S05_Fe = FromSigmaToAbundance(t_gyr, Sol_S05_Fe, dict_S05_Fe)
Abund_S05_O = FromSigmaToAbundance(t_gyr, Sol_S05_O, dict_S05_O)
Abund_S05_Si = FromSigmaToAbundance(t_gyr, Sol_S05_Si, dict_S05_Si)
#-------------------------------------------------------------------




# Totani+08 DTD
#------------------------------------------------------------
dict_T08 = { "AI" : 1., "tauI" : 1., "tau0": 0., "tau1" : 0.1, "tau2" : 10.}
for key, value in Model_and_infall_parameters.items(): dict_T08[key] = value; # Add the parameters of the model and the infall
dict_T08 = prepare_chemdict(dict_T08) # Check the format
# Compute CIa
dict_T08["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_T08, present_day_time=today)
# Now add the values of the elements, with zero initial density:
dict_T08_Fe = add_element(dict_T08, "Fe", 0.0);
dict_T08_O = add_element(dict_T08, "O", 0.0);
dict_T08_Si = add_element(dict_T08, "Si", 0.0);
del(dict_T08)

# Solve the models for all the elements simultaneously:
Sol_T08_Fe = SolveChemEvolModel( t_gyr, dict_T08_Fe)
Sol_T08_O = SolveChemEvolModel( t_gyr, dict_T08_O)
Sol_T08_Si = SolveChemEvolModel( t_gyr, dict_T08_Si)

# From sigma (surface density) to abundance:
Abund_T08_Fe = FromSigmaToAbundance(t_gyr, Sol_T08_Fe, dict_T08_Fe)
Abund_T08_O = FromSigmaToAbundance(t_gyr, Sol_T08_O, dict_T08_O)
Abund_T08_Si = FromSigmaToAbundance(t_gyr, Sol_T08_Si, dict_T08_Si)
#------------------------------------------------------------


# Fit of the Wide G05 DTD
#------------------------------------------------------------
# Load the fit of the Wide G05 DTD
dict_G05Wide = Load_G05Wide_dict();
for key, value in Model_and_infall_parameters.items(): dict_G05Wide[key] = value; # Add the parameters of the model and the infall

# Compute CIa
dict_G05Wide["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_G05Wide, present_day_time=today)

# Now add the values of the elements, with zero initial density:
dict_G05Wide_Fe = add_element(dict_G05Wide, "Fe", 0.0);
dict_G05Wide_O = add_element(dict_G05Wide, "O", 0.0);
dict_G05Wide_Si = add_element(dict_G05Wide, "Si", 0.0);
del(dict_G05Wide)

# Solve the models for all the elements simultaneously:
Sol_G05Wide_Fe = SolveChemEvolModel( t_gyr, dict_G05Wide_Fe)
Sol_G05Wide_O = SolveChemEvolModel( t_gyr, dict_G05Wide_O)
Sol_G05Wide_Si = SolveChemEvolModel( t_gyr, dict_G05Wide_Si)

# From sigma (surface density) to abundance:
Abund_G05Wide_Fe = FromSigmaToAbundance(t_gyr, Sol_G05Wide_Fe, dict_G05Wide_Fe)
Abund_G05Wide_O = FromSigmaToAbundance(t_gyr, Sol_G05Wide_O, dict_G05Wide_O)
Abund_G05Wide_Si = FromSigmaToAbundance(t_gyr, Sol_G05Wide_Si, dict_G05Wide_Si)
#-------------------------------------------------------------------



# Fit of the Close G05 DTD
#------------------------------------------------------------
# Load the fit of the Close G05 DTD
dict_G05Close = Load_G05Close_dict();
for key, value in Model_and_infall_parameters.items(): dict_G05Close[key] = value; # Add the parameters of the model and the infall

# Compute CIa
dict_G05Close["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_G05Close, present_day_time=today)

# Now add the values of the elements, with zero initial density:
dict_G05Close_Fe = add_element(dict_G05Close, "Fe", 0.0);
dict_G05Close_O = add_element(dict_G05Close, "O", 0.0);
dict_G05Close_Si = add_element(dict_G05Close, "Si", 0.0);
del(dict_G05Close)

# Solve the models for all the elements simultaneously:
Sol_G05Close_Fe = SolveChemEvolModel( t_gyr, dict_G05Close_Fe)
Sol_G05Close_O = SolveChemEvolModel( t_gyr, dict_G05Close_O)
Sol_G05Close_Si = SolveChemEvolModel( t_gyr, dict_G05Close_Si)

# From sigma (surface density) to abundance:
Abund_G05Close_Fe = FromSigmaToAbundance(t_gyr, Sol_G05Close_Fe, dict_G05Close_Fe)
Abund_G05Close_O = FromSigmaToAbundance(t_gyr, Sol_G05Close_O, dict_G05Close_O)
Abund_G05Close_Si = FromSigmaToAbundance(t_gyr, Sol_G05Close_Si, dict_G05Close_Si)
#-------------------------------------------------------------------



# Fit of the P08 DTD
#------------------------------------------------------------
# Load the fit of the P08 DTD
dict_P08 = Load_P08_dict();
for key, value in Model_and_infall_parameters.items(): dict_P08[key] = value; # Add the parameters of the model and the infall

# Compute CIa
dict_P08["CIa"] = Get_CIa(TypeIa_SNe_ratio, Area, dict_P08, present_day_time=today)

# Now add the values of the elements, with zero initial density:
dict_P08_Fe = add_element(dict_P08, "Fe", 0.0);
dict_P08_O = add_element(dict_P08, "O", 0.0);
dict_P08_Si = add_element(dict_P08, "Si", 0.0);
del(dict_P08)

# Solve the models for all the elements simultaneously:
Sol_P08_Fe = SolveChemEvolModel( t_gyr, dict_P08_Fe)
Sol_P08_O = SolveChemEvolModel( t_gyr, dict_P08_O)
Sol_P08_Si = SolveChemEvolModel( t_gyr, dict_P08_Si)

# From sigma (surface density) to abundance:
Abund_P08_Fe = FromSigmaToAbundance(t_gyr, Sol_P08_Fe, dict_P08_Fe)
Abund_P08_O = FromSigmaToAbundance(t_gyr, Sol_P08_O, dict_P08_O)
Abund_P08_Si = FromSigmaToAbundance(t_gyr, Sol_P08_Si, dict_P08_Si)
#-------------------------------------------------------------------



# Solve for the IRA model:
Model_and_infall_parameters["CIa"] = 0.0; # Because it is the IRA
Model_and_infall_parameters = prepare_chemdict(Model_and_infall_parameters) # Check the format
Model_and_infall_parameters = add_element(Model_and_infall_parameters, "Fe", 0.0);
Sol_IRA_FE = SolveChemEvolModel( t_gyr, Model_and_infall_parameters)


# Plotting ALL the DTDs
# -------------------------------------------------------------
# Figure
fig, ax = plt.subplots(1, 2, num="Oxygen and Silicon", figsize=(13.5, 6.75), sharey=True)
ax = ax.flatten()

# Plot [O/Fe]
ax[0].plot(Abund_MR01_Fe-Solar_values["FeH"]+0.125, Abund_MR01_O-Abund_MR01_Fe-Solar_values["OFe"], color="blue", linewidth=2.5, label="MR01", zorder=1)# MR+01 (0.125 approx log10(0.75)
ax[0].plot(Abund_M06_Fe-Solar_values["FeH"]+0.125, Abund_M06_O-Abund_M06_Fe-Solar_values["OFe"], color="forestgreen", linewidth=2.5, label="MVP06", zorder=3)# M+06 (Original)
ax[0].plot(Abund_S05_Fe-Solar_values["FeH"]+0.125, Abund_S05_O-Abund_S05_Fe-Solar_values["OFe"], color="orange", linewidth=2.5, label="S05", zorder=4)# S+05 (Original)
ax[0].plot(Abund_T08_Fe-Solar_values["FeH"]+0.125, Abund_T08_O-Abund_T08_Fe-Solar_values["OFe"], color="deepskyblue", linewidth=2.5, label="T08", zorder=2)# T+08 (Original)
ax[0].plot(Abund_G05Wide_Fe-Solar_values["FeH"]+0.125, Abund_G05Wide_O-Abund_G05Wide_Fe-Solar_values["OFe"], color="mediumorchid", linewidth=2.5, label="WIDE G05", zorder=10)# Wide G+05 (0.125 approx log10(0.75)
ax[0].plot(Abund_G05Close_Fe-Solar_values["FeH"]+0.125, Abund_G05Close_O-Abund_G05Close_Fe-Solar_values["OFe"], color="pink", linewidth=2.5, label="CLOSE G05", zorder=10)# Close G+05 (0.125 approx log10(0.75)
ax[0].plot(Abund_P08_Fe-Solar_values["FeH"]+0.125, Abund_P08_O-Abund_P08_Fe-Solar_values["OFe"], color="brown", linewidth=2.5, label="P08", zorder=10)# P+08 (0.125 approx log10(0.75)
# Axis format
for label in ['top','bottom','left','right']: ax[0].spines[label].set_linewidth(2) # Change all spines
ax[0].tick_params(width=2, length=8)# increase tick width
ax[0].set_xlabel(r"$\rm{[Fe/H]}$   (dex)", fontsize=fontsize)
ax[0].set_ylabel(r"$\rm{[X/Fe]}$   (dex)", fontsize=fontsize)
#ax[0].set_title(r"$\rm{[O/Fe]}$", fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=fontsize-4)
ax[0].legend(fontsize=fontsize-6, loc="lower left", framealpha=1, ncol=1).set_zorder(-1)
ax[0].axis([-2.5, .5, -0.2, 0.7])
ax[0].set_xticks(np.arange(ax[0].get_xlim()[0], ax[0].get_xlim()[1], .5))
ax[0].axvline(0, 0, 1, color="gray", linestyle="dashed", zorder=-2)# Vertical axis
ax[0].axhline(0, 0, 1, color="gray", linestyle="dashed", zorder=-2)# Horizontal axis
ax[0].text(0.25, 0.6, r"$\rm{O}$", fontsize=fontsize-4)


# Plot [Si/Fe]
ax[1].plot(Abund_MR01_Fe-Solar_values["FeH"]+0.125, Abund_MR01_Si-Abund_MR01_Fe-Solar_values["SiFe"], color="blue", linewidth=2.5, label="MR01", zorder=1)# MR+01 (0.125 approx log10(0.75)
ax[1].plot(Abund_G05Wide_Fe-Solar_values["FeH"]+0.125, Abund_G05Wide_Si-Abund_G05Wide_Fe-Solar_values["SiFe"], color="mediumorchid", linewidth=2.5, label="WIDE G05", zorder=10)# Wide G+05 (0.125 approx
ax[1].plot(Abund_G05Close_Fe-Solar_values["FeH"]+0.125, Abund_G05Close_Si-Abund_G05Close_Fe-Solar_values["SiFe"], color="pink", linewidth=2.5, label="CLOSE G05", zorder=10)# Close G+05 (0.125 approx log10(0.75)
ax[1].plot(Abund_M06_Fe-Solar_values["FeH"]+0.125, Abund_M06_Si-Abund_M06_Fe-Solar_values["SiFe"], color="forestgreen", linewidth=2.5, label="MVP06", zorder=3)# M+06 (Original)
ax[1].plot(Abund_S05_Fe-Solar_values["FeH"]+0.125, Abund_S05_Si-Abund_S05_Fe-Solar_values["SiFe"], color="orange", linewidth=2.5, label="S05", zorder=4)# S+05 (Original)
ax[1].plot(Abund_T08_Fe-Solar_values["FeH"]+0.125, Abund_T08_Si-Abund_T08_Fe-Solar_values["SiFe"], color="deepskyblue", linewidth=2.5, label="T08", zorder=2)# T+08 (Original)
ax[1].plot(Abund_P08_Fe-Solar_values["FeH"]+0.125, Abund_P08_Si-Abund_P08_Fe-Solar_values["SiFe"], color="brown", linewidth=2.5, label="P08", zorder=10)# G+05 (0.125 approx log10(0.75)
# Axis format
for label in ['top','bottom','left','right']: ax[1].spines[label].set_linewidth(2) # Change all spines
ax[1].tick_params(width=2, length=8)# increase tick width
ax[1].set_xlabel(r"$\rm{[Fe/H]}$   (dex)", fontsize=fontsize)
#x[1].set_title(r"$\rm{[Si/Fe]}$", fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=fontsize-4)
#ax[1].legend(fontsize=fontsize, loc="lower left", framealpha=1)
ax[1].axis([-2.5, .5, -0.2, 0.7])
ax[1].set_xticks(np.arange(ax[1].get_xlim()[0]+.5, ax[1].get_xlim()[1]+.5, .5))
ax[1].axvline(0, 0, 1, color="gray", linestyle="dashed")# Vertical axis
ax[1].axhline(0, 0, 1, color="gray", linestyle="dashed")# Horizontal axis
ax[1].text(0.25, 0.6, r"$\rm{Si}$", fontsize=fontsize-4)

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
plt.show()


fig.savefig("Sol_XFe_vsFeH_DTDs.png")


# Plotting ALL [Fe/H](t)
# -------------------------------------------------------------
fig, ax = plt.subplots(1, 1, num="[Fe/H](t)", figsize=(7,7))
ax = [ax]
ax[0].plot( today-t_gyr, Abund_MR01_Fe-Solar_values["FeH"]+0.125,"blue", linewidth=2.5, label="MR01")
ax[0].plot( today-t_gyr, Abund_G05Wide_Fe-Solar_values["FeH"]+0.125,"mediumorchid", linewidth=2.5, label="WIDE G05")
ax[0].plot( today-t_gyr, Abund_G05Close_Fe-Solar_values["FeH"]+0.125,"pink", linewidth=2.5, label="CLOSE G05")
ax[0].plot( today-t_gyr, Abund_M06_Fe-Solar_values["FeH"]+0.125,"forestgreen", linewidth=2.5, label="MVP06")
ax[0].plot( today-t_gyr, Abund_S05_Fe-Solar_values["FeH"]+0.125,"orange", linewidth=2.5, label="S05")
ax[0].plot( today-t_gyr, Abund_T08_Fe-Solar_values["FeH"]+0.125,"deepskyblue", linewidth=2.5, label="T08")
ax[0].plot( today-t_gyr, Abund_P08_Fe-Solar_values["FeH"]+0.125,"brown", linewidth=2.5, label="P08")
# Axis format
for label in ['top','bottom','left','right']: ax[0].spines[label].set_linewidth(2) # Change all spines
ax[0].tick_params(width=2, length=8)# increase tick width
ax[0].set_ylim([-1.5, 0.5])
ax[0].set_xlabel("Age (Gyr)", fontsize=fontsize)
ax[0].set_ylabel(r"$\rm{[Fe/H]}$   (dex)", fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=fontsize-4)
ax[0].set_xticks(np.arange(0,14.1,2) )
ax[0].invert_xaxis()
ax[0].legend(fontsize=fontsize-2, loc="lower right", framealpha=1)
plt.show()

plt.tight_layout()
fig.savefig("Sol_FeH_vs_t_DTDs.png")




# Plotting ALL sigma_Fe(t)
# -------------------------------------------------------------
fig, ax = plt.subplots(1, 1, num="sigma_Fe(t)", figsize=(7.8,7.8))
ax = [ax]
ax[0].plot( today-t_gyr, (Sol_MR01_Fe-Sol_IRA_FE)/Sol_MR01_Fe,"blue", linewidth=2.5, label="MR01")
ax[0].plot( today-t_gyr, (Sol_G05Wide_Fe-Sol_IRA_FE)/Sol_G05Wide_Fe,"mediumorchid", linewidth=2.5, label="WIDE G05")
ax[0].plot( today-t_gyr, (Sol_G05Close_Fe-Sol_IRA_FE)/Sol_G05Close_Fe,"pink", linewidth=2.5, label="CLOSE G05")
ax[0].plot( today-t_gyr, (Sol_M06_Fe-Sol_IRA_FE)/Sol_M06_Fe,"forestgreen", linewidth=2.5, label="MVP06")
ax[0].plot( today-t_gyr, (Sol_S05_Fe-Sol_IRA_FE)/Sol_S05_Fe,"orange", linewidth=2.5, label="S05")
ax[0].plot( today-t_gyr, (Sol_T08_Fe-Sol_IRA_FE)/Sol_T08_Fe,"deepskyblue", linewidth=2.5, label="T08")
ax[0].plot( today-t_gyr, (Sol_P08_Fe-Sol_IRA_FE)/Sol_P08_Fe,"brown", linewidth=2.5, label="P08")
#ax[0].plot(t_gyr, Sol_IRA_FE, "k", linewidth=2.5, label="None", ls="--")
# Axis format
for label in ['top','bottom','left','right']: ax[0].spines[label].set_linewidth(2) # Change all spines
ax[0].tick_params(width=2, length=8)# increase tick width
ax[0].set_xlabel("Age (Gyr)", fontsize=fontsize)
ax[0].set_ylabel(r"$\rm{\sigma_{Fe, Ia}/\sigma_{Fe}}$", fontsize=fontsize)
ax[0].tick_params(axis='both', which='both', labelsize=fontsize)
ax[0].set_xticks(np.arange(0,14.1,2) )
ax[0].invert_xaxis()
ax[0].legend(fontsize=fontsize-4, loc="lower right", framealpha=1, ncol=1)
ax2 = ax[0].twinx()

ax2.set_ylim([1,0.])
ax[0].set_ylim([0,1.])

ax2.tick_params(width=2, length=8)# increase tick width
ax2.set_xlabel("Age (Gyr)", fontsize=fontsize)
ax2.set_ylabel(r"$\rm{\sigma_{Fe, IRA}/\sigma_{Fe}}$", fontsize=fontsize)
ax2.tick_params(axis='both', which='both', labelsize=fontsize)

plt.show()

plt.tight_layout()
fig.savefig("Figure7.png")

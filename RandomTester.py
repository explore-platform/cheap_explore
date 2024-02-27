from CheapTools import *
import numpy as np
import matplotlib.pyplot as plt


# This code performs some tests with random (but numerically reasonable) parameters to
# check if the analytic solution is correct.
Ninfall = 3;

t = np.linspace(0, 14., 30000);

# Generate random model parameters:
kwargs = dict()
kwargs["sigmaX_0"] = np.random.rand()
kwargs["omega"] = np.random.rand()
kwargs["yx"] = np.random.rand()
kwargs["R"] = np.random.rand()
kwargs["nuL"] = np.random.rand()
kwargs["tauj"] = 0.2 + 10*np.random.rand(Ninfall)
kwargs["tj"] = np.sort(14*np.random.rand(Ninfall))
kwargs["Aj"] = np.random.rand(Ninfall)
kwargs["sigma_gas_0"] = np.random.rand()
kwargs["CIa"] = np.random.rand()
kwargs["mx1a"] = np.random.rand()
# --------------------------------------------------------


# Gaussian:
kwargs_G = kwargs.copy()
kwargs_G["AG"] = np.random.rand()
kwargs_G["taup"] = np.random.rand()
kwargs_G["sigma_p"] = 0.01+np.random.rand()# Avoid numerical instabilities due to small numerators
kwargs_G["tau1G"] = 0.1+5*np.random.rand()# Avoid numerical instabilities due to small numerators
kwargs_G["tau2G"] = kwargs_G["tau1G"] + 5*np.random.rand()

print(" Gaussian ")
solr_g, soll_g = LatexCheckEqA9a(t, kwargs_G) # Solve
print(np.max(np.abs(solr_g-soll_g)))
soln_g = R1a_numeric(t, kwargs_G)
print(np.max(np.abs(solr_g-soln_g)))
print(" -------- ")

# Exponential
kwargs_E = kwargs.copy()
kwargs_E["AE"] = np.random.rand()
kwargs_E["tauD"] = np.random.rand()
kwargs_E["tau1E"] = 0.1+5*np.random.rand()
kwargs_E["tau2E"] = kwargs_E["tau1E"] + 5*np.random.rand()
headtail = np.random.rand()<0.5; # True of false:
if headtail:
	kwargs_E["tauD"] = kwargs_E["tauj"][0] # Test this situation:
	print("tauD is tauj[0]")

print(" Exponential ")
solr_e, soll_e = LatexCheckEqA9b(t, kwargs_E) # Solve
print(np.max(np.abs(solr_e-soll_e)))
soln_e = R1a_numeric(t, kwargs_E)
print(np.max(np.abs(solr_e-soln_e)))
print(" -------- ")

# Inverse
kwargs_I = kwargs.copy()
kwargs_I["tauI"] = np.random.rand()
kwargs_I["AI"] = np.random.rand()
kwargs_I["tau1I"] = 0.1+5*np.random.rand()
kwargs_I["tau2I"] = kwargs_I["tau1I"] + (t[-1]-kwargs_I["tau1I"])*np.random.rand()
kwargs_I["tau0"] = 0.8*kwargs_I["tau1I"]*np.random.rand()# Tau0 between 0 and tau1

print(" Inverse ")
solr_i, soll_i = LatexCheckEqA9c(t, kwargs_I) # Solve
print(np.max(np.abs(solr_i-soll_i)))
soln_i = R1a_numeric(t, kwargs_I)
print(np.max(np.abs(solr_i-soln_i)))
print(" -------- ")


dt = np.diff(t)[0]

plt.figure()
plt.subplot(211)
plt.plot(t, solr_g-soll_g, "r", alpha=0.25)
plt.plot(t, solr_e-soll_e, "g", alpha=0.25)
plt.plot(t, solr_i-soll_i, "b", alpha=0.25)

plt.subplot(234)
plt.plot(t, (solr_g-soln_g))

#plt.plot(t, soln_g, "--")

plt.subplot(235)
plt.plot(t, (solr_e-soln_e))

#plt.plot(t, soln_e, "--")

plt.subplot(236)
plt.plot(t, (solr_i-soln_i))

#plt.plot(t, soln_i, "--")
plt.show()

# Now check the goodness of the solutions with ChemicalSolutionVerifier:
# ChemicalSolutionVerifier evaluates NUMERICALLY the derivative of the ANALYTIC 
# solution to check if it is close to the right-hand side of the chemical equation
#
# Note: The discrepancies found so far are attributed to the aproximation of the derivative,
# NOT to the analytic solution itself (they reduce with dt)
error_g = ChemicalSolutionVerifier(t, kwargs_G)
error_e = ChemicalSolutionVerifier(t, kwargs_E)
error_i = ChemicalSolutionVerifier(t, kwargs_I)

plt.figure("Solutions", figsize=(16,9))
plt.plot(t, error_g, "r.", label="Gauss")
plt.plot(t, error_e, "b.", label="Exponential")
plt.plot(t, error_i, "g.", label="Inverse")
plt.legend()
plt.show()

# ChEAP: Chemical Evolution Analytic Package
===============================================================

The **ChEAP (Chemical Evolution Analytic Package)** code implements the analytic solution to the Chemical Evolution Model with TypeIa SNe presented in Palicio et al. (accepted., https://arxiv.org/abs/2304.00042, hereafter P23).
The functions required to compute the solution are contained in the `CheapTools.py` file, which should be imported as a Python library.
We include also the `RandomTester.py` file to illustrate, with a random-parameter chemical evolution model, the accuracy of our analytic solution compared to the numerical integration.

*ChEAP is a Python code created by P.A. Palicio and included in the paper "Analytic solution of Chemical Evolution Models with Type Ia SNe" (P23). If you make use of ChEAP in your work, please consider including the proper citation to this paper. For any question about ChEAP, please do not hesitate to contact the author at __pedro.alonso-palicio(at)oca.eu__*

## 1. The Chemical Evolution equation
-----------------------------------------------------
In this section, we provide a brief overview of the equations we solved analytically. This is not intended to be an exhaustive description of the chemical equation, but rather a general introduction to justify the notation used in the following sections, which focus more on the solution itself. Furthermore, due to the limitations of the markup language, we cannot illustrate these equations with clear notation. Thus, we refer interested readers to the main paper of this code (P23), as well as to the review by Matteucci (2021) and Section 2 of Vincenzo et al. (2017), which motivated this work.

The equation for the evolution of the surface gas density for the X-element (sigma_X) reads:

(d sigma_X)/dt = -alpha*sigma_X(t) + <y_X>(1-R)*psi(t) + <m_XIa>R_{Ia}(t)

or in LateX notation:

 \dfrac{d\sigma_X(t)}{dt} = - \alpha \sigma_X(t)  + \left\langle y_X \right\rangle(1-R)\psi(t) + \left\langle m_{X,Ia} \right\rangle \mathcal{R}\_{Ia}(t)

where the left-most term refers to the variation of the surface density of the X-element per unit time. This variation must be driven by the rates of capture, release, creation, and destruction of the X-element, modeled by the right-hand side terms in the two previous equations:

- The term ~alpha*sigma_X(t) accounts for the loss of X-element due to the star formation and galactic wind, as well as the mass returned to the inter-stellar medium by a single estellar generation. 

- The second term in the right-hand side models the production of X-element, and its posterior contribution to the ISM.

- Similarly, the last term accounts for pollution of the ISM produced by Type Ia SNe. Compared to the previous term, where the pollution was assumed to happen immediately after stellar birth, this contribution requires the inclusion of a time delay that accounts for the formation of the progenitors (generally white dwarfs) and the posterior collapse of the binary system.


## 2. Quick start with ChEAP
-----------------------------------------------------
`CheapTools.py` can be imported as an usual Python library by typing:

**import** CheapTools.py

providing that the `CheapTools.py` file (or a link to it) is placed in the working directory.

Most of the functions in `CheapTools.py`, especially those important for implementing the analytic solution, use a Python dictionary as input.
This dictionary, hereafter referred as **chemdict**, contains the values of the model parameters summarised in Table 2 of P23, including the Galaxy, IRA and DTD(s) parameters.

### 2.1 Basic concepts of Python dictionaries:

Dictionaries are key-value pair structures. Each parameter in a dictionary is defined by a "tag" -- the **key**-- and the **value** associated with that "tag".
To create a new Python dictionary, one may proceed as follows:

***chemdict = dict()# This creates and empty dictionary.***

The **key** is always a string, while the **value** can be almost any type of variable (an integer, a float, an array, a list, another dictionary, etc.).</br>
For example, the following dictionary has a **key** named "omega" whose **value** is set to 0.4:

***chemdict["omega"] = 0.4***

Within the same dictionary, we can save as many as key-values pairs as we want:
***chemdict["tauj"] = [7.]***

The previous lines are equivalent to the following definition of a Python dictonary:
***chemdict = {"omega":0.4, "tauj":[7.]}***

**IMPORTANT**: Dictionaries should be copied as:

dict2 = dict1.copy()</br>
and not as </br>
dict2 = dict1 </br>
otherwise the modification of dict2 would change dict1.

For more information about python dictionaries check https://docs.python.org/3/tutorial/datastructures.html


### 2.2 Keys of the ChEAP chemical dictionary:
Most of the ChEAP functions use a Python dictionary as input for parameters.</br>
The nomenclature of the ***keys*** of this dictionary is summarised as follows:

|key|value type|Description (see Table 2 in P23)|
|---------------|-------------------------------|-------------------------------|
|tauj|array of floats| Time-scale of gas accretion for each infall|
|tj|array of floats| Starting time of each gas infall|
|Aj|array of floats| Amplitude of each infall. Tunes the amount of accreted gas.|
|---------------|-------------------------------|-------------------------------|
|omega|float|Wind loading factor.|
|nuL|float|Star formation efficiency.|
|sigmaX_0|float|initial content of element X: sigma\_X(0)|
|sigma_gas_0|float|initial amount of gas|
|R|float|Recycling fraction|
|yx|float|Yield per stellar generation of the element X.|
|mx1a|float|Average amount of X synthesized by each single Type Ia SN event.|
|CIa|float|Normalisation constant for the TypeIa SNe rate.|
|---------------|-------------------------------|-------------------------------|
|taup|array of floats| Offset of the Gaussian DTD (tau' in Table 2)|
|sigma_p|array of floats| Width of the Gaussian DTD (sigma' in Table 2)|
|AG|array of floats| Amplitudes of the Gaussian DTD (A\_G in Table 2)|
|tau1G|array of floats| Lower time limit for a TypeIa SNe assuming the Gaussian DTD|
|tau2G|array of floats| Upper time limit for a TypeIa SNe assuming the Gaussian DTD (tau2G>tau1G).|
|---------------|-------------------------------|-------------------------------|
|tauD|array of floats|Time-scale of the exponential DTD.|
|AE|array of floats|Amplitudes of the exponential DTD (A\_E in Table 2)|
|tau1E|array of floats|Lower time limit for a TypeIa SNe assuming the exponential DTD|
|tau2E|array of floats|Upper time limit for a TypeIa SNe assuming the exponential DTD (tau2E>tau1E)|
|---------------|-------------------------------|-------------------------------|
|tau0|array of floats|Offset of the 1/t DTD: 1/(t-tau0). It must be lower than tau1I.|
|tauI|array of floats|Characteristic time (just for the units). Set to 1 Gyr.|
|AI|array of floats|Amplitudes of the 1/t DTD (A\_I in Table 2).|
|tau1I|array of floats|Lower time limit for a TypeIa SNe assuming the 1/t DTD.|
|tau2I|array of floats|Upper time limit for a TypeIa SNe assuming the 1/t DTD (tau2I>tau1I)|

The numeric values of the last three groups of parameters will be available at CDS for all the DTDs considered in P23.

### 2.3 Computing the exact solution
Once the parameters of the model are "encapsulated" in the chemical dictionary **chemdict**, one has to specify the time values (in Gyr) for evaluating the solution (i.e., the **time nodes array**).</br>
Generally, this array is a uniform sampling of the time interval of interest like t = np.linspace(0, 13.8, 1000) or t = np.arange(0, 13.8, 0.1), but the functions of **CheapTools** accept non-uniform (and even non-sorted) input times.</br>
Assuming "t" is the array containing the input time and **chemdict** the chemical dictionary, the solution sigma_X (density of the X-element) is computed by the function **SolveChemEvolModel** using:

Sigma_X = **SolveChemEvolModel(** t, chemdict **)**

However, it is more natural to express the amount of metals referred to the hydrogen in log. scale as \[X/H\]:

X_over_H = **FromSigmaToAbundance(**t_, Sigma_X, chemdict **)**  + 0.125 

where Sigma_X is the output of **SolveChemEvolModel** and the additional term **0.125** corresponds to log10(0.75), since the hydrogen constitutes approximately the 75% of the primordial gas.</br>
Finally, an additional term has to be added so that the solar values make X_over_H = 0. This **scaling factor** --named "SV" in Section 3.1 of P23-- varies from one element to another and can be computed using the solar values reported in Table 1 of Asplund et al. (2009).


### 2.4 Verification of the analytic solution
`CheapTools.py` provides the function **ChemicalSolutionVerifier** to evaluate numerically the analytic solution in the model equation.</br>
Basically, this function compares the left and right-hand side terms of the differential equation for sigma_X (Eq.1 in this file and Eq. 9 in P23). The time derivative of sigma_X (the left-hand side term) is approximated by the finite-difference derivate implemented in the numpy function "gradient". This returns a reasonably good approximation for the exact derivative when the input time-node array is dense enough, while we must exclude some specific nodes where that approximation is expected to be worse:

* The first node (no left neighbour node available).
* The last node (no right neighbour node available).
* The left and right neighbour nodes of each infall (the derivative is not continuous at "tj").

The function **ChemicalSolutionVerifier** requires to inputs:

1. The **time-node array** (in Gyrs), ideally with a small step size so the numerical derivative is accurate enough.
2. The **chemical dictionary** that defines the model. It will be used internally by **ChemicalSolutionVerifier** to compute the exact solution and evaluate it in the right-hand side of Eq. 9.

If the analytic solution is correct, the output is an array of values very small number (in absolute value). Up to this point, all relatively significant discrepancies relative to the zero value encountered can be ascribed to an excessively large time step size.


### 2.5 Shortcuts
* Once the chemical dictionary is constructed, we recommend to use the **prepare\_chemdict** function before the first use to check if there are some inconsistencies: **chemdict = prepare\_chemdict( chemdict )**
* If all the DTDs have the same tau1G, tau1E and tau1I, then we can simply use the key "tau1" for all of them. Then, **prepare_\_chemdict** will fill tau1G, tau1E and tau1I with the provided tau1. The same stands for tau2.
* There are some functions for loading the parameters of the more complex DTDs:
	* **Load\_MR01\_dict()** for the **Matteucci and Recchi (2001)** DTD.
	* **Load\_G05Wide\_dict()** for the **WIDE Greggio (2005)** DTD.
	* **Load\_G05Close\_dict()** for the **CLOSE Greggio (2005)** DTD.
	* **Load\_P08\_dict()** for the **Pritchet et al. (2008)** DTD.
* If a dictionary is provided as input for these functions, the output dictionary will contain the input key-value pairs.
* The value of the parameter "CIa" can be defined by the present-day TypeIa SN rate by using the **Get_CIa** function.


## 3. Examples of ChEAP usage:
-----------------------------------------------------
Run the file `QuickTest.py` for an example with the evolution of iron, and `Example_Fig7.py` for reproducing Fig. 7 in P23.


## 4. References:
-----------------------------------------------------
1. Palicio et al. (2023, accepted) [P23]: https://ui.adsabs.harvard.edu/abs/2023arXiv230400042P/abstract
2. Matteucci (2021): https://ui.adsabs.harvard.edu/abs/2021A%26ARv..29....5M/abstract
3. Vincenzo et al. (2017): https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.2939V/abstract
4. The table with the numeric values of the parameters of the DTD fitting (Table 1 in P23) is available at CDS with larger decimal precision than the printed version: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/678/A61#/browse 
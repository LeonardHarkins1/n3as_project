# n3as_project

## Task 1, warm up
Calculate integral 1
Notes: This integral diverges as $\Lambda\to \infty$.

* Chose several values of $\Lambda >> m$, also keep $E << \Lambda$.
* Chose several values of $\epsilon << m$. In principle we want to take the limit $\epsilon\to0$, you might want to check what happens if you take $\epsilon=0$ and try to perform the integral.
* Try to see if you find something in common among integrals with different $\Lambda$, i.e. something independent of $\Lambda$.
* (April 16) Show the results for several values of $\Lambda$ in the same plot, have plots for real and imaginary parts of the integral.

### Possible follow ups (April 16)
* Try to explain the behavior of the plots directly from the analytic expression, e.g. is there a qualitative difference in the integrand for values of $E$ above/below the cusp that we saw around $E\approx 2$?
* Play around with the options of the NIntegrate function of mathematica, do changes to Accuracy/Recursion/Points/Precision alleviate the fuzziness?
* Can you think a way to calculate the $\epsilon\to0$ limit of the integral?

* Also if you have some time, begin writing a summary of the literature I've shared, you can put this in the overleaf, have it after a new "\section{Literature review}".

### Follow up from April 23
* Fix the factors of 2 in the analytic derivation, verify your analytic answer by plotting it together with the imaginary part of the integral to verify if they agree.
* Add another row of plots where you remove the "real" part of the square root from the real part of the integral to check if the cusp dissapears.
* Read more of last the paper I sent: RevModPhys.90.025001.pdf 


## Task 2
Play around with the integral

$$
I(E_f,E_i,q_0,m)=
\int \frac{ dk_0 d^3\vec{k} }{(2\pi)^4}
\left[\frac{1}{(E_f-k_0)^2-\omega_k^2+i\epsilon}\right]
\left[\frac{1}{(q_0+E_f-k_0)^2-\omega_k^2+i\epsilon}\right]
\left[\frac{1}{(E_i-k_0)^2-\omega_k^2+i\epsilon}\right]
\left[\frac{1}{k_0^2-\omega_k^2+i\epsilon}\right]
$$

Note that in the previous case you had only $E$ as an independent variable, in this case you have $E_f$, $E_i$ and $q_0$.
However, in this case the integral converges, and there is no need to introduce an artificial $\Lambda$ cut off for the integration region.

### Follow up from May 9
* You can perform the contour integral following the procedure in the notes sent on the email today, you should get 4 terms.
* Check which of these terms appear to have poles where the $\delta$, in the $\frac{1}{x-i\epsilon}=PV\frac{1}{x} + i\pi\delta(x)$, will be relevant.
* Probably some plots to illustrate the point above would be a good idea.

### Follow up from May 23
* Split up the integral into the three conserving symmetry terms ($E_f \leftrightarrow E_i, q_0 \rightarrow q_0 + E_f - E_i$).
* Check if each singularity on the integral comes from different terms, or if all singularities show up for each of them.
* Plot the denominators of each term to see which will require the Principal Value plus delta Dirac. Keep in mind the relevant kinematic regions, e.g. $4m> Ef > 1.8m$, etc.

### Follow up from June 4
* Perform the integral for each of the symmetry conserving terms.
  * Begin wtih the term where only one of the denominators requiring the principal value prescription.
  * Then try the one with the two PV denominators.
  * Finally try the term for which all three denominators vanish for some values of the integral
(Example of numerical integral in mathematica with principal value (epsilon=0))
* Begin reading chapters 1 and 3 of the dissertation

### June 18
NOTE: Can you get samples of the honor thesis proposal?
* Double check the derivations, we expect the real part of the numerical integral to be equal to the result of the delta diracs.

### July 3
* Get the Real part of the integral for the case of the other diagram we have not considered, this would amount to replace $q_0$ with $E_f$ in the last propagator.
* As a bonus, to better understand the spikes in the numerical integral, you can try and plot the lines where the double delta diracs would appear.
* A few plots to observe the dependence of $E_f$ in the integral would be useful, keep plots as a function of $E_i$, but also choose four or five values of $E_f$.

### July 9
* Give a physical interpretation to the singularities you observe around the physical region, with $E_i$ and $E_f$ in between $\sim 1.8$ and 4.
* Probably this would mean isolating different terms in your final result, so that you observe the different singularities associated with each term.

### July 16
* Summer school lectures, scattering
  
### July 23
* Identify singularties of the integral with specific propagators diverging in the integrand (cuts in the diagram).
* Sketch in a $E_i,E_f$ plane, with fixed $q_0$ the location of each of these singularities.


### October 30
* Finish up learning about relativity (up to Chapter 12.2.4 in Griffiths)
* Begin plotting form factors as a function of $Q^2$, begin with the axial form factor.
* Tighter constraints in variable $p_h^2$.

### November 21
* Finish the implementation of F_A, with creating the variable sample only once, instead of once per $Q^2$.
* Problem 12.36 of Griffiths electrodynamics. Add extra difficulty by comparing results from boosting to the CM frame.
* Begin looking at the implementation of $F_V$, $A$, $B$, $C$, and $R$.



Mathematica useful commands:
Save plots with Legends
SetDirectory[NotebookDirectory[]]
Export["plot.pdf", plot]

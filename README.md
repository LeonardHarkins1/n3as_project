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


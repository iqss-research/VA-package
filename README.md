# VA: R package for verbal autopsies

Gary King and Ying Lu

## About
VA is an easy-to-use R program that automates the analysis of verbal autopsy data. These data are widely used for estimating cause-specific mortality in areas without medical death certification.

Data on symptoms reported by caregivers along with the cause of death are collected from a medical facility, and the cause-of-death distribution is estimated in the population where only symptom data are available. Current approaches analyze only one cause at a time, involve assumptions judged difficult or impossible to satisfy, and require expensive, time consuming, or unreliable physician reviews, expert algorithms, or parametric statistical models. By generalizing current approaches to analyze multiple causes, King and Lu (2008) show how most of difficult assumptions underlying existing methods can be dropped. These generalizations, which we implement here, also make physician review, expert algorithms, and parametric statistical assumptions unnecessary. While no method of analyzing verbal autopsy data can give accurate estimates in all circumstances, the procedure offered is conceptually simpler, less expensive, more general, as or more replicable, and easier to use in practice.

More generally, the software takes as input a multicategory variable D, and a set of dichotomous variables S (cause of Death and Symptoms, respectively, in verbal autopsy applications). Both variables exist in one data set (a hospital in the application) but only S exists in the population of interest. The goal of the procedure is to estimate the probability distribution (or histogram) of D in the population of interest.

For more information, see Gary King and Ying Lu. Verbal Autopsy Methods with Multiple Causes of Death (Abstract: <a href="https://garyking.org/files/abs/vamc-abs.shtml">HTML</a>  | Paper: <a href="https://garyking.org/files/gking/files/vamc.pdf">PDF</a>) 

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/us/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/3.0/us/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/us/">Creative Commons Attribution-NonCommercial-NoDerivs 3.0 United States License</a>.

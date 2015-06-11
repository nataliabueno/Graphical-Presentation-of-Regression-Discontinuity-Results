# Graphical-Presentation-of-Regression-Discontinuity-Results
R Code - "Graphical Presentation of Regression Discontinuity Results"  (Bueno &amp; Tuñón, 2015)

This repository has the documented R code for presenting RDD results in a simple graphical way that shows the sensitivity of estimates to a wide range of possible windows.

Abstract: During the last decade, an increasing number of political scientists have turned to regression-discontinuity (RD) designs to estimate causal effects. One of the key issues in RD designs involves selection of the window or bandwidth as results are often sensitive to their size. We propose a simple graphical way of reporting RD results that shows the sensitivity of estimates to a wide range of possible windows. The advantages of these plots over conventional ways of presenting RD results are threefold. First, they reduce the opportunities for the intentional or unintentional selection of windows that yield significant findings while also communicating the robustness of results to variations in window size. Second, they provide a simple visual comparison of how different estimators (difference-of-means, local linear, and polynomial regression models) perform as a function of different windows. Third, these plots can be used for placebo tests to examine the sensitivity of balance to different windows, reporting the magnitude of the difference between treatment and control and its confidence interval. We illustrate how a researcher can use these plots to present RD results in a transparent and systematic way with a regression discontinuity design examining incumbency advantage in the U.S. House of Representatives.

Published in [*The Political Methodologist*] (http://thepoliticalmethodologist.com/) (blog and print edition)

\[Note on 06/10/2015\]: We corrected a bug on line 86 of the code. We fixed the RunningCenTreat object in the local_linear function. Thanks to Daniel Masterson for pointing it out to us. 

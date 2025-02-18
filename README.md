# DVSG
A Python module to calculate the $\Delta V_{\star-g}$ (pronounced 'DVSG') value of a galaxy.

## Overview
Following **Powley, Smethurst, and Lintott (in prep.)**, a galaxy's $\Delta V_{\star-g}$ value is defined as:

$\Delta V_{\star-g} = \frac{\sum_{j} \left| {V_{\star,\text{norm}}^{j} - V^{j}_{g,\text{norm}}} \right|}{N}$

where:
- $V_{\star,\text{norm}}$ is the **normalised stellar velocity map**,
- $V_{g,\text{norm}}$ is the **normalised gas velocity map**,
- $\sum_{j} \left| {V_{\star,\text{norm}}^{j} - V^{j}_{g,\text{norm}}} \right|$ is the **sum** over all spaxels, $j$, of the **absolute difference** between the normalised stellar and gas velocity maps
- $N$ is the **number of spaxels** contributing towards the sum

A $\Delta V_{\star-g}$ value of 0 would imply no difference in the kinematics of stellar and gas velocity, whereas a $\Delta V_{\star-g}$ value of 1 suggests the largest offset between the stellar and gas velocity. The general trend is that the more kinematically disturbed galaxies tend to have larger $\Delta V_{\star-g}$ values. For more information about $\Delta V_{\star-g}$, please refer to Powley, Smethurst and Lintott (in prep.).

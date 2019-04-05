# OXFORD_ASL 

A command line tool for quantification of perfusion from ASL data

`oxford_asl` is part of the BASIL toolbox within FSL for the analysis of 
perfusion ASL data. `oxford_asl` provides a single means to quantify CBF from ASL 
data, including kinetic-model inversion, absolute quantification via a 
calibration image and registration of the data. It provides most of the common 
options that someone who has raw ASL data might like to perform to extract 
perfusion images and thus is the normal place to begin for most users who want 
a command line tool. 

`oxford_asl` is the main underlying tool behind the 
equivalent graphical user interface: `Asl_gui`. To use `oxford_asl` you may find 
you first need to pre-process your data using `asl_file` (see the tutorial for 
examples where this is relevant). More advanced users wishing to do customised 
kinetic analysis might want to use the basil command line tool - this is the 
very core of `oxford_asl` and thus the `BASIL` toolbox overall.

For full documentation on the `oxford_asl` suite of tools see:

https://asl-docs.readthedocs.io/

## Referencing

If you use oxford_asl in your research, please reference the article below, plus
any others that specifically relate to the analysis you have performed:

*Chappell MA, Groves AR, Whitcher B, Woolrich MW. Variational Bayesian inference 
for a non-linear forward model. IEEE Transactions on Signal Processing 
57(1):223-236, 2009.*

If you employ spatial priors you should ideally reference this article too.

*A.R. Groves, M.A. Chappell, M.W. Woolrich, Combined Spatial and Non-Spatial 
Prior for Inference on MRI Time-Series , NeuroImage 45(3) 795-809, 2009.*

If you fit the macrovascular (arterial) contribution you should reference this 
article too.

*Chappell MA, MacIntosh BJ, Donahue MJ, Gunther M, Jezzard P, Woolrich MW. 
Separation of Intravascular Signal in Multi-Inversion Time Arterial Spin 
Labelling MRI. Magn Reson Med 63(5):1357-1365, 2010.*

If you employ the partial volume correction method then you should reference 
this article too.

*Chappell MA, MacIntosh BJ, Donahue MJ,Jezzard P, Woolrich MW. Partial volume 
correction of multiple inversion time arterial spin labeling MRI data, Magn 
Reson Med, 65(4):1173-1183, 2011.*

If you perform model-based analysis of QUASAR ASL data then you should reference 
this article too.

*Chappell, M. A., Woolrich, M. W., Petersen, E. T., Golay, X., & Payne, S. J. 
(2012). Comparing model-based and model-free analysis methods for QUASAR 
arterial spin labeling perfusion quantification. doi:10.1002/mrm.24372*

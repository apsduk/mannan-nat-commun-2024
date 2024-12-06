README
-----------------------------------------------------------------------------------
MATLAB CODE FOR FINDING OPTIMAL CHANGES IN GENE EXPRESSION  (STATIC DESIGNS)


This folder contains all the MATLAB code (developing in MATLAB R2023b) for
solving the multiobjective optimisation problems maximising:
  		-  objectives (A): growth rate and synthesis flux,
  		-  objectives (B): volumetric productivity and yield.

 These generate the plots in
 		-  Figure 1
 		-  Supplementary Figure 1
-----------------------------------------------------------------------------------


The following code scripts do the following:

 -  "A1_RUN.m"
 	This script runs "A1_MultiobjOpt.m" to regenerate all plots in Figure 1 of 
 	Results subsection 2.1.

 -  "A1_MultiobjOpt.m"
 	Script that is doing all analyses, i.e.:
 	 -  defining the bounds on the scaling of the transcription rates of Ep, Tp, E
 	 -  solving the multiobjective optimisation with objectives (A)
 	 	 -  calculating the vol. productivity and yields for all optimal solutions
 	 		of optimisation (A)
 	 	 -  plotting the Pareto optimal solutions of optimisation (A)
 	 -  solving the multiobjective optimisation with objectives (B)
 	 	 -  calculating the growth rate and synthesis flux for all optimal solutions 
 	 		of optimisation (B)
 	 	 -  plotting the Pareto optimal solutions of optimisation (B)


 -  All other scripts with names beginning "m_... .m" are used in "A1_MultiobjOpt.m"
 	or within the other scripts described here:
 	
 	 -  "model_full.m"
 	 	the full host aware model
 	 
 	 -  "model_params.m"
 	 	the defined model parameters
 	 
 	 -  "m_steadystate.m"
 	 	script to simulate the single cell system till it reaches steady state, to
 	 	within a tolerance of 10^{-5}

 	 -  "m_ssFor0pwExp.m"
 	 	script to simulate the single cell system without the heterologous synthesis
 	 	pathway, till steady state

 	 -  "m_BatchCultModel.m"
 	 	script to simulate batch culture of a growing population till all initial
 	 	set concentration of feed reaches 0, for calculating the vol. productivity
 	 	and yield

 	 -  "m_MultiObjs_1.m"
 	 	script calculating objectives (A)

 	 -  "m_MultiObjs_2.m"
 	 	script calculating objectives (B)

 	 -  "m_Event_SubsRunOut.m"
 	 	events function defining that the batch culture simulations should stop when
 	 	all feed substrate is consumed to a concentration ≤ 1nM

 	 -  "m_SingleObjs_2.m"
 	 	script to calculate the maximum value of each of the objectives, one at a 
 	 	time, either vol. productivity or yield

 	 -  "m_TXreg.m"
 	 	script defining the regulation/expression of the transporters
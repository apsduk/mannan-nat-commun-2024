% ======================================================================= %
% EXAMPLE CODE FOR REVIEW                                                 %
% ======================================================================= %
%
% ======================================================================= %
% Ahmad A. Mannan and Alexander P.S. Darlington                           %
% ======================================================================= %

% ----------------------------------------------------------------------- %
% Scripts                                                                 %
% ----------------------------------------------------------------------- %

>> RUN_model_example.m - Simulates the biotechnological model.
 
>> RUN_multiobj_example.m - Sets up and solves the multiobjective
optimisation problem and saves the Pareto fronts (as objVals) and 
designs (as pOpts).

% ----------------------------------------------------------------------- %
% Functions                                                               %
% ----------------------------------------------------------------------- %
>> model_params.m
    Description: Contains the parameters for the host and circuit.
    Inputs:
        - sS0 (energy yield per iS molecule)
        - vX0 (vmax of promiscuous host exporter)
        - KmX0 (K_M of promiscuous host exporter)
        - VolCult0 (Volume of cutlure)
    Outputs:
        - h_params (host parameters)
        - x_params (circuit/pathway parameters)

>> writejacobianfile.m
    Description: Symbolic constructs the jacobian of the system described 
    by the inputs and writes it to the file "fname" which can be supplied 
    to ode15s for enhanced accuracy and efficiency.
    Inputs:
        - odesourcename (name of the ODE file as string)
        - fname (name of the output file as string)
        - T (vector to specifify size of T, e.g. 1x1 double)
        - Y (vector to specifify size of Y, e.g. 1xlength(Y0) double)
        - hPR (host parameters)
        - xPR (circuit/pathway parameters)
        - y0scales (size of Y0scales variable needed by ode function)
        - SysTopol (system topology as a vector)
    Outputs:
        - fname (name of the output file as string)
        - J (symbolic Jacobian)

>> m_BatchCultModel_DC.m
    Description: ODE function
    Inputs:
        - T (current time point for ode solver)
        - Y (current Y vector for ode solver)
        - hPR (host parameters)
        - xPR (circuit/pathway parameters)
        - y0scales (not used)
        - SysTopol (system topology as a vector)
    Outputs:
        - dY (vector of ODEs)
        - lambda (current growth rate)
        - TXrates (vector of current transcription rates)
        - TLrates (vector of current translation rates)
        - gammaX (current peptide elongation rate)
        - protmass (vector of current proteome in terms of mass)
        - summass (current cell mass)
        - ribomassfrac (current mass fraction of ribosomes)
        - fluxes (vector of current fluxes)

>> m_BatchCultModel_SS.m
    Description: Cell only model where dN/dt = dxS/dt = dxP/dt = dxI/dt = 0.
    Inputs: See m_BatchCultModel_DC.m.
    Outputs: See m_BatchCultModel_DC.m.

>> m_BuildBounds.m
    Description: Builds vectors for bounds using in optimsiations.
    Inputs: 
        - SysTopol (topology)
        - sTX_Z (vector of bounds for transcription scaling 
          parameters, [lb ub])
        - K_Z (vector of bounds for K threshold values, [lb ub])
        - t_ind (vector of bounds for induction time, [lb ub])
    Outputs:
        - lb (vector of lower bounds)
        - ub (vector of upper bounds)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - xnames (names of parameters which are varied)
        - topodescr (description of topology in human readable format) 

>> m_defSysAndTopols.m
    Description: Make matrix of all topologies to consider for optimization
    Inputs: None.
    Outputs:
        - SysTopols (matrix of topologies, each row is a topology to be 
          tested)

>> m_MultiObj_v2.m
    Description: Cost/objective function for the 2D optimisation.
    Inputs:
        - p (current parameter set to simulate, equiv. to X in multiga)
        - hPR0 (host parameters)
        - xPR0 (circuit/pathway parameters)
        - Y0scales (not used)
        - SysTopol (topology)
        - tmax (maximium value of time for simulation)
        - sfac (max values of objectives used to scale optimisation values 
          between 0 and 1, see Methods)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - odesolver (choice of ODE solver)
        - odeoptions (options for ODE solver)
    Outputs:
        - scaledObjs (objective values scaled by sfac)
        - Objs (objective values for parameter set p)
        - T (time vector from ODE solver)
        - Y (Y vector from ODE solver)
        - flag (output flag)
        - tsol (time vector from ODE solver)
        - ysol (Y vector from ODE solver)

>> m_ProdYieldObj_v2.m
    Description: Objective function for Yield optimisation.
    Inputs:
        - p (current parameter set to simulate, equiv. to X in multiga)
        - hPR0 (host parameters)
        - xPR0 (circuit/pathway parameters)
        - Y0scales (not used)
        - SysTopol (topology)
        - tmax (maximium value of time for simulation)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - odesolver (choice of ODE solver)
        - odeoptions (options for ODE solver)
    Outputs:
        - Objs (yield output for parameter set p)
        - T (time vector from ODE solver)
        - Y (Y vector from ODE solver)
        - tsol (time vector from ODE solver)
        - ysol (Y vector from ODE solver)

>> m_SimBatchCult_ProdPerf.m
    Description: Simulates the batch culture of system with dynamic
    controller and measures performance metrics.
    Inputs:
        - hPR (host parameters)
        - xPR (circuit/pathway parameters)
        - y0scales (not used)
        - SysTopol (topology)
        - tind (induction time)
        - tmax (maximium value of time for simulation)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - odesolver (choice of ODE solver)
        - odeoptions (options for ODE solver)
    Outputs:
        - ProdPerf (performance output, [vProd, pYield])
        - T (time vector from ODE solver)
        - Y (Y vector from ODE solver)
        - tsol (time vector from ODE solver)
        - ysol (Y vector from ODE solver)

>> m_SimBatchCult_SystDynamics.m
    Description: Simulates the batch culture of system with dynamic
    controller and measures performance metrics. Also returns all internal
    dynamics.
    Inputs:
        - hPR (host parameters)
        - xPR (circuit/pathway parameters)
        - y0scales (not used)
        - SysTopol (topology)
        - tind (induction time)
        - tmax (maximium value of time for simulation)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - odesolver (choice of ODE solver)
        - odeoptions (options for ODE solver)
    Outputs:
        - ProdPerf (performance output, [vProd, pYield])
        - lambda (current growth rate)
        - TXrates (vector of current transcription rates)
        - TLrates (vector of current translation rates)
        - gammaX (current peptide elongation rate)
        - protmass (vector of current proteome in terms of mass)
        - fluxes (vector of current fluxes)
        - summass (current cell mass)
        - T (time vector from ODE solver)
        - Y (Y vector from ODE solver)
        - dY_by_dt (derivatives over time)
        - flag (quality of simulation)

>> m_TXreg.m
    Description: Models the impact of the transcription factor pTF on the 
    regulated gene.
    Inputs:
        - params (parameter values)
        - pTF (transcription factor molecules per cell)
    Outputs:
        - r_g2m (transcription rate)

>> m_VolProdObj_v2.m
    Description: Objective function for Vol Prod. optimisation.
    Inputs:
        - p (current parameter set to simulate, equiv. to X in multiga)
        - hPR0 (host parameters)
        - xPR0 (circuit/pathway parameters)
        - Y0scales (not used)
        - SysTopol (topology)
        - tmax (maximium value of time for simulation)
        - xidx (idx of parameters which are varied, set by 
          topology SysTopol)
        - odesolver (choice of ODE solver)
        - odeoptions (options for ODE solver)
    Outputs:
        - Objs (yield output for parameter set p)
        - T (time vector from ODE solver)
        - Y (Y vector from ODE solver)
        - tsol (time vector from ODE solver)
        - ysol (Y vector from ODE solver)

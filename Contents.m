%Computational simulation model of the early visual system
%
%DESCRIPTION
%
%This framework uses four core functions to run a leaky integrate-and-fire
%neuron model simulation of (at this moment) monkey V1 and V2. Work is
%underway to widen its applicability. The four functions are:  
%
%1) runBuildConnectivity.m
%2) runBuildStimulation.m 
%3a) runSimulation.m AND 3b) runSimulationIndInp.m
%4) runSimBatchAggregator.m
%
%The program runBuildConnectivity.m (1) creates a .mat file containing the
%synaptic connectivity and simulation parameters. It in turn uses the
%subfunction buildConnectivity() to create the actual connectivity matrix.
%You can set the simulation parameters and model connectivity by editing
%the fields of the sConnParams structure in this file. See the function
%buildConnectivity() for more information. This function effectively
%creates the "anatomy" of your model network, which you can then use to run
%a simulation.
%
%The second script that is required to run the simulation is
%runBuildStimulation (2), which creates the input that will feed into the
%network through a simplified retina/LGN filter layer. This function can
%create visual stimuli in retinal degrees, just as you would for an actual
%visual neuroscience experiment.
%
%The main workhorse is accessed through runSimulation.m (3a) or
%runSimulationIndInp.m (3b). These functions use string-based inputs with
%parameter delimited by commas so compiled version can also be started from
%the command line. The difference between the two functions is that (3a)
%will simulate a single retina, which means that shared noise will
%propagate to cortex and induce information-limiting correlations. If you
%want to remove this source of shared fluctuations, you will have to use
%(3b). The drawback of (3b) is that simulating independent inputs for all
%V1 cells is computationally very intensive, so it will noticeably slow
%down the simulation speed.
%
%Finally, if you run the simulation on a cluster with independent workers
%in a massively parallel manner, you can use runSimBatchAggregator.m (4) to
%combine the separate simulation runs from a single connectivity structure
%into one data file. See the runSimBatchAggregator.m file for more
%information.
%
%VERSION HISTORY
%
%V0.9.0.20190204 
%These functions were written by Jorrit Montijn in Alexandre Pouget's lab
%in Geneva, over the course of 2016-2018. 
%
%V0.9.1.20190212 
%Started to update help functions and rename functions to things that are
%more concise, and/or make more sense [by JM]. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mutual repressor switch simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
la = 0;   % TALA-KRAB leakage fraction (0-1) 
lb = 0;   % TALb_KRAB leakage fraction (0-1)
n1 = 1.0; % Cooperativity of binding to construct 1
n2 = 1.0; % Cooperativity of binding to construct 2
n3 = 1.0; % Cooperativity of binding to construct 3
n4 = 1.0; % Cooperativity of binding to construct 4

% Switching pattern is denoted by e.g. [1 0 0; 0 0 0],
% where time is divided in three frames (t1, t2, t3) and 
% the first row represents BFP activation and second mCITRINE.

% Example
% Switching from BFP in t1 to mCITRINE in t2
% [1 0 0; 0 1 0] stands for 
%	[BFP(t1), 0 , 0 ; 
%	0 , MCITRINE(t1) ,0]

% Run the simulations
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [1 0 0; 0 0 0]); save -ascii output/MutRep_10.txt A;
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [0 0 0; 1 0 0]); save -ascii output/MutRep_20.txt A;
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [1 0 0; 0 1 0]); save -ascii output/MutRep_12.txt A;
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [0 1 0; 1 0 0]); save -ascii output/MutRep_21.txt A;
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [0 0 0; 0 0 0]); save -ascii output/MutRep_000.txt A;
A = MutualRepressorSimulation(n1,n2,n3,n4,la,lb, [1 0 0; 1 0 0]); save -ascii output/MutRep_122.txt A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Positive feedback loop switch simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [1 0 0; 0 0 0]); save -ascii output/PFS_10.txt A;
A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [0 0 0; 1 0 0]); save -ascii output/PFS_20.txt A;
A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [1 0 0; 0 1 0]); save -ascii output/PFS_12.txt A;
A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [0 1 0; 1 0 0]); save -ascii output/PFS_21.txt A;
A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [0 0 0; 0 0 0]); save -ascii output/PFS_000.txt A;
A = PositiveLoopSwitchSimulation(n1,n1,n1,n1,n1,n1,la,lb, [1 0 0; 1 0 0]); save -ascii output/PFS_122.txt A;

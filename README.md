========
SWITCH-IT 
=========

The repository provides scripts, used for modeling a bistable switch based on designed DNA binding elements in mammalian cells.

Reference:
Lebar, Tina, Urban Bezeljak, Anja Golob, Miha Jerala, Lucija Kadunc, Boštjan Pirš, Martin Stražar, Roman Jerala et al. "A bistable genetic switch based on designable DNA-binding domains." Nature communications 5 (2014).


DETERMINISTIC MODEL
===================

Prerequisites
-------------
The scripts for the deterministic ODE model are written in MATLAB/Octave.

Quick start
-----------
We provide models for classic mutual repressor switch (toggle) and the positive feedback loop switch.

Run the simulations by running the script
	```octave
	run.m
	```

inside appropriate MATLAB/Octave interpreter.

Results are saved the the folder 'output' and are store in tab-separated files storing the concentrations of both reporters:
	```
	TIME	BFP	mCITRINE
	```

STOCHASTIC MODEL
================

Prerequisites
-------------
The currect version is written for Unix/Linux/Mac OS platforms.
To enable running on other OS, code in ssatest.c shall be modified when referecing files.

In order to compile the files, gcc or other compiler is required.


Quick start
-----------
To run a stochastic model, cd to the stochastic directory
	```
	cd stochastic
	```

Compile all models
	```
	chmod 755 make.sh
	./make.sh
	```

Run all models
	```
	chmod 755 run.sh
	```

The results will be stored in the folder output, containing the tab-separated files with rows of reporter concentrations at each iteration and columns defined as follows:
	```
	BFP	mCITRINE	REACTION_NUMBER	TIME
	```

Modelling
---------
All models are defined as hardcoded C arrays containing the stoichiometry of chemical species. Chemical reactions and species are asigned to unique rows and columns of the stoichiometry matrix, respectively.

The following models are implemented:
	```
	model_minfeed.h    (Bistable switch with positive feedback loop w/ minimal promoters)
	model_cmvfeed.h    (Bistable switch with positive feedback loop w/ constitutive promoters)
	model_nocomp.h	   (Bistable switch with positive feedback loop w/ minimal promoters, no competition for binding sites)
	model_noloop.h	   (Bistable switch with positive feedback loop w/ constitutive promoters, no positive feedback)
	model_toggle.h     (Classic toggle switch with no positive feedback)
	```


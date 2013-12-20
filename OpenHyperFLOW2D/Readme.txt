                     OpenHyperFLOW2D

       Copyright (C)  1995-2013 by Serge A. Suchkov
               Copyright policy: LGPL V3

CFD code for simulation 2D (flat/axisymmetrical) transient viscous
compressible multicomponent sub/trans/supersonic reacting gas flows.

include:
	- BC of I, II and III type for any of the dependent variables;
	- no-slip BC on walls;
	- multicomponent flow (4 base components);
	- temperature dependence properties of components;
	- chemical reactions with infinity speed (Zeldovich model);
	- turbulence (RANS models):
	    + Zero-equation models:
		* Prandtl;
		* van Driest;
		* Escudier;
		* Klebanoff;
	    + One-equation models:
		* Spalart-Allmaras model;
	    + Two-equations models:
		* Standard (Spalding) k-eps model;
		* Chien k-eps model;
		* Jones-Launder k-eps model;
		* Launder and Sharma, with Yap-correction k-eps model;
		* RNG k-eps model;
	- Tecplot (ascii) 2D output format;
	- Save&Restart simulation state;
	- Serial&Parallel(OpenMP, MPI) version;
                                                                                                                                                                                                       - Serial and Parallel (OpenMP,MPI) solver versions
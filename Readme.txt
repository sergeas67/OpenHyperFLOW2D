                      OpenHyperFLOW2D-CUDA
                   (Experimental CUDA branch)

Copyright (C)  1995-2014 by Serge A. Suchkov (sergeas67@gmail.com)
                Copyright policy: LGPL V3

CFD code for simulation 2D (flat/axisymmetrical) transient viscous
compressible multicomponent sub/trans/supersonic reacting gas flows.

Basic features
==============

	- BC of I, II and III type for any of the dependent variables;
	- no-slip BC on walls;
	- multicomponent flow (4 base components);
	- temperature dependence properties of components(-);
	- chemical reactions with infinity speed (Zeldovich model);
	- turbulence (RANS/URANS models):
	    + Zero-equation models:
		* Prandtl;
		* van Driest;
		* Escudier;
		* Klebanoff;
		* Smagorinsky-Lilly;
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

Experimental features
======================
1. Additional blending factor functions (BFF), very unstable !
2. Separate local blending factor (LBF) for each equation, may be unstable in some cases.
3. CUDA based computational core.

*Prebuilded binary executble compiled for sm 3.0

OpenHyperFLOW2D it is open source clone of HyperFLOW2D (in-house CFD code)
with reduced functionality, and can be used without restriction for educational
and academic purpose. Reference at the publication is obligatory.

----------------------------------------------------------------------------------------

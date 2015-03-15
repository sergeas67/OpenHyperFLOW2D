# Open-Source CFD simulation package **OpenHyperFLOW2D** #

|![http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/doc/OpenHyperFLOW2D-Logo.png](http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/doc/OpenHyperFLOW2D-Logo.png) | Copyright (C)  1995-2014 by Serge A. Suchkov |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------|
|  | Copyright policy: LGPL V3 |

**OpenHyperFLOW2D** it is parallel (C++/MPI/CUDA/OpenMP) research CFD code for simulation 2D (flat/axisymmetrical) transient viscous compressible multicomponent sub/trans/supersonic reacting gas flows.

## Include: ##
  * Own simple text-based preprocessor
  * Automatically generated cartesian mesh (does not need extra mesh)
  * Boundary conditions of I, II and III type for any of the dependent variables;
  * No-slip boundary conditions on walls;
  * Multicomponent flow (4 base components);
  * Temperature dependence properties of components;
  * Chemical reactions with infinity speed (_Zeldovich_ model);
  * Turbulence (RANS/URANS models):
    * Zero-equation models:
      * _Prandtl_;
      * _van Driest_;
      * _Escudier_;
      * _Klebanoff_;
      * _Smagorinsky-Lilly_;
    * One-equation models:
      * _Spalart-Allmaras_ model;
    * Two-equations models:
      * Standard (_Spalding_) k-eps model;
      * _Chien_ k-eps model;
      * _Jones-Launder_ k-eps model;
      * _Launder_ and _Sharma_, with _Yap_-correction k-eps model;
      * RNG k-eps model;

  * Tecplot (ascii) 2D output format;
  * Save&Restart simulation state;
  * Serial&Parallel(OpenMP, MPI, CUDA<sup><code>*</code></sup>) version;


## Test Cases: ##

|  **Test Case**  | **Source**  |  **Input data file**  |
|:----------------|:------------|:----------------------|
| Oblique Shock | [![](http://openhyperflow2d.googlecode.com/svn/branches/OpenHyperFLOW2D-CUDA/doc/ObliqueShock.png)](http://www.aero.polimi.it/freecase/?OpenFOAM_%2B_Code_Aster:Aerodynamic_problems:Oblique_shock) | http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/TestCases/ObliqueShock.dat |
| 15 deg. Wedge | <img src='http://openhyperflow2d.googlecode.com/svn/branches/OpenHyperFLOW2D-CUDA/doc/Wedge.png' width='600'> <table><thead><th> <a href='http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/TestCases/Wedge.dat'>http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/TestCases/Wedge.dat</a> </th></thead><tbody>
<tr><td> Mach 3 supersonic flow over a forward-facing step </td><td><a href='http://www.youtube.com/watch?feature=player_embedded&v=S7BRdFXaG04' target='_blank'><img src='http://img.youtube.com/vi/S7BRdFXaG04/0.jpg' width='425' height=344 /></a> </td><td><a href='http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/TestCases/Step.dat'>http://openhyperflow2d.googlecode.com/svn/trunk/OpenHyperFLOW2D/TestCases/Step.dat</a> </td></tr></tbody></table>

<b>OpenHyperFLOW2D</b> it is open source clone of <b>HyperFLOW2D</b> (in-house CFD code) with reduced functionality, and can be used without restriction for educational and academic purpose. Reference at the publication is obligatory.<br>
<br>
<br>
<hr><br>
<br>
<br>
<blockquote><sup><code>*</code></sup> - Experimental <b><i>OpenHyperFLOW2D-CUDA</i></b> branch.
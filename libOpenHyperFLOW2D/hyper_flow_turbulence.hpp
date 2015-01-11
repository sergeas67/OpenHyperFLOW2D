/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 11/01/2015                                                    *
*******************************************************************************/
#ifndef _hyper_flow_turbulence_hpp_
#define _hyper_flow_turbulence_hpp_

#define  NUM_COMPONENTS 3 // 0,1,2,3 ... NUM_COMPONENTS - number of _ADDITIONAL_ (!) components,
//  e.g. if you have 3 components, NUM_COMPONENTS=2 (1(base)+2(additional)=3).

#define    i2d_k       4+NUM_COMPONENTS   // k - index (k-eps, k-omega)
#define    i2d_nu_t    4+NUM_COMPONENTS   // eddy visc. index (SA)
#define    i2d_eps     4+NUM_COMPONENTS+1 // eps index (k-eps)
#define    i2d_omega   4+NUM_COMPONENTS+1 // omega index (k-omega, k-omega SST)

enum     TurbulenceCondType2D {
    TCT_No_Turbulence_2D            = 0x0,        // 0
// k-eps BC                         
    TCT_k_CONST_2D                  = 0x01,       // 1  k
    TCT_eps_CONST_2D                = 0x02,       // 2  eps
                                    
    TCT_dkdx_NULL_2D                = 0x04,       // 3
    TCT_depsdx_NULL_2D              = 0x08,       // 4
                                    
    TCT_dkdy_NULL_2D                = 0x010,      // 5
    TCT_depsdy_NULL_2D              = 0x020,      // 6

    TCT_d2kdx2_NULL_2D              = 0x040,      // 7
    TCT_d2epsdx2_NULL_2D            = 0x080,      // 8
                                    
    TCT_d2kdy2_NULL_2D              = 0x0100,     // 9
    TCT_d2epsdy2_NULL_2D            = 0x0200,     // 10
// k-eps model
    TCT_k_eps_Model_2D              = 0x0400,     // 11 - k-eps model
// Zero Equation models    
    TCT_Prandtl_Model_2D            = 0x0800,     // 12  Zero equation models (Prandtl types and it's modifications)
// Integral Model f(Re)    
    TCT_Integral_Model_2D           = 0x01000,    // 13
// Special BC on Wall for k-eps model
    TCT_eps_mud2kdx2_WALL_2D        = 0x02000,    // 14
    TCT_eps_mud2kdy2_WALL_2D        = 0x04000,    // 15
    TCT_eps_Cmk2kXn_WALL_2D         = 0x08000,    // 16
// SA model    
    TCT_Spalart_Allmaras_Model_2D   = 0x010000,   // 17
// k-omega model    
    TCT_k_omega_Model_2D            = 0x020000,   // 18
// k-omega SST model    
    TCT_k_omega_SST_Model_2D        = 0x040000,   // 19
// Baldwin-Lomax model    
    TCT_Baldwin_Lomax_Model_2D      = 0x080000,   // 20
// nut_92 (Sekundov) model          
    TCT_nut_92_Model_2D             = 0x0100000,  // 21
// Smagorinsky-Lilly model          
    TCT_Smagorinsky_Model_2D        = 0x0200000,  // 22
};                                                 

enum TurbulenceExtendedModel {
     TEM_Prandtl,         // 0 - zero equation Prandtl model
     TEM_vanDriest,       // 1 - zero equation van Driest
     TEM_Escudier,        // 2 - zero equation Escudier
     TEM_Klebanoff,       // 3 - zero equation Klebanoff
     TEM_k_eps_Std,       // 4 - Standart k-eps model
     TEM_k_eps_Chien,     // 5 - Chien k-eps model
     TEM_k_eps_JL,        // 6 - Jones-Launder k-eps model
     TEM_k_eps_LSY,       // 7 - Launder and Sharma, with Yap-correction
     TEM_k_eps_RNG,       // 8 - RNG k-eps model
     TEM_k_eps_Realisable,// 9 - Realisable k-eps model  (TODO)
     TEM_Spalart_Allmaras,// 10 - SA model               
     TEM_Baldwin_Lomax,   // 11 - BL model               (TODO)
     TEM_nut_92_Sekundov, // 12 - nut_92 model           (TODO)
     TEM_k_omega_Wilcox,  // 14 - k-omega model          (TODO)
     TEM_k_omega_SST,     // 15 - k-omega SST model      (TODO)
     TEM_Smagorinsky,     // 16 - Smagorinsky-Lilly model
};

// omega as eps alias 
#define TCT_omega_CONST_2D     TCT_eps_CONST_2D
#define TCT_domegadx_NULL_2D   TCT_depsdx_NULL_2D
#define TCT_domegady_NULL_2D   TCT_depsdy_NULL_2D
#define TCT_d2omegadx2_NULL_2D TCT_d2epsdx2_NULL_2D
#define TCT_d2omegady2_NULL_2D TCT_d2epsdy2_NULL_2D
// eddy viscosity as k alias
#define TCT_nu_t_CONST_2D      TCT_k_CONST_2D
#define TCT_dnu_t_dx_NULL_2D   TCT_dkdx_NULL_2D
#define TCT_dnu_t_dy_NULL_2D   TCT_dkdy_NULL_2D

// Macro turbulence condition types as combination of TurbulenceCond bit flags
enum TurbulenceNodeType2D {
     TNT_UNDEF_2D = 0,
     TNT_FC_2D    = TCT_k_CONST_2D | TCT_eps_CONST_2D,
     TNT_D0X_2D   = TCT_dkdx_NULL_2D | TCT_depsdx_NULL_2D,
     TNT_D0Y_2D   = TCT_dkdy_NULL_2D | TCT_depsdy_NULL_2D
};

template <class T, int a>
struct FlowNodeCore2D;

template <class T, int a>
class FlowNode2D;

template <class T, int a>
class FlowNodeTurbulence2D {

    friend struct FlowNodeCore2D<T,a>;
    friend class  FlowNode2D<T,a>;
    T          I;                    // turbulence intensity (%)
public:
    ulong      TurbType;             // 0x0   - no turbulence, 
                                     // 0x01  - Integral model, 
                                     // 0x02  - zero equations model, 
                                     // 0x04  - k-eps  model
    T          l_min;                // Minimal lenght from node to wall
    T          y_plus;               // y+=l_min * U_t/mu  
    T          Re_local;             // Re local.
    T          mu_t,lam_t;           // turbulence viscosity and turbulence heat conductivity
    T          dkdx,dkdy,            // dk/dx, dk/dy
               depsdx,depsdy;        // deps/dx, deps/dy
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
    inline ulong  isTurbulenceCond2D(ulong);
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
    inline void   SetTurbulenceIntensity2D(T i=0.005);
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
    inline void   SetTurbulenceCond2D(ulong tct) {
        TurbType = TurbType | tct;
    }
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
    inline void   CleanTurbulenceCond2D(ulong tct) {
        TurbType = (TurbType^tct)&TurbType;
    }
};

//template <class T, int a>
//T FlowNodeTurbulence2D<T,a>::I=0.005;  // turbulence intensity (%)

template <class T, int a>
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
inline ulong FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(ulong tct) {
    return ((TurbType & tct) == tct);  
}
template <class T, int a>
inline void FlowNodeTurbulence2D<T,a>::SetTurbulenceIntensity2D(T i) {
    FlowNodeTurbulence2D<T,a>::I = i;
}
#endif //_hyper_flow_turbulence_hpp_ 

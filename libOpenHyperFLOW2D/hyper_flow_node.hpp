/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   last update: 01/02/2015                                                    *
*******************************************************************************/

#ifndef  _hyper_flow_node_hpp
#define  _hyper_flow_node_hpp

#include "libFlow/flow2d.hpp"

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif

#include "libOpenHyperFLOW2D/hyper_flow_turbulence.hpp"

inline FP  Sgn(FP a)
{
    if ( a != 0. ) return fabs(a)/a;
    else           return 0;
}

#define  NUM_COMPONENTS 3 // 0,1,2,3 ... NUM_COMPONENTS - number of _ADDITIONAL_ (!) components,
//  e.g. if you have 3 components, NUM_COMPONENTS=2 (1(base)+2(additional)=3).

#ifdef _CUDA_
#define NUM_EQ 6+NUM_COMPONENTS
#endif // _CUDA_

#define    h_fu      0
#define    h_ox      1
#define    h_cp      2
#define    h_air     3

enum  SolverMode {
      SM_EULER = 0,
      SM_NS
};

enum     FlowType {
         FT_FLAT,         // 2D-Flat flow
         FT_AXISYMMETRIC, // 2D-Axisymmetric flow
};

// <------------- 2D --------------->

#define    i2d_Rho      0
#define    i2d_RhoU     1
#define    i2d_RhoV     2
#define    i2d_RhoE     3

#define    i2d_Yfu     4+NUM_COMPONENTS-3
#define    i2d_Yox     4+NUM_COMPONENTS-2
#define    i2d_Ycp     4+NUM_COMPONENTS-1

// 2D Condition types
enum CondType2D   {                  // Condition type (bit flags)
    CT_NO_COND_2D      = 0x0,        // 0  - no BC
    CT_Rho_CONST_2D    = 0x01,       // 1  - Rho constant (Dirichlet BC)
    CT_U_CONST_2D      = 0x02,       // 2  - RhoU constant (Dirichlet BC)
    CT_V_CONST_2D      = 0x04,       // 3  - RhoV constant (Dirichlet BC)
    CT_T_CONST_2D      = 0x08,       // 4  - RhoE constant (Dirichlet BC)
    CT_Y_CONST_2D      = 0x010,      // 5  - RhoY1...RhoYn constant (Dirichlet BC)
    CT_dRhodx_NULL_2D  = 0x020,      // 6  - dRho/dx = 0 (Neumann BC)
    CT_dUdx_NULL_2D    = 0x040,      // 7  - dRhoU/dx = 0 (Neumann BC)
    CT_dVdx_NULL_2D    = 0x080,      // 8  - dRhoV/dx = 0 (Neumann BC)
    CT_dTdx_NULL_2D    = 0x0100,     // 9  - dRhoE/dx = 0 (Neumann BC)
    CT_dYdx_NULL_2D    = 0x0200,     // 10 - dRhoY1...RhoYn/dx = 0 (Neumann BC)
    CT_dRhody_NULL_2D  = 0x0400,     // 11 - dRho/dy = 0 (Neumann BC)
    CT_dUdy_NULL_2D    = 0x0800,     // 12 - dRhoU/dy = 0 (Neumann BC)
    CT_dVdy_NULL_2D    = 0x01000,    // 13 - dRhoV/dy = 0 (Neumann BC)
    CT_dTdy_NULL_2D    = 0x02000,    // 14 - dRhoE/dy = 0 (Neumann BC)
    CT_dYdy_NULL_2D    = 0x04000,    // 15 - dRhoY1...RhoYn/dy = 0 (Neumann BC)
    CT_d2Rhodx2_NULL_2D= 0x08000,    // 16 - d2Rho/dx2 = 0 (Cauchy BC)
    CT_d2Udx2_NULL_2D  = 0x010000,   // 17 - d2RhoU/dx2 = 0 (Cauchy BC)
    CT_d2Vdx2_NULL_2D  = 0x020000,   // 18 - d2RhoV/dx2 = 0 (Cauchy BC)
    CT_d2Tdx2_NULL_2D  = 0x040000,   // 19 - d2RhoE/dx2 = 0 (Cauchy BC)
    CT_d2Ydx2_NULL_2D  = 0x080000,   // 20 - d2RhoY1...RhoYn/dx2 = 0 (Cauchy BC)
    CT_d2Rhody2_NULL_2D= 0x0100000,  // 21 - d2Rho/dy2 = 0 (Cauchy BC)
    CT_d2Udy2_NULL_2D  = 0x0200000,  // 22 - d2RhoU/dy2 = 0 (Cauchy BC)
    CT_d2Vdy2_NULL_2D  = 0x0400000,  // 23 - d2RhoV/dy2 = 0 (Cauchy BC)
    CT_d2Tdy2_NULL_2D  = 0x0800000,  // 24 - d2RhoE/dy2 = 0 (Cauchy BC)
    CT_d2Ydy2_NULL_2D  = 0x01000000, // 25 - d2RhoY1...RhoYn/dy2 = 0 (Cauchy BC)
    CT_TIME_DEPEND_2D  = 0x02000000, // 26 - time depended Dirichlet BC / F=f(t) /
    CT_WALL_NO_SLIP_2D = 0x04000000, // 27 - (if this bit set this node is WALL whis no-slip condition)
    CT_WALL_LAW_2D     = 0x08000000, // 28 - (if this bit set this node is WALL law  condition)
    CT_GAS_2D          = 0x010000000,// 29 - (if this bit set this node is GAS)
    CT_BL_REFINEMENT_2D= 0x020000000,// 30 - (if this bit set this node content additional array of boundary layer nodes)
    CT_SOLID_2D        = 0x040000000,// 31 - (if this bit set this node is SOLID)
    CT_NODE_IS_SET_2D  = 0x080000000,// 32 - (if this bit set, node alredy initialized, use for seek uninitialized nodes)
    CT_LIQUID_2D       = 0x0100000000,//33 - node is LIQUID
    CT_NONREFLECTED_2D = 0x0200000000,//34 - non-reflected BC
    CT_4TH_ORDER_X_2D  = 0x0400000000,//35 - 4-th order X bound
    CT_4TH_ORDER_Y_2D  = 0x0800000000,//36 - 4-th order Y bound
};


// Macro bound types as combination of CondType bit flags
enum     NodeType2D {
    NT_UNDEF_2D      = 0,                                            // Undefinded node type
    NT_FC_2D         = CT_Rho_CONST_2D | CT_U_CONST_2D | CT_V_CONST_2D | CT_Y_CONST_2D | CT_T_CONST_2D  |
                       CT_NODE_IS_SET_2D,                            // External gas flow
    NT_D0X_2D        = CT_NODE_IS_SET_2D | CT_dRhodx_NULL_2D | CT_dUdx_NULL_2D | CT_dVdx_NULL_2D | CT_dTdx_NULL_2D  |
                       CT_dYdx_NULL_2D,                              // Nongradient condition  (dF/dx = 0) in x direction
    NT_D2X_2D        = CT_NODE_IS_SET_2D | CT_d2Rhodx2_NULL_2D | CT_d2Udx2_NULL_2D | CT_d2Vdx2_NULL_2D | CT_d2Tdx2_NULL_2D  |
                       CT_d2Ydx2_NULL_2D,                            // Soft condition  (d2F/dx2 = 0) in x direction

    NT_D0Y_2D        = CT_NODE_IS_SET_2D | CT_dRhody_NULL_2D | CT_dUdy_NULL_2D | CT_dVdy_NULL_2D   | CT_dTdy_NULL_2D |
                       CT_dYdy_NULL_2D,                              // Nongradient condition  (dF/dy = 0) in y direction

    NT_D2Y_2D        = CT_NODE_IS_SET_2D | CT_d2Rhody2_NULL_2D | CT_d2Udy2_NULL_2D | CT_d2Vdy2_NULL_2D | CT_d2Tdy2_NULL_2D |
                       CT_d2Ydy2_NULL_2D,                            // Soft condition  (d2F/dy2 = 0) in y direction

    NT_AY_2D         = CT_NODE_IS_SET_2D | NT_D0X_2D | CT_U_CONST_2D,// AXIAL condition (axi has Y direction)
    NT_AX_2D         = CT_NODE_IS_SET_2D | NT_D0Y_2D | CT_V_CONST_2D,// AXIAL condition (axi has X direction)
    NT_WALL_LAW_2D   = CT_NODE_IS_SET_2D | CT_WALL_LAW_2D,           // Wall slip condition
    NT_WNS_2D        = CT_NODE_IS_SET_2D | CT_WALL_NO_SLIP_2D |      // Wall no-slip
                       CT_U_CONST_2D | CT_V_CONST_2D,                // condition
    NT_S_2D          = CT_SOLID_2D   | CT_NODE_IS_SET_2D,            // Solid body
    NT_F_2D          = !CT_SOLID_2D  | CT_NODE_IS_SET_2D,            // Internal gas area
    NT_FC_TIME_DEPEND_2D = CT_Rho_CONST_2D | CT_U_CONST_2D | CT_V_CONST_2D | CT_Y_CONST_2D | CT_T_CONST_2D  |
                      CT_TIME_DEPEND_2D | CT_NODE_IS_SET_2D,         // Time depend External gas flow
    NT_FARFIELD_2D   = NT_FC_2D | CT_NONREFLECTED_2D                 // Nonreflected far-field BC
};


template <class T, int a>
struct FlowNodeCore2D {
    T S[6+a];   // Ro, Ro*U, Ro*V, Ro*e,k*Ro, eps*Ro, Ro*Y1...Ro*.Yk;
    T dSdx[6+a];// dRo/dx, d(Ro*U)/dx,d(Ro*V)/dx, d(Ro*e)/dx, d(Ro*k)/dx, d(Ro*eps)/dx, d(Ro*Y1)/dx...d(Ro*.Yk)/dx;
    T dSdy[6+a];// dRo/dy, d(Ro*U)/dy,d(Ro*V)/dy, d(Ro*e)/dy, d(Ro*k)/dy, d(Ro*eps)/dy, d(Ro*Y1)/dy...d(Ro*.Yk)/dy;
};

// 2D-FlowNode class
template <class T, int a>
class FlowNode2D: public  FlowNodeCore2D<T,a>,
                  public  FlowNodeTurbulence2D<T,a> {

public:
    static const int NumEq;        // Number of equations in system (6 + num components - 1)
    static FlowType  FT;           // Flow type (flat or axisymmetric)
    static T         Hu[a+1];      // Hu[a]  Specific heat of formation
    static int       isSrcAdd;     // is additional Src[] present ?
#ifdef    _UNIFORM_MESH_
    static
#endif // _UNIFORM_MESH_
    T                dx,dy;       // dx, dy(dr)
    T                x,y;         // x, y(r)
    int              ix,iy;       // i_x,i_y
// Neighboring nodes
    FlowNode2D<T,a>*            UpNode;
    FlowNode2D<T,a>*            DownNode;
    FlowNode2D<T,a>*            LeftNode;
    FlowNode2D<T,a>*            RightNode;
    
    T          p;                    // pressure 
    int        idXl;                 // is left node present ? (0 or 1)
    int        idYu;                 // is up node present ? (0 or 1)
    int        idXr;                 // is right node present ? (0 or 1)
    int        idYd;                 // is down node present ? (0 or 1)
    int        NGX;                  // dF/dx integer coeff. 0, 1 or -1
    int        NGY;                  // dF/dy integer coeff. 0, 1 or -1
    //int        NodeID;             // Material ID (for solid Nodes)
    ulong      CT;                   // Condition type  (bit flags combination)
    int        i_wall,j_wall;        // neast wall coordinates
    T          beta[6+a];            // superlocal blending factor (SLBF).
    T          Q_conv;               // Convective heat flux 
    T          time;                 // Global time.
    T          k,R,lam,mu,CP,Diff;   // Cp/Cv, gas constant (R/m), lam, mu, Cp, mu/Cp.
    T          Tf;                   // ignition temperature
    T          A[6+a];               // A
    T          B[6+a];               // B
    T          F[6+a];               // F (only for axisymmetric flow)
    T          RX[6+a];              // RX
    T          RY[6+a];              // RY
    T          Src[6+a];             // Sources members
    T          SrcAdd[6+a];          // Additional virtual sources on the wall
    T          Tg,U,V,               // Temperature, U and V velocity components of gas
               Y[a+1];               // Y[a+1] for gas ...
    T          Uw,Vw;                // U and V wall velocity components
    T          droYdx[a+1];          // d(Ro*Y)/dx
    T          droYdy[a+1];          // d(Ro*Y)/dy
    T          dUdx,dUdy,            // dU/dx,dU/dy
               dVdx,dVdy,            // dV/dx,dV/dy
               dTdx,dTdy;            // dT/dx,dT/dy
    T          BGX;                  // dF/dx angle coeff. (cos(Ayz))
    T          BGY;                  // dF/dy angle coeff. (cos(Axz))
    FlowNode2D(FlowType ft = FT_FLAT); // flat model - default for 2D
    FlowNode2D(T RO,
             T U,
             T V,
             T P,
             T K,
             T RR,
             T LAM,
             T MU,
             T* Y,
             FlowType    FT= FT_FLAT,
             CondType2D  cT= CT_NODE_IS_SET_2D,
             T GX =(T)(1.),
             T GY =(T)(1.));
    FlowNode2D(FlowNode2D<T,a>&);
   ~FlowNode2D() {
        ;
    }

    inline void   CopyFlowNodeCore2D(FlowNodeCore2D<T,a>& fnc);
#ifdef _CUDA_
__host__ __device__
   inline void FillNode2D(int is_mu_t, 
                          int is_init, 
                          T sig_w, 
                          T sig_f, 
                          TurbulenceExtendedModel tem, 
                          T delta,
                          T _dx, T _dy,
                          T* _Hu,
                          int _isSrcAdd, 
                          SolverMode sm,
                          FlowType _FT);
__host__ __device__
  inline void  TurbModRANS2D(int is_mu_t,
                             int is_init,
                             TurbulenceExtendedModel tem,
                             T delta,
                             T _dx, T _dy,
                             FlowType _FT,
                             T _I);
__host__ __device__
    inline ulong  isCond2D(ulong);
#else    
    inline void   FillNode2D(int is_mu_t=0, int i=0, T sig_w=0.0, T sig_f=0.0, 
                             TurbulenceExtendedModel tem = TEM_k_eps_Std,
                             T delta = 0., SolverMode sm = SM_NS);
    
    inline void TurbModRANS2D(int is_mu_t, int is_init, TurbulenceExtendedModel tem, T delta = 0.);
    inline ulong  isCond2D(ulong);
#endif //_CUDA_    
    inline FlowNode2D<T,a>& operator = (FlowNode2D<T,a>& fn);
    inline FlowNode2D<T,a>& operator = (FlowNodeCore2D<T,a>& fc);
    FlowNode2D<T,a>& operator = (Flow& f);
    FlowNode2D<T,a>& operator = (Flow2D& f);
    inline FlowNode2D<T,a>& operator *= (T& m);
    inline Flow2D*  MakeFlow2D(Flow2D*);

    inline FlowType GetFlowType() {
        return FT;
    }
    inline void SetFlowType(FlowType ft) {
        FT=ft;
    }
    inline void SetCond2D(ulong ct) {
        CT = CT | ct ;
    }
    inline void CleanCond2D(ulong ct) {
        CT = (CT^ct)&CT;
    }
    // Turbulence models functions
#ifdef _CUDA_
 __host__ __device__
     void TurbulenceAxisymmAddOn(int is_init, FlowType _FT);
#else
     void TurbulenceAxisymmAddOn(int is_init);
#endif //_CUDA_ 
    int  CollectSubNodes2D();
    int    ExpandSubNodes2D();
};

template <class T, int a>
const int FlowNode2D<T,a>::NumEq=6+a; // Number of equations in 2D system (6 + num components - 1)

#ifndef _CUDA_
template <class T, int a>
FlowType FlowNode2D<T,a>::FT=FT_FLAT; // Flow type (flat or axisymmetric)
#endif //_CUDA_ 

template <class T, int a>
T  FlowNode2D<T,a>::Hu[a+1];          // Hu[a]  Specific heat of formation

#ifndef _CUDA_
template <class T, int a>
int FlowNode2D<T,a>::isSrcAdd=0;
#endif //_CUDA_ 

#ifdef    _UNIFORM_MESH_

template <class T, int a>
T  FlowNode2D<T,a>::dx;               // dx

template <class T, int a>
T  FlowNode2D<T,a>::dy;               // dy

#endif // _UNIFORM_MESH_

template <class T, int a>
#ifdef _CUDA_
__host__ __device__
#endif // _CUDA_
inline ulong FlowNode2D<T,a>::isCond2D(ulong ct) {
    return ((CT & ct) == ct);
}

template <class T, int a>
FlowNode2D<T,a>::FlowNode2D(FlowType ft) {
    memset(this,0,sizeof(FlowNode2D<T,a>));
    Y[a]  = 1.;
    FT    = ft;
    BGX   = 1.;
    BGY   = 1.;
    NGX   = 1;
    NGY   = 1;
    idXl  = 1;
    idYu  = 1;
    idXr  = 1;
    idYd  = 1;
    CT    = CT_NODE_IS_SET_2D;
}

template <class T, int a>
FlowNode2D<T,a>::FlowNode2D(T  RO,
                            T  U,
                            T  V,
                            T  P,
                            T  K,
                            T  RR,
                            T  LAM,
                            T  MU,
                            T* y,
                            FlowType    _FT,
                            CondType2D   cT,
                            T  GX,
                            T  GY) {
#ifdef __ICC
    __declspec(align(_ALIGN)) int      i;
    __declspec(align(_ALIGN)) T        Tmp1,Tmp3=0.;
#else
    int      i __attribute__ ((aligned (_ALIGN)));
    T        Tmp1 __attribute__ ((aligned (_ALIGN)));
    T        Tmp3 __attribute__ ((aligned (_ALIGN))) =0.;
#endif //__ICC

    FlowNodeCore2D<T,a>::S[i2d_Rho]  = RO;
    FlowNodeCore2D<T,a>::S[i2d_RhoU] = RO*U;
    FlowNodeCore2D<T,a>::S[i2d_RhoV] = RO*V;
    FT   = _FT;
    CT   =  cT | CT_NODE_IS_SET_2D;
    R    =  RR;
    lam  =  LAM;
    mu   =  MU;
    k    =  K;
    CP   =  k*R/(k-1);
    p    =  P;

    Diff   =  (lam+FlowNodeTurbulence2D<T,a>::lam_t)/CP;
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Y[i] = y[i];
        FlowNodeCore2D<T,a>::S[i+4] = Y[i]*FlowNodeCore2D<T,a>::S[i2d_Rho];
    }

    Tmp1 = FlowNodeCore2D<T,a>::S[i2d_Rho];
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Tmp3+=Hu[i]*FlowNodeCore2D<T,a>::S[i+4];
        Tmp1-= FlowNodeCore2D<T,a>::S[i+4];
    }

    Tmp3+=Hu[a]*Tmp1;

    FlowNodeCore2D<T,a>::S[i2d_RhoE] = P/(k-1.)+RO*(U*U+V*V)*0.5+Tmp3;

    BGX =1.;
    BGY =1.;

    NGX =1;
    NGY =1;

    idXl=1;
    idYu=1;
    idXr=1;
    idYd=1;
    //FillNode2D(0,1);
}

template <class T, int a>
FlowNode2D<T,a>::FlowNode2D(FlowNode2D<T,a>& fn) {
    memcpy(this,&fn,sizeof(FlowNode2D<T,a>));
}

#ifdef _CUDA_
template <class T, int a>
__host__ __device__
inline void FlowNode2D<T,a>::FillNode2D(int is_mu_t, 
                                        int is_init, 
                                        T sig_w, 
                                        T sig_f, 
                                        TurbulenceExtendedModel tem, 
                                        T delta,
                                        T _dx, T _dy,
                                        T* _Hu,
                                        int _isSrcAdd,
                                        SolverMode sm,
                                        FlowType _FT) {

    
    if(isCond2D(CT_SOLID_2D)) return;

    if(FlowNodeCore2D<T,a>::S[i2d_Rho]==0)
        return;
    else if(k < 1)
        return;
#ifdef __ICC
    __declspec(align(_ALIGN)) unsigned int i;
    __declspec(align(_ALIGN)) T  sxx,txy,syy,qx,qy;   // viscous stresses & diffusion fluxes
    __declspec(align(_ALIGN)) T L, _mu, _lam, t00;    //, G;
    __declspec(align(_ALIGN)) T Tmp1,Tmp2,Tmp3=0.;
#else
    unsigned int i __attribute__ ((aligned (_ALIGN)));
     // viscous stresses & diffusion fluxes
    T  sxx   __attribute__ ((aligned (_ALIGN)));
    T  txy   __attribute__ ((aligned (_ALIGN)));
    T  syy   __attribute__ ((aligned (_ALIGN)));
    T  qx    __attribute__ ((aligned (_ALIGN)));
    T  qy    __attribute__ ((aligned (_ALIGN)));
    T  L     __attribute__ ((aligned (_ALIGN)));
    T  _mu   __attribute__ ((aligned (_ALIGN)));
    T  _lam  __attribute__ ((aligned (_ALIGN)));
    T  t00   __attribute__ ((aligned (_ALIGN)));
   // T  G     __attribute__ ((aligned (_ALIGN)));
    T  Tmp1  __attribute__ ((aligned (_ALIGN)));
    T  Tmp2  __attribute__ ((aligned (_ALIGN)));
    T  Tmp3  __attribute__ ((aligned (_ALIGN))) =0.;
#endif //__ICC

    k = CP/(CP-R);
    
    if(isCond2D(CT_U_CONST_2D))
        FlowNodeCore2D<T,a>::S[i2d_RhoU] = U*FlowNodeCore2D<T,a>::S[i2d_Rho];
    else
        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];

    if(isCond2D(CT_V_CONST_2D))
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = V*FlowNodeCore2D<T,a>::S[i2d_Rho];
    else
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];
    
    if (sm == SM_NS) {
        
        _mu = _lam = 0.;
        
        if(is_init)
           FlowNodeTurbulence2D<T,a>::mu_t = FlowNodeTurbulence2D<T,a>::lam_t = 0.;

        if(FlowNodeTurbulence2D<T,a>::TurbType > 0)
           TurbModRANS2D(is_mu_t,is_init,tem,delta,_dx,_dy,_FT,0.005);
    }

    Tmp1   = FlowNodeCore2D<T,a>::S[i2d_Rho];

    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Tmp3+= _Hu[i]*FlowNodeCore2D<T,a>::S[i+4];
        Tmp1-= FlowNodeCore2D<T,a>::S[i+4];
    }

    Tmp3+=_Hu[a]*Tmp1;

    if(isCond2D(CT_WALL_LAW_2D)) {
        // Simple  model (TODO)
        Tmp1 = sqrt(U*U+V*V+1.e-30);

        FlowNodeCore2D<T,a>::S[i2d_RhoU] = Tmp1*BGX; 
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = Tmp1*BGY;

        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        // MKT-Model (removed)
    } else if(isCond2D(CT_WALL_NO_SLIP_2D)) {

        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];

        if( _isSrcAdd ) {
            SrcAdd[i2d_Rho]  = BGX*(U-Uw)*FlowNodeCore2D<T,a>::S[i2d_Rho]/_dx +
                               BGY*(V-Vw)*FlowNodeCore2D<T,a>::S[i2d_Rho]/_dy;
            SrcAdd[i2d_RhoU] = BGX*(U-Uw)*FlowNodeCore2D<T,a>::S[i2d_Rho];
            SrcAdd[i2d_RhoV] = BGY*(V-Vw)*FlowNodeCore2D<T,a>::S[i2d_Rho];
            SrcAdd[i2d_RhoE] = 0.;
        } else {
            SrcAdd[i2d_RhoU] = SrcAdd[i2d_RhoV] = SrcAdd[i2d_RhoE] = 0.;
        }

        U = Uw; // Gas move
        V = Vw; // together with wall
        
        
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(i=0;i<a;i++) {
            if(_isSrcAdd) {
                SrcAdd[4+i] = SrcAdd[i2d_Rho]*Y[i];
            } else {
                SrcAdd[4+i] = 0.;
            }
        }
        FlowNodeCore2D<T,a>::S[i2d_RhoU] = U * FlowNodeCore2D<T,a>::S[i2d_Rho];
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = V * FlowNodeCore2D<T,a>::S[i2d_Rho];
    } else {
        
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(int i=0;i<NUM_EQ;i++)
            SrcAdd[i] = 0.;
    }
    
    p  = (FlowNodeCore2D<T,a>::S[i2d_RhoE]-FlowNodeCore2D<T,a>::S[i2d_Rho]*(U*U+V*V)*0.5-Tmp3)*(k-1.0);

    Tg = p/R/FlowNodeCore2D<T,a>::S[i2d_Rho];

    if (sm == SM_NS) {
        
        FlowNodeTurbulence2D<T,a>::lam_t  = FlowNodeTurbulence2D<T,a>::mu_t*CP;

        if(is_mu_t) {
            if(isCond2D(CT_WALL_NO_SLIP_2D) || isCond2D(CT_WALL_LAW_2D)){
                _mu    = max(0,(mu+FlowNodeTurbulence2D<T,a>::mu_t*sig_w));
                _lam   = max(0,(lam+FlowNodeTurbulence2D<T,a>::lam_t*sig_w)); 
            } else {
                _mu  =  max(0,(mu+FlowNodeTurbulence2D<T,a>::mu_t*sig_f));
                _lam =  max(0,(lam+FlowNodeTurbulence2D<T,a>::lam_t*sig_f)); 
            }
        } else {
            _mu = mu;
            _lam = lam;
        }

        Diff   = _lam/CP;
        L=(2./3.)*_mu;                  // 2-nd viscosity (dilatation)

        if(_FT == FT_AXISYMMETRIC)
           Tmp2 = L*(dUdx+dVdy+_FT*V/y); // L*dilatation (2D)
        else
           Tmp2 = L*(dUdx+dVdy);        // L*dilatation (2D)
    }

    A[i2d_Rho]  = FlowNodeCore2D<T,a>::S[i2d_RhoU];
    A[i2d_RhoU] = p + FlowNodeCore2D<T,a>::S[i2d_RhoU]*U;
    A[i2d_RhoV] = FlowNodeCore2D<T,a>::S[i2d_RhoV]*U;
    A[i2d_RhoE] = (FlowNodeCore2D<T,a>::S[i2d_RhoE]+p)*U;

    B[i2d_Rho]  = FlowNodeCore2D<T,a>::S[i2d_RhoV];
    B[i2d_RhoU] = A[i2d_RhoV];
    B[i2d_RhoV] = p + FlowNodeCore2D<T,a>::S[i2d_RhoV]*V;
    B[i2d_RhoE] = (FlowNodeCore2D<T,a>::S[i2d_RhoE]+p)*V;
    
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i = 4; i < 4+a; i++) {
        B[i]=FlowNodeCore2D<T,a>::S[i]*V;
        A[i]=FlowNodeCore2D<T,a>::S[i]*U;
    }

    if( _FT == 1 ) {
        F[i2d_Rho]  = _FT*B[i2d_Rho];  
        F[i2d_RhoU] = _FT*A[i2d_RhoV];
        F[i2d_RhoV] = _FT*F[i2d_Rho]*V;
        F[i2d_RhoE] = _FT*B[i2d_RhoE]; 
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(i = 4; i < 4+a;i++)
            F[i]=_FT*B[i]; 
    }
    
    if (sm == SM_NS) {
        
        sxx   = 2.*_mu*dUdx - Tmp2; 
        syy   = 2.*_mu*dVdy - Tmp2;

        txy   = _mu*(dUdy+dVdx);

        qx    = _lam*dTdx;
        qy    = _lam*dTdy;

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
        for(i = 0; i < a+1; i++) {
            qx += Diff*(CP*Tg+_Hu[i])*droYdx[i];
            qy += Diff*(CP*Tg+_Hu[i])*droYdy[i];
        }


        RX[i2d_Rho]  = 0.;
        RX[i2d_RhoU] = sxx;
        RX[i2d_RhoV] = txy;
        RX[i2d_RhoE] = U*sxx+V*txy+qx;

        RY[i2d_Rho]  = 0;
        RY[i2d_RhoU] = txy;
        RY[i2d_RhoV] = syy;
        RY[i2d_RhoE] = U*txy+V*syy+qy;

        A[i2d_RhoU] = A[i2d_RhoU]-RX[i2d_RhoU];
        A[i2d_RhoV] = A[i2d_RhoV]-RX[i2d_RhoV];
        A[i2d_RhoE] = A[i2d_RhoE]-RX[i2d_RhoE];

        B[i2d_RhoU] = B[i2d_RhoU]-RY[i2d_RhoU];
        B[i2d_RhoV] = B[i2d_RhoV]-RY[i2d_RhoV];
        B[i2d_RhoE] = B[i2d_RhoE]-RY[i2d_RhoE];

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
        for(i = 4; i < 4+a; i++) {

            RX[i]=Diff*droYdx[i-4];
            RY[i]=Diff*droYdy[i-4];

            A[i]=A[i]-RX[i];
            B[i]=B[i]-RY[i];
        }

        if( _FT == FT_AXISYMMETRIC ) {
            t00   = 2*_mu*V/y - Tmp2;
            F[i2d_RhoU] -= RY[i2d_RhoU];
            F[i2d_RhoV] -= RY[i2d_RhoV]+t00;
            F[i2d_RhoE] -= RY[i2d_RhoE];

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
            for(i = 4; i < 4+a;i++)
                F[i]-=RY[i];
        } else {

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
            for(int i=0;i<NUM_EQ;i++)
                F[i]=0.;
        }
    }
}

template <class T, int a>
#ifdef _CUDA_
__host__ __device__
void  FlowNode2D<T,a>::TurbModRANS2D(int is_mu_t,
                                     int is_init,
                                     TurbulenceExtendedModel tem,
                                     T delta,
                                     T _dx, T _dy,
                                     FlowType _FT,
                                     T _I) {

 #ifdef __ICC
    __declspec(align(_ALIGN)) T l = max(_dy,max((FlowNodeTurbulence2D<T,a>::l_min),_dx)) * 0.41;
 #else
   T l  __attribute__ ((aligned (_ALIGN))) = max(_dy,max((FlowNodeTurbulence2D<T,a>::l_min),_dx)) * 0.41;
 #endif // __ICC
    if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Prandtl_Model_2D)) { 
 #ifdef __ICC
        __declspec(align(_ALIGN)) T n_0, A_p=26.;
 #else
       T n_0  __attribute__ ((aligned (_ALIGN)));
       T A_p  __attribute__ ((aligned (_ALIGN))) = 26.;
 #endif // __ICC
        n_0 = FlowNodeTurbulence2D<T,a>::l_min * 0.41;  // Basic Prandtl equation

        if(tem == TEM_Prandtl) {
            l = n_0;
        } else if (tem == TEM_vanDriest) {
            l = n_0 * (1.-exp(-FlowNodeTurbulence2D<T,a>::y_plus/A_p));
        } else if (tem == TEM_Escudier) {
            if(delta > 0.)
               l = min(n_0,0.09*delta);
            else
               l = n_0;
        } else if (tem == TEM_Klebanoff) {
            if(delta > 0.)
               l = n_0/sqrt(1+5.5*pow(FlowNodeTurbulence2D<T,a>::l_min/delta,6));
            else
               l = n_0;
        }

        FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));
        FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;

    } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_eps_Model_2D)) { 
#ifdef __ICC
        __declspec(align(_ALIGN)) T C1eps=1.44,C2eps=1.92,C_mu=0.09,sig_k=1.0,sig_eps=1.3;
        __declspec(align(_ALIGN)) T f1=1.,f2=1.,f_mu = 1.;                   // dumping function
        __declspec(align(_ALIGN)) T Rt, nu;                                  // turbulent Re
        __declspec(align(_ALIGN)) T G,L_k,L_eps,Mt=0;
#else
        T C1eps __attribute__ ((aligned (_ALIGN))) =1.44;
        T C2eps __attribute__ ((aligned (_ALIGN))) =1.92;
        T C_mu  __attribute__ ((aligned (_ALIGN))) =0.09;
        T sig_k __attribute__ ((aligned (_ALIGN))) =1.0;
        T sig_eps __attribute__ ((aligned (_ALIGN))) =1.3;
        T f1    __attribute__ ((aligned (_ALIGN))) =1.;
        T f2 __attribute__ ((aligned (_ALIGN))) =1.;
        T f_mu __attribute__ ((aligned (_ALIGN))) =1.;                   // dumping function
        T Rt __attribute__ ((aligned (_ALIGN)));
        T nu __attribute__ ((aligned (_ALIGN)));                         // turbulent Re
        T G __attribute__ ((aligned (_ALIGN)));
        T L_k __attribute__ ((aligned (_ALIGN)));
        T L_eps __attribute__ ((aligned (_ALIGN)));
        T Mt __attribute__ ((aligned (_ALIGN))) =0;
#endif // __ICC
        L_k   = L_eps = 0; 

        sig_k   = 1;
        sig_eps = 1.3;

#ifdef __ICC
        __declspec(align(_ALIGN)) T Tmp1    = dUdy + dVdx;
        __declspec(align(_ALIGN)) T Tmp2    = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        __declspec(align(_ALIGN)) T Tmp3    = dUdx*dUdx + dVdy*dVdy;
#else
        T Tmp1 __attribute__ ((aligned (_ALIGN)))   = dUdy + dVdx;
        T Tmp2 __attribute__ ((aligned (_ALIGN)))   = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        T Tmp3 __attribute__ ((aligned (_ALIGN)))   = dUdx*dUdx + dVdy*dVdy;
#endif // __ICC

        if(_FT)
          Tmp3  += U/y;

        if(FlowNodeTurbulence2D<T,a>::mu_t == 0)
           FlowNodeTurbulence2D<T,a>::mu_t = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));

        G     = FlowNodeTurbulence2D<T,a>::mu_t*(Tmp1*Tmp1+2*Tmp3);

        // Turbulent Reynolds number
        if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.0 && mu != 0.0)
          Rt =  FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]/mu; 
        else
          Rt = 0;

        if(tem == TEM_k_eps_Std) {
            // Std k-eps model
            // Calc various dumbing functions here
             f_mu  = 1.0;
             C1eps = 1.44;
             C2eps = 1.92;
             C_mu  = 0.09;
             sig_k = 1.0;
             sig_eps = 1.3;
         } else if(tem == TEM_k_eps_Chien) {
             // Chien k-eps model
             // Calc various dumbing functions here
             C1eps   = 1.35; 
             C2eps   = 1.8; 
             C_mu    = 0.09; 
             sig_k   = 1.0;  
             sig_eps = 1.3;
             f2      = 1.0 - 0.4/1.8* exp(-(Rt*Rt)/36.0);
             f_mu    = 1.0 - exp(-0.0115*FlowNodeTurbulence2D<T,a>::y_plus);
             L_k     = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_k]/(Tmp2*Tmp2);
             L_eps   = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_eps]/(Tmp2*Tmp2)*exp(-FlowNodeTurbulence2D<T,a>::y_plus/2.0);
             Mt      = 1.5 * FlowNodeCore2D<T,a>::S[i2d_k]/k/p; // Compressibility Correction  
         } else if(tem == TEM_k_eps_JL) { 
             // Jones-Launder k-eps model (JL)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-2.5/(1.0+Rt/50));
         } else if(tem == TEM_k_eps_LSY) {
            // Launder and Sharma, with Yap-correction (LSY)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-3.4/(1.0+Rt/50.0)/(1.0+Rt/50.0));
         }else if(tem == TEM_k_eps_RNG) {
            // RNG
             // Calc various dumbing functions here
             T nu_0 = 4.38;

             if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.)
               nu = sqrt(G)*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps];
             else
               nu = 0;

             C_mu  = 0.0845;
             C1eps = 1.42;
             C2eps = 1.68 + C_mu*nu*nu*nu*(1.0-nu/nu_0)/(1.+ 0.012*nu*nu*nu);
             sig_k = 0.7194;
             sig_eps = 0.7194;
             f_mu = 1.;
         }
      // eps calc from Prandtl model
      // l = l_min * 0.41
      // mu_t= S[i2d_Rho]*l*l*(dUdy + dVdx);
      // mu_t= 0.09 * k * k/eps
      // eps = 0.09 * k * k / (S[i2d_Rho] * l * l * (dUdy + dVdx))  
        if(is_init) {
#ifdef __ICC        
            __declspec(align(_ALIGN)) T TmpI = _I*sqrt(U*U+V*V+1.e-30);
#else
           T TmpI __attribute__ ((aligned (_ALIGN))) = _I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC           
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow((double)C_mu,3./4.)*pow((double)(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho]),3./2.)/l;
           if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0)
              FlowNodeTurbulence2D<T,a>::mu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);
        }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T TmpI = _I*sqrt(U*U+V*V+1.e-30);
#else
            T TmpI __attribute__ ((aligned (_ALIGN))) = _I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D)) {
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow((double)C_mu,3./4.)*pow((double)(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho]),3./2.)/l;
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_Cmk2kXn_WALL_2D) ) 
           FlowNodeCore2D<T,a>::S[i2d_eps] = pow((double)C_mu,3./4.)*pow((double)(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho]),3./2.)/l;

        if(is_mu_t && FlowNodeCore2D<T,a>::S[i2d_eps] != 0) {
            T nu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);   
            FlowNodeTurbulence2D<T,a>::mu_t = min(nu_t,(FlowNodeTurbulence2D<T,a>::mu_t));
        }
        if (!is_init) {

            A[i2d_k]  = FlowNodeCore2D<T,a>::S[i2d_k]  *U;
            A[i2d_eps]= FlowNodeCore2D<T,a>::S[i2d_eps]*U;

            B[i2d_k]   = FlowNodeCore2D<T,a>::S[i2d_k]  *V;
            B[i2d_eps] = FlowNodeCore2D<T,a>::S[i2d_eps]*V;

            RX[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdx); // +FlowNodeTurbulence2D<T,a>::dkdy
            RX[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdx); //+FlowNodeTurbulence2D<T,a>::depsdy

            RY[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdy); //+FlowNodeTurbulence2D<T,a>::dkdx
            RY[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdy); //+FlowNodeTurbulence2D<T,a>::depsdx

            A[i2d_k]   = A[i2d_k]   - RX[i2d_k]  ;
            A[i2d_eps] = A[i2d_eps] - RX[i2d_eps];  

            B[i2d_k]   = B[i2d_k]   - RY[i2d_k]  ;
            B[i2d_eps]  =B[i2d_eps] - RY[i2d_eps];

            if(FlowNodeCore2D<T,a>::S[i2d_k] != 0.) {


              SrcAdd[i2d_k] = SrcAdd[i2d_eps] = 0.;

              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D))
                  FlowNode2D<T,a>::Src[i2d_k] = (G  - FlowNodeCore2D<T,a>::S[i2d_eps]*(1+Mt)
                                                + L_k*FlowNodeCore2D<T,a>::S[i2d_Rho]);
              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D) && FlowNodeCore2D<T,a>::S[i2d_k] != 0)
                  FlowNode2D<T,a>::Src[i2d_eps] = (C1eps * f1 * FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k] * G -
                                                   C2eps * f2 * FlowNodeCore2D<T,a>::S[i2d_eps]*FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k]
                                                 + L_eps * FlowNodeCore2D<T,a>::S[i2d_Rho]);
            }
           //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init, _FT);
        }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T Cb1, Cb2, sig, _k, Cw1, Cw2, Cw3, Cv1, /*Ct1,*/ Ct2, /*Ct3,*/ Ct4, nu, ksi, fv1, fv2, ft2, fw, g, r, Omega, S_hat;
            __declspec(align(_ALIGN)) T Wxy, Div_nu;
#else
           T Cb1 __attribute__ ((aligned (_ALIGN)));
           T Cb2 __attribute__ ((aligned (_ALIGN)));
           T sig __attribute__ ((aligned (_ALIGN)));
           T _k __attribute__ ((aligned (_ALIGN)));
           T Cw1 __attribute__ ((aligned (_ALIGN)));
           T Cw2 __attribute__ ((aligned (_ALIGN)));
           T Cw3 __attribute__ ((aligned (_ALIGN)));
           T Cv1 __attribute__ ((aligned (_ALIGN))); /*Ct1,*/
           T Ct2 __attribute__ ((aligned (_ALIGN)));
           //T Ct3 __attribute__ ((aligned (_ALIGN)));
           T Ct4 __attribute__ ((aligned (_ALIGN)));
           T nu __attribute__ ((aligned (_ALIGN)));
           T ksi __attribute__ ((aligned (_ALIGN)));
           T fv1 __attribute__ ((aligned (_ALIGN)));
           T fv2 __attribute__ ((aligned (_ALIGN)));
           T ft2 __attribute__ ((aligned (_ALIGN)));
           T fw __attribute__ ((aligned (_ALIGN)));
           T g __attribute__ ((aligned (_ALIGN)));
           T r __attribute__ ((aligned (_ALIGN)));
           T Omega __attribute__ ((aligned (_ALIGN)));
           T S_hat __attribute__ ((aligned (_ALIGN)));
           T Wxy __attribute__ ((aligned (_ALIGN)));
           T Div_nu __attribute__ ((aligned (_ALIGN)));
#endif //__ICC
            fv1 = 1.0;

            if(is_init) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]/100.0; // initial turbulence visc.
            } else if (isCond2D(CT_WALL_NO_SLIP_2D) || isCond2D(CT_WALL_LAW_2D) || 
                       FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = 0.; 
            } else if (isCond2D(NT_FC_2D) ) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]*5.0;
            } else { // if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) 
                Cb1 = 0.1355;
                Cb2 = 0.622;
                sig = 2.0/3.0;
                _k  = 0.41;
                Cw1 = Cb1/(_k*_k) + (1 + Cb2)/sig;
                Cw2 = 0.3;
                Cw3 = 2.0;
                Cv1 = 7.1;
                Ct2 = 2.0;
                //Ct3 = 1.2; // orig 1.1
                Ct4 = 0.5;   // orig 2
                nu  = mu/FlowNodeCore2D<T,a>::S[i2d_Rho];
                ksi = FlowNodeCore2D<T,a>::S[i2d_nu_t]/nu;
                fv1 = ksi*ksi*ksi/(ksi*ksi*ksi + Cv1*Cv1*Cv1);
                fv2 = 1.0 - ksi/(1.0 + ksi * fv1);
                Wxy = 0.5*(dVdx - dUdy);
                Omega = sqrt(2.0*Wxy*Wxy);
                S_hat = Omega + FlowNodeCore2D<T,a>::S[i2d_nu_t]/(_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min) * fv2;

                if(S_hat < 0.3 * Omega)
                   S_hat = 0.3 * Omega;

                r   = min((FlowNodeCore2D<T,a>::S[i2d_nu_t]/(S_hat*_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min)),10.0);
                g   = r + Cw2*(pow((double)r,6.0) - r);

                fw  = g * pow((double)(1.0 + pow((double)Cw3,6.0))/(pow((double)g,6.0)+pow((double)Cw3,6.0)),1.0/6.0);
                ft2 = Ct2 * exp(-Ct4*ksi*ksi);

                A[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*U;
                B[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*V;

                Div_nu        =  (FlowNodeTurbulence2D<T,a>::dkdx + FlowNodeTurbulence2D<T,a>::dkdy );

                RX[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdx)/sig;
                RY[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdy)/sig;

                A[i2d_nu_t]   = A[i2d_nu_t]  - RX[i2d_nu_t];
                B[i2d_nu_t]   = B[i2d_nu_t]  - RY[i2d_nu_t];

                Src[i2d_nu_t] = Cb1*(1.0 - ft2)*S_hat*FlowNodeCore2D<T,a>::S[i2d_nu_t] - (Cw1*fw - Cb1/(_k*_k)*ft2)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min) + (Cb2*Div_nu*Div_nu)/sig;      
            }
            //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init, _FT);

            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho] * FlowNodeCore2D<T,a>::S[i2d_nu_t] * fv1;
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D)) {

        } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Integral_Model_2D) && mu != 0.0) {
            FlowNodeTurbulence2D<T,a>::Re_local = FlowNodeTurbulence2D<T,a>::l_min*sqrt(FlowNodeCore2D<T,a>::S[i2d_RhoU]*FlowNodeCore2D<T,a>::S[i2d_RhoU]+
                                                                                        FlowNodeCore2D<T,a>::S[i2d_RhoV]*FlowNodeCore2D<T,a>::S[i2d_RhoV]+
                                                                                        1.e-30)/mu; 
        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D) && !is_init) {

        } else if ( FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Smagorinsky_Model_2D)) {
#ifdef __ICC
          __declspec(align(_ALIGN)) T Omega;
          __declspec(align(_ALIGN)) T Wxy;
          __declspec(align(_ALIGN)) T Cs = 0.1; // Smagorinsky const. Need change for DSGS
#else
          T Omega __attribute__ ((aligned (_ALIGN)));
          T Wxy   __attribute__ ((aligned (_ALIGN)));
          T Cs    __attribute__ ((aligned (_ALIGN))) = 0.1; // Smagorinsky const. Need change for DSGS 
#endif // __ICC
#ifdef    _UNIFORM_MESH_
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( _dx*_dy);
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( _dx*_dy);
 #endif // __ICC
#else
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( _dx * _dy );
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( _dx * _dy );    
 #endif // __ICC
#endif // _UNIFORM_MESH_
          Wxy = 0.5*(dVdx - dUdy);
          Omega = sqrt(2.0*Wxy*Wxy);
            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho] * (Cs * _delta) * (Cs * _delta) * Omega;
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }
        }
}
#else
template <class T, int a>
void  FlowNode2D<T,a>::TurbModRANS2D(int is_mu_t,
                                     int is_init,
                                     TurbulenceExtendedModel tem,
                                     T delta ) {

 #ifdef __ICC
    __declspec(align(_ALIGN)) T l = max(dy,max((FlowNodeTurbulence2D<T,a>::l_min),dx)) * 0.41;
 #else
   T l  __attribute__ ((aligned (_ALIGN))) = max(dy,max((FlowNodeTurbulence2D<T,a>::l_min),dx)) * 0.41;
 #endif // __ICC
    if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Prandtl_Model_2D)) { 
 #ifdef __ICC
        __declspec(align(_ALIGN)) T n_0, A_p=26.;
 #else
       T n_0  __attribute__ ((aligned (_ALIGN)));
       T A_p  __attribute__ ((aligned (_ALIGN))) = 26.;
 #endif // __ICC
        n_0 = FlowNodeTurbulence2D<T,a>::l_min * 0.41;  // Basic Prandtl equation

        if(tem == TEM_Prandtl) {
            l = n_0;
        } else if (tem == TEM_vanDriest) {
            l = n_0 * (1.-exp(-FlowNodeTurbulence2D<T,a>::y_plus/A_p));
        } else if (tem == TEM_Escudier) {
            if(delta > 0.)
               l = min(n_0,0.09*delta);
            else
               l = n_0;
        } else if (tem == TEM_Klebanoff) {
            if(delta > 0.)
               l = n_0/sqrt(1+5.5*pow(FlowNodeTurbulence2D<T,a>::l_min/delta,6));
            else
               l = n_0;
        }

        FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));
        FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;

    } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_eps_Model_2D)) { 
#ifdef __ICC
        __declspec(align(_ALIGN)) T C1eps=1.44,C2eps=1.92,C_mu=0.09,sig_k=1.0,sig_eps=1.3;
        __declspec(align(_ALIGN)) T f1=1.,f2=1.,f_mu = 1.;                   // dumping function
        __declspec(align(_ALIGN)) T Rt, nu;                                  // turbulent Re
        __declspec(align(_ALIGN)) T G,L_k,L_eps,Mt=0;
#else
        T C1eps __attribute__ ((aligned (_ALIGN))) =1.44;
        T C2eps __attribute__ ((aligned (_ALIGN))) =1.92;
        T C_mu  __attribute__ ((aligned (_ALIGN))) =0.09;
        T sig_k __attribute__ ((aligned (_ALIGN))) =1.0;
        T sig_eps __attribute__ ((aligned (_ALIGN))) =1.3;
        T f1    __attribute__ ((aligned (_ALIGN))) =1.;
        T f2 __attribute__ ((aligned (_ALIGN))) =1.;
        T f_mu __attribute__ ((aligned (_ALIGN))) =1.;                   // dumping function
        T Rt __attribute__ ((aligned (_ALIGN)));
        T nu __attribute__ ((aligned (_ALIGN)));                         // turbulent Re
        T G __attribute__ ((aligned (_ALIGN)));
        T L_k __attribute__ ((aligned (_ALIGN)));
        T L_eps __attribute__ ((aligned (_ALIGN)));
        T Mt __attribute__ ((aligned (_ALIGN))) =0;
#endif // __ICC
        L_k   = L_eps = 0; 

        sig_k   = 1;
        sig_eps = 1.3;

#ifdef __ICC
        __declspec(align(_ALIGN)) T Tmp1    = dUdy + dVdx;
        __declspec(align(_ALIGN)) T Tmp2    = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        __declspec(align(_ALIGN)) T Tmp3    = dUdx*dUdx + dVdy*dVdy;
#else
        T Tmp1 __attribute__ ((aligned (_ALIGN)))   = dUdy + dVdx;
        T Tmp2 __attribute__ ((aligned (_ALIGN)))   = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        T Tmp3 __attribute__ ((aligned (_ALIGN)))   = dUdx*dUdx + dVdy*dVdy;
#endif // __ICC

        if(FT)
          Tmp3  += U/y;

        if(FlowNodeTurbulence2D<T,a>::mu_t == 0)
           FlowNodeTurbulence2D<T,a>::mu_t = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));

        G     = FlowNodeTurbulence2D<T,a>::mu_t*(Tmp1*Tmp1+2*Tmp3);

        // Turbulent Reynolds number
        if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.0 && mu != 0.0)
          Rt =  FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]/mu; 
        else
          Rt = 0;

        if(tem == TEM_k_eps_Std) {
            // Std k-eps model
            // Calc various dumbing functions here
             f_mu  = 1.0;
             C1eps = 1.44;
             C2eps = 1.92;
             C_mu  = 0.09;
             sig_k = 1.0;
             sig_eps = 1.3;
         } else if(tem == TEM_k_eps_Chien) {
             // Chien k-eps model
             // Calc various dumbing functions here
             C1eps   = 1.35; 
             C2eps   = 1.8; 
             C_mu    = 0.09; 
             sig_k   = 1.0;  
             sig_eps = 1.3;
             f2      = 1.0 - 0.4/1.8* exp(-(Rt*Rt)/36.0);
             f_mu    = 1.0 - exp(-0.0115*FlowNodeTurbulence2D<T,a>::y_plus);
             L_k     = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_k]/(Tmp2*Tmp2);
             L_eps   = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_eps]/(Tmp2*Tmp2)*exp(-FlowNodeTurbulence2D<T,a>::y_plus/2.0);
             Mt      = 1.5 * FlowNodeCore2D<T,a>::S[i2d_k]/k/p; // Compressibility Correction  
         } else if(tem == TEM_k_eps_JL) { 
             // Jones-Launder k-eps model (JL)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-2.5/(1.0+Rt/50));
         } else if(tem == TEM_k_eps_LSY) {
            // Launder and Sharma, with Yap-correction (LSY)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-3.4/(1.0+Rt/50.0)/(1.0+Rt/50.0));
         }else if(tem == TEM_k_eps_RNG) {
            // RNG
             // Calc various dumbing functions here
             T nu_0 = 4.38;

             if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.)
               nu = sqrt(G)*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps];
             else
               nu = 0;

             C_mu  = 0.0845;
             C1eps = 1.42;
             C2eps = 1.68 + C_mu*nu*nu*nu*(1.0-nu/nu_0)/(1.+ 0.012*nu*nu*nu);
             sig_k = 0.7194;
             sig_eps = 0.7194;
             f_mu = 1.;
         }
      // eps calc from Prandtl model
      // l = l_min * 0.41
      // mu_t= S[i2d_Rho]*l*l*(dUdy + dVdx);
      // mu_t= 0.09 * k * k/eps
      // eps = 0.09 * k * k / (S[i2d_Rho] * l * l * (dUdy + dVdx))  
        if(is_init) {
#ifdef __ICC        
            __declspec(align(_ALIGN)) T TmpI = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#else
           T TmpI __attribute__ ((aligned (_ALIGN))) = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC           
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;
           if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0)
              FlowNodeTurbulence2D<T,a>::mu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);
        }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T TmpI = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#else
            T TmpI __attribute__ ((aligned (_ALIGN))) = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D)) {
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_Cmk2kXn_WALL_2D) ) 
           FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;

        if(is_mu_t && FlowNodeCore2D<T,a>::S[i2d_eps] != 0) {
            T nu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNoSolverModedeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);   
            FlowNodeTurbulence2D<T,a>::mu_t = min(nu_t,(FlowNodeTurbulence2D<T,a>::mu_t));
        }
        if (!is_init) {

            A[i2d_k]  = FlowNodeCore2D<T,a>::S[i2d_k]  *U;
            A[i2d_eps]= FlowNodeCore2D<T,a>::S[i2d_eps]*U;

            B[i2d_k]   = FlowNodeCore2D<T,a>::S[i2d_k]  *V;
            B[i2d_eps] = FlowNodeCore2D<T,a>::S[i2d_eps]*V;

            RX[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdx); // +FlowNodeTurbulence2D<T,a>::dkdy
            RX[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdx); //+FlowNodeTurbulence2D<T,a>::depsdy

            RY[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdy); //+FlowNodeTurbulence2D<T,a>::dkdx
            RY[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdy); //+FlowNodeTurbulence2D<T,a>::depsdx

            A[i2d_k]   = A[i2d_k]   - RX[i2d_k]  ;
            A[i2d_eps] = A[i2d_eps] - RX[i2d_eps];  

            B[i2d_k]   = B[i2d_k]   - RY[i2d_k]  ;
            B[i2d_eps]  =B[i2d_eps] - RY[i2d_eps];

            if(FlowNodeCore2D<T,a>::S[i2d_k] != 0.) {


              SrcAdd[i2d_k] = SrcAdd[i2d_eps] = 0.;

              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D))
                  FlowNode2D<T,a>::Src[i2d_k] = (G  - FlowNodeCore2D<T,a>::S[i2d_eps]*(1+Mt)
                                                + L_k*FlowNodeCore2D<T,a>::S[i2d_Rho]);
              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D) && FlowNodeCore2D<T,a>::S[i2d_k] != 0)
                  FlowNode2D<T,a>::Src[i2d_eps] = (C1eps * f1 * FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k] * G -
                                                   C2eps * f2 * FlowNodeCore2D<T,a>::S[i2d_eps]*FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k]
                                                 + L_eps * FlowNodeCore2D<T,a>::S[i2d_Rho]);
            }
           //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init);
        }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T Cb1, Cb2, sig, _k, Cw1, Cw2, Cw3, Cv1, /*Ct1,*/ Ct2, Ct3, Ct4, nu, ksi, fv1, fv2, ft2, fw, g, r, Omega, S_hat;
            __declspec(align(_ALIGN)) T Wxy, Div_nu;
#else
           T Cb1 __attribute__ ((aligned (_ALIGN)));
           T Cb2 __attribute__ ((aligned (_ALIGN)));
           T sig __attribute__ ((aligned (_ALIGN)));
           T _k __attribute__ ((aligned (_ALIGN)));
           T Cw1 __attribute__ ((aligned (_ALIGN)));
           T Cw2 __attribute__ ((aligned (_ALIGN)));
           T Cw3 __attribute__ ((aligned (_ALIGN)));
           T Cv1 __attribute__ ((aligned (_ALIGN))); /*Ct1,*/
           T Ct2 __attribute__ ((aligned (_ALIGN)));
           T Ct3 __attribute__ ((aligned (_ALIGN)));
           T Ct4 __attribute__ ((aligned (_ALIGN)));
           T nu __attribute__ ((aligned (_ALIGN)));
           T ksi __attribute__ ((aligned (_ALIGN)));
           T fv1 __attribute__ ((aligned (_ALIGN)));
           T fv2 __attribute__ ((aligned (_ALIGN)));
           T ft2 __attribute__ ((aligned (_ALIGN)));
           T fw __attribute__ ((aligned (_ALIGN)));
           T g __attribute__ ((aligned (_ALIGN)));
           T r __attribute__ ((aligned (_ALIGN)));
           T Omega __attribute__ ((aligned (_ALIGN)));
           T S_hat __attribute__ ((aligned (_ALIGN)));
           T Wxy __attribute__ ((aligned (_ALIGN)));
           T Div_nu __attribute__ ((aligned (_ALIGN)));
#endif //__ICC
            fv1 = 1.0;

            if(is_init) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]/100.0; // initial turbulence visc.
            } else if (isCond2D(CT_WALL_NO_SLIP_2D) || isCond2D(CT_WALL_LAW_2D) || 
                       FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = 0.; 
            } else if (isCond2D(NT_FC_2D) ) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]*5.0;
            } else { // if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) 
                Cb1 = 0.1355;
                Cb2 = 0.622;
                sig = 2.0/3.0;
                _k  = 0.41;
                Cw1 = Cb1/(_k*_k) + (1 + Cb2)/sig;
                Cw2 = 0.3;
                Cw3 = 2.0;
                Cv1 = 7.1;
                Ct2 = 2.0;
                Ct3 = 1.2; // orig 1.1
                Ct4 = 0.5; // orig 2
                nu  = mu/FlowNodeCore2D<T,a>::S[i2d_Rho];
                ksi = FlowNodeCore2D<T,a>::S[i2d_nu_t]/nu;
                fv1 = ksi*ksi*ksi/(ksi*ksi*ksi + Cv1*Cv1*Cv1);
                fv2 = 1.0 - ksi/(1.0 + ksi * fv1);
                Wxy = 0.5*(dVdx - dUdy);
                Omega = sqrt(2.0*Wxy*Wxy);
                S_hat = Omega + FlowNodeCore2D<T,a>::S[i2d_nu_t]/(_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min) * fv2;

                if(S_hat < 0.3 * Omega)
                   S_hat = 0.3 * Omega;

                r   = min((FlowNodeCore2D<T,a>::S[i2d_nu_t]/(S_hat*_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min)),10.0);
                g   = r + Cw2*(pow(r,6.0) - r);

                fw  = g * pow((1.0 + pow(Cw3,6.0))/(pow(g,6.0)+pow(Cw3,6.0)),1.0/6.0);
                ft2 = Ct2 * exp(-Ct4*ksi*ksi);

                A[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*U;
                B[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*V;

                Div_nu        =  (FlowNodeTurbulence2D<T,a>::dkdx + FlowNodeTurbulence2D<T,a>::dkdy );

                RX[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdx)/sig;
                RY[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdy)/sig;

                A[i2d_nu_t]   = A[i2d_nu_t]  - RX[i2d_nu_t];
                B[i2d_nu_t]   = B[i2d_nu_t]  - RY[i2d_nu_t];

                Src[i2d_nu_t] = Cb1*(1.0 - ft2)*S_hat*FlowNodeCore2D<T,a>::S[i2d_nu_t] - (Cw1*fw - Cb1/(_k*_k)*ft2)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min) + (Cb2*Div_nu*Div_nu)/sig;      
            }
            //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init);

            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho] * FlowNodeCore2D<T,a>::S[i2d_nu_t] * fv1;
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D)) {

        } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Integral_Model_2D) && mu != 0.0) {
            FlowNodeTurbulence2D<T,a>::Re_local = FlowNodeTurbulence2D<T,a>::l_min*sqrt(FlowNodeCore2D<T,a>::S[i2d_RhoU]*FlowNodeCore2D<T,a>::S[i2d_RhoU]+
                                                                                        FlowNodeCore2D<T,a>::S[i2d_RhoV]*FlowNodeCore2D<T,a>::S[i2d_RhoV]+
                                                                                        1.e-30)/mu; 
        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D) && !is_init) {

        } else if ( FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Smagorinsky_Model_2D)) {
#ifdef __ICC
          __declspec(align(_ALIGN)) T Omega;
          __declspec(align(_ALIGN)) T Wxy;
          __declspec(align(_ALIGN)) T Cs = 0.1; // Smagorinsky const. Need change for DSGS
#else
          T Omega __attribute__ ((aligned (_ALIGN)));
          T Wxy   __attribute__ ((aligned (_ALIGN)));
          T Cs    __attribute__ ((aligned (_ALIGN))) = 0.1; // Smagorinsky const. Need change for DSGS 
#endif // __ICC
#ifdef    _UNIFORM_MESH_
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( FlowNode2D<T,a>::dx*FlowNode2D<T,a>::dy);
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( FlowNode2D<T,a>::dx*FlowNode2D<T,a>::dy);
 #endif // __ICC
#else
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( dx * dy );
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( dx * dy );    
 #endif // __ICC
#endif // _UNIFORM_MESH_
          Wxy = 0.5*(dVdx - dUdy);
          Omega = sqrt(2.0*Wxy*Wxy);
            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho] * (Cs * _delta) * (Cs * _delta) * Omega;
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }
        }
}

#endif // _CUDA_




template <class T, int a>
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
void FlowNode2D<T,a>::TurbulenceAxisymmAddOn(int is_init, FlowType _FT) {
      if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_eps_Model_2D) && !is_init) {
        F[i2d_k]   = _FT*(mu+FlowNodeTurbulence2D<T,a>::mu_t)*FlowNodeTurbulence2D<T,a>::dkdy;
        F[i2d_eps] = _FT*(mu+FlowNodeTurbulence2D<T,a>::mu_t/1.3)*FlowNodeTurbulence2D<T,a>::depsdy;
      } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) && !is_init) {
        F[i2d_nu_t]   = _FT*(mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdy;
      } else if(is_init) {
        F[i2d_nu_t]   = F[i2d_eps] = 0.;
        Src[i2d_nu_t] = Src[i2d_eps] = 0.;
      }
}

#else
// General FlowNode2D class function
template <class T, int a>
inline void FlowNode2D<T,a>::FillNode2D(int is_mu_t, 
                                        int is_init, 
                                        T sig_w, 
                                        T sig_f, 
                                        TurbulenceExtendedModel tem, 
                                        T delta,
                                        SolverMode sm) {

    if(isCond2D(CT_SOLID_2D)) return;

    if(FlowNodeCore2D<T,a>::S[i2d_Rho]==0)
        return;
    else if(k < 1)
        return;
#ifdef __ICC
    __declspec(align(_ALIGN)) unsigned int i;
    __declspec(align(_ALIGN)) T  sxx,txy,syy,qx,qy;   // viscous stresses & diffusion fluxes
    __declspec(align(_ALIGN)) T L, _mu, _lam, t00, G;
    __declspec(align(_ALIGN)) T Tmp1,Tmp2,Tmp3=0.;
#else
    unsigned int i __attribute__ ((aligned (_ALIGN)));
     // viscous stresses & diffusion fluxes
    T  sxx   __attribute__ ((aligned (_ALIGN)));
    T  txy   __attribute__ ((aligned (_ALIGN)));
    T  syy   __attribute__ ((aligned (_ALIGN)));
    T  qx    __attribute__ ((aligned (_ALIGN)));
    T  qy    __attribute__ ((aligned (_ALIGN)));
    T  L     __attribute__ ((aligned (_ALIGN)));
    T  _mu   __attribute__ ((aligned (_ALIGN)));
    T  _lam  __attribute__ ((aligned (_ALIGN)));
    T  t00   __attribute__ ((aligned (_ALIGN)));
    T  G     __attribute__ ((aligned (_ALIGN)));
    T  Tmp1  __attribute__ ((aligned (_ALIGN)));
    T  Tmp2  __attribute__ ((aligned (_ALIGN)));
    T  Tmp3  __attribute__ ((aligned (_ALIGN))) =0.;
#endif //__ICC

    k = CP/(CP-R);
    

    if(isCond2D(CT_U_CONST_2D))
        FlowNodeCore2D<T,a>::S[i2d_RhoU] = U*FlowNodeCore2D<T,a>::S[i2d_Rho];
    else
        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];

    if(isCond2D(CT_V_CONST_2D))
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = V*FlowNodeCore2D<T,a>::S[i2d_Rho];
    else
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];
    
    Tmp1   = FlowNodeCore2D<T,a>::S[i2d_Rho];
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Tmp3+=Hu[i]*FlowNodeCore2D<T,a>::S[i+4];
        Tmp1-= FlowNodeCore2D<T,a>::S[i+4];
    }

    Tmp3+=Hu[a]*Tmp1;

    if(isCond2D(CT_WALL_LAW_2D)) {
        // Simple  model
        Tmp1 = sqrt(U*U+V*V+1.e-30);

        FlowNodeCore2D<T,a>::S[i2d_RhoU] = Tmp1*BGX; 
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = Tmp1*BGY;

        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        // MKT-Model (removed)
    } else if(isCond2D(CT_WALL_NO_SLIP_2D)) {

        U    = FlowNodeCore2D<T,a>::S[i2d_RhoU]/FlowNodeCore2D<T,a>::S[i2d_Rho];
        V    = FlowNodeCore2D<T,a>::S[i2d_RhoV]/FlowNodeCore2D<T,a>::S[i2d_Rho];

        if( isSrcAdd ) {
            SrcAdd[i2d_Rho]  = BGX*(U-Uw)*FlowNodeCore2D<T,a>::S[i2d_Rho]/dx +
                              BGY*(V-Vw)*FlowNodeCore2D<T,a>::S[i2d_Rho]/dy;
            SrcAdd[i2d_RhoU] = BGX*(U-Uw)*FlowNodeCore2D<T,a>::S[i2d_Rho];
            SrcAdd[i2d_RhoV] = BGY*(V-Vw)*FlowNodeCore2D<T,a>::S[i2d_Rho];
            SrcAdd[i2d_RhoE] = 0.;
        } else {
            SrcAdd[i2d_RhoU] = SrcAdd[i2d_RhoV] = SrcAdd[i2d_RhoE] = 0.;
        }

        U = Uw; // Gas move
        V = Vw; // together with wall
        
        
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(i=0;i<a;i++) {
            if(FlowNode2D<T,a>::isSrcAdd) {
                SrcAdd[4+i] = SrcAdd[i2d_Rho]*Y[i];
            } else {
                SrcAdd[4+i] = 0.;
            }
        }
        FlowNodeCore2D<T,a>::S[i2d_RhoU] = U * FlowNodeCore2D<T,a>::S[i2d_Rho];
        FlowNodeCore2D<T,a>::S[i2d_RhoV] = V * FlowNodeCore2D<T,a>::S[i2d_Rho];
    } else {
        
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(int i=0;i<NumEq;i++)
            SrcAdd[i] = 0.;
    }

    p  = (k-1.)*(FlowNodeCore2D<T,a>::S[i2d_RhoE]-FlowNodeCore2D<T,a>::S[i2d_Rho]*(U*U+V*V)*0.5-Tmp3);
    Tg = p/R/FlowNodeCore2D<T,a>::S[i2d_Rho];

    A[i2d_Rho]  = FlowNodeCore2D<T,a>::S[i2d_RhoU];
    A[i2d_RhoU] = p + FlowNodeCore2D<T,a>::S[i2d_RhoU]*U;
    A[i2d_RhoV] = FlowNodeCore2D<T,a>::S[i2d_RhoV]*U;
    A[i2d_RhoE] = (FlowNodeCore2D<T,a>::S[i2d_RhoE]+p)*U;

    B[i2d_Rho]  = FlowNodeCore2D<T,a>::S[i2d_RhoV];
    B[i2d_RhoU] = A[i2d_RhoV];
    B[i2d_RhoV] = p + FlowNodeCore2D<T,a>::S[i2d_RhoV]*V;
    B[i2d_RhoE] = (FlowNodeCore2D<T,a>::S[i2d_RhoE]+p)*V;
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i = 4; i < 4+a; i++) {
        B[i]=FlowNodeCore2D<T,a>::S[i]*V;
        A[i]=FlowNodeCore2D<T,a>::S[i]*U;
    }

    if( FT == FT_AXISYMMETRIC ) {
        F[i2d_Rho]  = FT*B[i2d_Rho];  
        F[i2d_RhoU] = FT*A[i2d_RhoV];
        F[i2d_RhoV] = FT*F[i2d_Rho]*V;
        F[i2d_RhoE] = FT*B[i2d_RhoE]; 
        
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
        for(i = 4; i < 4+a;i++)
            F[i]=FT*B[i]; 
    } else if( sm == SM_EULER) {
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
            for(int i=0;i<NumEq;i++)
                F[i]=0.;
    }

    if (sm == SM_NS) {

        _mu = _lam = G = 0.;
        if(is_init)
           FlowNodeTurbulence2D<T,a>::mu_t = FlowNodeTurbulence2D<T,a>::lam_t = 0.;

        TurbModRANS2D(is_mu_t,is_init,tem,delta);

        FlowNodeTurbulence2D<T,a>::lam_t  = FlowNodeTurbulence2D<T,a>::mu_t*CP;

        if(is_mu_t) {
            if(isCond2D(CT_WALL_NO_SLIP_2D) || isCond2D(CT_WALL_LAW_2D)){
                _mu    = max(0,(mu+FlowNodeTurbulence2D<T,a>::mu_t*sig_w));
                _lam   = max(0,(lam+FlowNodeTurbulence2D<T,a>::lam_t*sig_w)); 
            } else {
                _mu  =  max(0,(mu+FlowNodeTurbulence2D<T,a>::mu_t*sig_f));
                _lam =  max(0,(lam+FlowNodeTurbulence2D<T,a>::lam_t*sig_f)); 
            }
        } else {
            _mu = mu;
            _lam = lam;
        }

        Diff   = _lam/CP;
        L=(2./3.)*_mu;                  // 2-nd viscosity (dilatation)

        if(FT == FT_AXISYMMETRIC)
           Tmp2 = L*(dUdx+dVdy+FT*V/y); // L*dilatation (2D)
        else
           Tmp2 = L*(dUdx+dVdy);        // L*dilatation (2D)
        
        sxx   = 2.*_mu*dUdx - Tmp2; 
        syy   = 2.*_mu*dVdy - Tmp2;

        txy   = _mu*(dUdy+dVdx);

        qx    = _lam*dTdx;
        qy    = _lam*dTdy;

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
        for(i = 0; i < a+1; i++) {
            qx += Diff*(CP*Tg+Hu[i])*droYdx[i];
            qy += Diff*(CP*Tg+Hu[i])*droYdy[i];
        }

        RX[i2d_Rho]  = 0.;
        RX[i2d_RhoU] = sxx;
        RX[i2d_RhoV] = txy;
        RX[i2d_RhoE] = U*sxx+V*txy+qx;

        RY[i2d_Rho]  = 0;
        RY[i2d_RhoU] = txy;
        RY[i2d_RhoV] = syy;
        RY[i2d_RhoE] = U*txy+V*syy+qy;

        A[i2d_RhoU] = A[i2d_RhoU]-RX[i2d_RhoU];
        A[i2d_RhoV] = A[i2d_RhoV]-RX[i2d_RhoV];
        A[i2d_RhoE] = A[i2d_RhoE]-RX[i2d_RhoE];

        B[i2d_RhoU] = B[i2d_RhoU]-RY[i2d_RhoU];
        B[i2d_RhoV] = B[i2d_RhoV]-RY[i2d_RhoV];
        B[i2d_RhoE] = B[i2d_RhoE]-RY[i2d_RhoE];


    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
        for(i = 4; i < 4+a; i++) {

            RX[i]=Diff*droYdx[i-4];
            RY[i]=Diff*droYdy[i-4];

            A[i]=A[i]-RX[i];
            B[i]=B[i]-RY[i];
        }

        if( FT == FT_AXISYMMETRIC ) {
            t00   = 2*_mu*V/y - Tmp2;
            F[i2d_RhoU] -= RY[i2d_RhoU];
            F[i2d_RhoV] -= RY[i2d_RhoV]+t00;
            F[i2d_RhoE] -= RY[i2d_RhoE];

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
            for(i = 4; i < 4+a;i++)
                F[i]-=RY[i];
        } else {

    #ifdef _CUDA_
    #pragma unroll
    #endif // _CUDA_
            for(int i=0;i<NumEq;i++)
                F[i]=0.;
        }

    }
}

template <class T, int a>
void  FlowNode2D<T,a>::TurbModRANS2D(int is_mu_t,
                                     int is_init,
                                     TurbulenceExtendedModel tem,
                                     T delta ) {

 #ifdef __ICC
    __declspec(align(_ALIGN)) T l = max(dy,max((FlowNodeTurbulence2D<T,a>::l_min),dx)) * 0.41;
 #else
   T l  __attribute__ ((aligned (_ALIGN))) = max(dy,max((FlowNodeTurbulence2D<T,a>::l_min),dx)) * 0.41;
 #endif // __ICC
    if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Prandtl_Model_2D)) { 
 #ifdef __ICC
        __declspec(align(_ALIGN)) T n_0, A_p=26.;
 #else
       T n_0  __attribute__ ((aligned (_ALIGN)));
       T A_p  __attribute__ ((aligned (_ALIGN))) = 26.;
 #endif // __ICC
        n_0 = FlowNodeTurbulence2D<T,a>::l_min * 0.41;  // Basic Prandtl equation

        if(tem == TEM_Prandtl) {
            l = n_0;
        } else if (tem == TEM_vanDriest) {
            l = n_0 * (1.-exp(-FlowNodeTurbulence2D<T,a>::y_plus/A_p));
        } else if (tem == TEM_Escudier) {
            if(delta > 0.)
               l = min(n_0,0.09*delta);
            else
               l = n_0;
        } else if (tem == TEM_Klebanoff) {
            if(delta > 0.)
               l = n_0/sqrt(1+5.5*pow(FlowNodeTurbulence2D<T,a>::l_min/delta,6));
            else
               l = n_0;
        }

        FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));
        FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;

    } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_eps_Model_2D)) { 
#ifdef __ICC
        __declspec(align(_ALIGN)) T C1eps=1.44,C2eps=1.92,C_mu=0.09,sig_k=1.0,sig_eps=1.3;
        __declspec(align(_ALIGN)) T f1=1.,f2=1.,f_mu = 1.;                   // dumping function
        __declspec(align(_ALIGN)) T Rt, nu;                                  // turbulent Re
        __declspec(align(_ALIGN)) T G,L_k,L_eps,Mt=0;
#else
        T C1eps __attribute__ ((aligned (_ALIGN))) =1.44;
        T C2eps __attribute__ ((aligned (_ALIGN))) =1.92;
        T C_mu  __attribute__ ((aligned (_ALIGN))) =0.09;
        T sig_k __attribute__ ((aligned (_ALIGN))) =1.0;
        T sig_eps __attribute__ ((aligned (_ALIGN))) =1.3;
        T f1    __attribute__ ((aligned (_ALIGN))) =1.;
        T f2 __attribute__ ((aligned (_ALIGN))) =1.;
        T f_mu __attribute__ ((aligned (_ALIGN))) =1.;                   // dumping function
        T Rt __attribute__ ((aligned (_ALIGN)));
        T nu __attribute__ ((aligned (_ALIGN)));                         // turbulent Re
        T G __attribute__ ((aligned (_ALIGN)));
        T L_k __attribute__ ((aligned (_ALIGN)));
        T L_eps __attribute__ ((aligned (_ALIGN)));
        T Mt __attribute__ ((aligned (_ALIGN))) =0;
#endif // __ICC
        L_k   = L_eps = 0; 

        sig_k   = 1;
        sig_eps = 1.3;

#ifdef __ICC
        __declspec(align(_ALIGN)) T Tmp1    = dUdy + dVdx;
        __declspec(align(_ALIGN)) T Tmp2    = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        __declspec(align(_ALIGN)) T Tmp3    = dUdx*dUdx + dVdy*dVdy;
#else
        T Tmp1 __attribute__ ((aligned (_ALIGN)))   = dUdy + dVdx;
        T Tmp2 __attribute__ ((aligned (_ALIGN)))   = FlowNodeCore2D<T,a>::S[i2d_Rho]*l;
        T Tmp3 __attribute__ ((aligned (_ALIGN)))   = dUdx*dUdx + dVdy*dVdy;
#endif // __ICC

        if(FT)
          Tmp3  += U/y;

        if(FlowNodeTurbulence2D<T,a>::mu_t == 0)
           FlowNodeTurbulence2D<T,a>::mu_t = FlowNodeCore2D<T,a>::S[i2d_Rho]*l*l*max(fabs(dUdy),fabs(dVdx));

        G     = FlowNodeTurbulence2D<T,a>::mu_t*(Tmp1*Tmp1+2*Tmp3);

        // Turbulent Reynolds number
        if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.0 && mu != 0.0)
          Rt =  FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]/mu; 
        else
          Rt = 0;

        if(tem == TEM_k_eps_Std) {
            // Std k-eps model
            // Calc various dumbing functions here
             f_mu  = 1.0;
             C1eps = 1.44;
             C2eps = 1.92;
             C_mu  = 0.09;
             sig_k = 1.0;
             sig_eps = 1.3;
         } else if(tem == TEM_k_eps_Chien) {
             // Chien k-eps model
             // Calc various dumbing functions here
             C1eps   = 1.35; 
             C2eps   = 1.8; 
             C_mu    = 0.09; 
             sig_k   = 1.0;  
             sig_eps = 1.3;
             f2      = 1.0 - 0.4/1.8* exp(-(Rt*Rt)/36.0);
             f_mu    = 1.0 - exp(-0.0115*FlowNodeTurbulence2D<T,a>::y_plus);
             L_k     = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_k]/(Tmp2*Tmp2);
             L_eps   = -2.0 * mu * FlowNodeCore2D<T,a>::S[i2d_eps]/(Tmp2*Tmp2)*exp(-FlowNodeTurbulence2D<T,a>::y_plus/2.0);
             Mt      = 1.5 * FlowNodeCore2D<T,a>::S[i2d_k]/k/p; // Compressibility Correction  
         } else if(tem == TEM_k_eps_JL) { 
             // Jones-Launder k-eps model (JL)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-2.5/(1.0+Rt/50));
         } else if(tem == TEM_k_eps_LSY) {
            // Launder and Sharma, with Yap-correction (LSY)
             // Calc various dumbing functions here
             C1eps = 1.44; 
             C2eps = 1.92; 
             C_mu  = 0.09; 
             sig_k = 1.0;  
             sig_eps = 1.3;
             f_mu = exp(-3.4/(1.0+Rt/50.0)/(1.0+Rt/50.0));
         }else if(tem == TEM_k_eps_RNG) {
            // RNG
             // Calc various dumbing functions here
             T nu_0 = 4.38;

             if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0.)
               nu = sqrt(G)*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps];
             else
               nu = 0;

             C_mu  = 0.0845;
             C1eps = 1.42;
             C2eps = 1.68 + C_mu*nu*nu*nu*(1.0-nu/nu_0)/(1.+ 0.012*nu*nu*nu);
             sig_k = 0.7194;
             sig_eps = 0.7194;
             f_mu = 1.;
         }
      // eps calc from Prandtl model
      // l = l_min * 0.41
      // mu_t= S[i2d_Rho]*l*l*(dUdy + dVdx);
      // mu_t= 0.09 * k * k/eps
      // eps = 0.09 * k * k / (S[i2d_Rho] * l * l * (dUdy + dVdx))  
        if(is_init) {
#ifdef __ICC        
            __declspec(align(_ALIGN)) T TmpI = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#else
           T TmpI __attribute__ ((aligned (_ALIGN))) = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC           
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;
           if(FlowNodeCore2D<T,a>::S[i2d_eps] != 0)
              FlowNodeTurbulence2D<T,a>::mu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);
        }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T TmpI = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#else
            T TmpI __attribute__ ((aligned (_ALIGN))) = FlowNodeTurbulence2D<T,a>::I*sqrt(U*U+V*V+1.e-30);
#endif //__ICC
            FlowNodeCore2D<T,a>::S[i2d_k]   = 1.5*TmpI*TmpI*FlowNodeCore2D<T,a>::S[i2d_Rho];
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D)) {
            FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;
         }

        if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_Cmk2kXn_WALL_2D) ) 
           FlowNodeCore2D<T,a>::S[i2d_eps] = pow(C_mu,3./4.)*pow(FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_Rho],3./2.)/l;

        if(is_mu_t && FlowNodeCore2D<T,a>::S[i2d_eps] != 0) {
            T nu_t = fabs(C_mu * f_mu * FlowNodeCore2D<T,a>::S[i2d_k]*FlowNodeCore2D<T,a>::S[i2d_k]/FlowNodeCore2D<T,a>::S[i2d_eps]);   
            FlowNodeTurbulence2D<T,a>::mu_t = min(nu_t,(FlowNodeTurbulence2D<T,a>::mu_t));
        }
        if (!is_init) {

            A[i2d_k]  = FlowNodeCore2D<T,a>::S[i2d_k]  *U;
            A[i2d_eps]= FlowNodeCore2D<T,a>::S[i2d_eps]*U;

            B[i2d_k]   = FlowNodeCore2D<T,a>::S[i2d_k]  *V;
            B[i2d_eps] = FlowNodeCore2D<T,a>::S[i2d_eps]*V;

            RX[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdx); // +FlowNodeTurbulence2D<T,a>::dkdy
            RX[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdx); //+FlowNodeTurbulence2D<T,a>::depsdy

            RY[i2d_k]  = (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_k)*(FlowNodeTurbulence2D<T,a>::dkdy); //+FlowNodeTurbulence2D<T,a>::dkdx
            RY[i2d_eps]= (mu+FlowNodeTurbulence2D<T,a>::mu_t/sig_eps)*(FlowNodeTurbulence2D<T,a>::depsdy); //+FlowNodeTurbulence2D<T,a>::depsdx

            A[i2d_k]   = A[i2d_k]   - RX[i2d_k]  ;
            A[i2d_eps] = A[i2d_eps] - RX[i2d_eps];  

            B[i2d_k]   = B[i2d_k]   - RY[i2d_k]  ;
            B[i2d_eps]  =B[i2d_eps] - RY[i2d_eps];

            if(FlowNodeCore2D<T,a>::S[i2d_k] != 0.) {


              SrcAdd[i2d_k] = SrcAdd[i2d_eps] = 0.;

              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_CONST_2D))
                  FlowNode2D<T,a>::Src[i2d_k] = (G  - FlowNodeCore2D<T,a>::S[i2d_eps]*(1+Mt)
                                                + L_k*FlowNodeCore2D<T,a>::S[i2d_Rho]);
              if(!FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_eps_CONST_2D) && FlowNodeCore2D<T,a>::S[i2d_k] != 0)
                  FlowNode2D<T,a>::Src[i2d_eps] = (C1eps * f1 * FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k] * G -
                                                   C2eps * f2 * FlowNodeCore2D<T,a>::S[i2d_eps]*FlowNodeCore2D<T,a>::S[i2d_eps]/FlowNodeCore2D<T,a>::S[i2d_k]
                                                 + L_eps * FlowNodeCore2D<T,a>::S[i2d_Rho]);
            }
           //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init);
        }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
#ifdef __ICC
            __declspec(align(_ALIGN)) T Cb1, Cb2, sig, _k, Cw1, Cw2, Cw3, Cv1, /*Ct1,*/ Ct2, Ct3, Ct4, nu, ksi, fv1, fv2, ft2, fw, g, r, Omega, S_hat;
            __declspec(align(_ALIGN)) T Wxy, Div_nu;
#else
           T Cb1 __attribute__ ((aligned (_ALIGN)));
           T Cb2 __attribute__ ((aligned (_ALIGN)));
           T sig __attribute__ ((aligned (_ALIGN)));
           T _k __attribute__ ((aligned (_ALIGN)));
           T Cw1 __attribute__ ((aligned (_ALIGN)));
           T Cw2 __attribute__ ((aligned (_ALIGN)));
           T Cw3 __attribute__ ((aligned (_ALIGN)));
           T Cv1 __attribute__ ((aligned (_ALIGN))); /*Ct1,*/
           T Ct2 __attribute__ ((aligned (_ALIGN)));
           T Ct3 __attribute__ ((aligned (_ALIGN)));
           T Ct4 __attribute__ ((aligned (_ALIGN)));
           T nu __attribute__ ((aligned (_ALIGN)));
           T ksi __attribute__ ((aligned (_ALIGN)));
           T fv1 __attribute__ ((aligned (_ALIGN)));
           T fv2 __attribute__ ((aligned (_ALIGN)));
           T ft2 __attribute__ ((aligned (_ALIGN)));
           T fw __attribute__ ((aligned (_ALIGN)));
           T g __attribute__ ((aligned (_ALIGN)));
           T r __attribute__ ((aligned (_ALIGN)));
           T Omega __attribute__ ((aligned (_ALIGN)));
           T S_hat __attribute__ ((aligned (_ALIGN)));
           T Wxy __attribute__ ((aligned (_ALIGN)));
           T Div_nu __attribute__ ((aligned (_ALIGN)));
#endif //__ICC
            fv1 = 1.0;

            if(is_init) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]/100.0; // initial turbulence visc.
            } else if (isCond2D(CT_WALL_NO_SLIP_2D) || isCond2D(CT_WALL_LAW_2D) || 
                       FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = 0.; 
            } else if (isCond2D(NT_FC_2D) ) {
               FlowNodeCore2D<T,a>::S[i2d_nu_t] = mu/FlowNodeCore2D<T,a>::S[i2d_Rho]*5.0;
            } else { // if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_nu_t_CONST_2D)) 
                Cb1 = 0.1355;
                Cb2 = 0.622;
                sig = 2.0/3.0;
                _k  = 0.41;
                Cw1 = Cb1/(_k*_k) + (1 + Cb2)/sig;
                Cw2 = 0.3;
                Cw3 = 2.0;
                Cv1 = 7.1;
                Ct2 = 2.0;
                Ct3 = 1.2; // orig 1.1
                Ct4 = 0.5; // orig 2
                nu  = mu/FlowNodeCore2D<T,a>::S[i2d_Rho];
                ksi = FlowNodeCore2D<T,a>::S[i2d_nu_t]/nu;
                fv1 = ksi*ksi*ksi/(ksi*ksi*ksi + Cv1*Cv1*Cv1);
                fv2 = 1.0 - ksi/(1.0 + ksi * fv1);
                Wxy = 0.5*(dVdx - dUdy);
                Omega = sqrt(2.0*Wxy*Wxy);
                S_hat = Omega + FlowNodeCore2D<T,a>::S[i2d_nu_t]/(_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min) * fv2;

                if(S_hat < 0.3 * Omega)
                   S_hat = 0.3 * Omega;

                r   = min((FlowNodeCore2D<T,a>::S[i2d_nu_t]/(S_hat*_k*_k * FlowNodeTurbulence2D<T,a>::l_min * FlowNodeTurbulence2D<T,a>::l_min)),10.0);
                g   = r + Cw2*(pow(r,6.0) - r);

                fw  = g * pow((1.0 + pow(Cw3,6.0))/(pow(g,6.0)+pow(Cw3,6.0)),1.0/6.0);
                ft2 = Ct2 * exp(-Ct4*ksi*ksi);

                A[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*U;
                B[i2d_nu_t]   = FlowNodeCore2D<T,a>::S[i2d_nu_t]*V;

                Div_nu        =  (FlowNodeTurbulence2D<T,a>::dkdx + FlowNodeTurbulence2D<T,a>::dkdy );

                RX[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdx)/sig;
                RY[i2d_nu_t]  = ((mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdy)/sig;

                A[i2d_nu_t]   = A[i2d_nu_t]  - RX[i2d_nu_t];
                B[i2d_nu_t]   = B[i2d_nu_t]  - RY[i2d_nu_t];

                Src[i2d_nu_t] = Cb1*(1.0 - ft2)*S_hat*FlowNodeCore2D<T,a>::S[i2d_nu_t] - (Cw1*fw - Cb1/(_k*_k)*ft2)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min)*(FlowNodeCore2D<T,a>::S[i2d_nu_t]/FlowNodeTurbulence2D<T,a>::l_min) + (Cb2*Div_nu*Div_nu)/sig;      
            }
            //Axisymmetric turbulence add-on
             TurbulenceAxisymmAddOn(is_init);

            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = FlowNodeCore2D<T,a>::S[i2d_Rho] * FlowNodeCore2D<T,a>::S[i2d_nu_t] * fv1;
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }

        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D)) {

        } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Integral_Model_2D) && mu != 0.0) {
            FlowNodeTurbulence2D<T,a>::Re_local = FlowNodeTurbulence2D<T,a>::l_min*sqrt(FlowNodeCore2D<T,a>::S[i2d_RhoU]*FlowNodeCore2D<T,a>::S[i2d_RhoU]+
                                                                                        FlowNodeCore2D<T,a>::S[i2d_RhoV]*FlowNodeCore2D<T,a>::S[i2d_RhoV]+
                                                                                        1.e-30)/mu; 
        } else if (FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_omega_SST_Model_2D) && !is_init) {

        } else if ( FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Smagorinsky_Model_2D)) {
#ifdef __ICC
          __declspec(align(_ALIGN)) T Omega;
          __declspec(align(_ALIGN)) T Wxy;
          __declspec(align(_ALIGN)) T Cs = 0.1; // Smagorinsky const. Need change for DSGS
#else
          T Omega __attribute__ ((aligned (_ALIGN)));
          T Wxy   __attribute__ ((aligned (_ALIGN)));
          T Cs    __attribute__ ((aligned (_ALIGN))) = 0.1; // Smagorinsky const. Need change for DSGS 
#endif // __ICC
#ifdef    _UNIFORM_MESH_
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( FlowNode2D<T,a>::dx*FlowNode2D<T,a>::dy);
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( FlowNode2D<T,a>::dx*FlowNode2D<T,a>::dy);
 #endif // __ICC
#else
 #ifdef __ICC
           __declspec(align(_ALIGN)) T _delta =sqrt( dx * dy );
 #else
           T _delta __attribute__ ((aligned (_ALIGN))) =sqrt( dx * dy );
 #endif // __ICC
#endif // _UNIFORM_MESH_
          Wxy = 0.5*(dVdx - dUdy);
          Omega = sqrt(2.0*Wxy*Wxy);
            if(is_mu_t) {
              FlowNodeTurbulence2D<T,a>::mu_t  = fabs(FlowNodeCore2D<T,a>::S[i2d_Rho] * (Cs * _delta) * (Cs * _delta) * Omega);
              FlowNodeTurbulence2D<T,a>::lam_t = FlowNodeTurbulence2D<T,a>::mu_t * CP;
            }
        }
}

template <class T, int a>
void FlowNode2D<T,a>::TurbulenceAxisymmAddOn(int is_init) {
      if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_k_eps_Model_2D) && !is_init) {
        F[i2d_k]   = FT*(mu+FlowNodeTurbulence2D<T,a>::mu_t)*FlowNodeTurbulence2D<T,a>::dkdy;
        F[i2d_eps] = FT*(mu+FlowNodeTurbulence2D<T,a>::mu_t/1.3)*FlowNodeTurbulence2D<T,a>::depsdy;
      } else if(FlowNodeTurbulence2D<T,a>::isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) && !is_init) {
        F[i2d_nu_t]   = FT*(mu/FlowNodeCore2D<T,a>::S[i2d_Rho]+FlowNodeCore2D<T,a>::S[i2d_nu_t])*FlowNodeTurbulence2D<T,a>::dkdy;
      } else if(is_init) {
        F[i2d_nu_t]   = F[i2d_eps] = 0.;
        Src[i2d_nu_t] = Src[i2d_eps] = 0.;
      }
}

#endif // _CUDA_

template <class T, int a>
inline FlowNode2D<T,a>& FlowNode2D<T,a>::operator = (FlowNode2D<T,a>& fn) {
    memcpy(this,&fn,sizeof(FlowNode2D<T,a>));return *this;
}

template <class T, int a>
inline FlowNode2D<T,a>& FlowNode2D<T,a>::operator = (FlowNodeCore2D<T,a>& fc) {
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(register int i=0;i<NumEq;i++) {
        FlowNodeCore2D<T,a>::S[i]=fc.FlowNodeCore2D<T,a>::S[i];
        FlowNodeCore2D<T,a>::dSdx[i]=fc.FlowNodeCore2D<T,a>::dSdx[i];
        FlowNodeCore2D<T,a>::dSdy[i]=fc.FlowNodeCore2D<T,a>::dSdy[i];
    }
    //FillNode2D(0,1);
    return *this;
}


template <class T, int a>
FlowNode2D<T,a>& FlowNode2D<T,a>::operator = (Flow& f) {
    register T   Tmp1,Tmp3=0.;
    register int i;

    FlowNodeCore2D<T,a>::S[i2d_RhoU] = f.ROG();
    FlowNodeCore2D<T,a>::S[i2d_RhoV] = f.ROG()*f.Wg();
    FlowNodeCore2D<T,a>::S[i2d_RhoE] = FlowNodeCore2D<T,a>::S[a+4] = FlowNodeCore2D<T,a>::S[a+5] = 0.;

    p    = f.P0();
    R    = f.Rg();
    lam  = f.lam;
    mu   = f.mu;
    CP   = f.C;
    k = CP/(CP-R);
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++)
        FlowNodeCore2D<T,a>::S[i+4]=FlowNodeCore2D<T,a>::S[i2d_Rho]*Y[i];

    Tmp1 = FlowNodeCore2D<T,a>::S[i2d_Rho];
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Tmp3+=Hu[i]*FlowNodeCore2D<T,a>::S[i+4];
        Tmp1-= FlowNodeCore2D<T,a>::S[i+4];
    }

    Tmp3+=Hu[a]*Tmp1;
    FlowNodeCore2D<T,a>::S[i2d_Rho]  = p/R/Tg;
    FlowNodeCore2D<T,a>::S[i2d_RhoE] = p/(k-1.)+FlowNodeCore2D<T,a>::S[i2d_Rho]*(f.Wg()*f.Wg())*0.5+Tmp3;
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<NumEq;i++)
        Src[i]=0;

    //FillNode2D(0,1);
    return *this;
}


template <class T, int a>
FlowNode2D<T,a>& FlowNode2D<T,a>::operator = (Flow2D& f) {
    register T Tmp1,Tmp3=0.;
    register int i;

    FlowNodeCore2D<T,a>::S[i2d_Rho] = f.ROG();
    FlowNodeCore2D<T,a>::S[i2d_RhoU] = f.ROG()*f.U();
    FlowNodeCore2D<T,a>::S[i2d_RhoV] = f.ROG()*f.V();
    FlowNodeCore2D<T,a>::S[i2d_k] = FlowNodeCore2D<T,a>::S[i2d_eps] = 0.;
    U    = f.U();
    V    = f.V();
    p    = f.Pg();
    R    = f.Rg();
    lam  = f.lam;
    mu   = f.mu;
    Tg   = f.Tg();
    CP   = f.C;
    k = CP/(CP-R);

    Diff   =  (lam+FlowNodeTurbulence2D<T,a>::lam_t)/CP;

    FlowNodeCore2D<T,a>::S[i2d_Rho] = f.Pg()/R/f.Tg();
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++)
        FlowNodeCore2D<T,a>::S[i+4] = Y[i]*FlowNodeCore2D<T,a>::S[i2d_Rho];

    Tmp1 = FlowNodeCore2D<T,a>::S[i2d_Rho];
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<a;i++) {
        Tmp3+=Hu[i]*FlowNodeCore2D<T,a>::S[i+4];
        Tmp1-= FlowNodeCore2D<T,a>::S[i+4];
    }

    Tmp3+=Hu[a]*Tmp1;
    FlowNodeCore2D<T,a>::S[i2d_RhoE] = p/(k-1)+FlowNodeCore2D<T,a>::S[i2d_Rho]*(f.U()*f.U()+f.V()*f.V())*0.5+Tmp3;
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(i=0;i<NumEq;i++)
        Src[i]=0;

    //FillNode2D(0,1);
    return *this;
}

template <class T, int a>
inline FlowNode2D<T,a>& FlowNode2D<T,a>::operator *= (T& m) {
    
#ifdef _CUDA_
#pragma unroll
#endif // _CUDA_
    for(register int i=i2d_RhoU;i<i2d_RhoE;i++) {
        FlowNodeCore2D<T,a>::S[i]=FlowNodeCore2D<T,a>::S[i]*m;
        FlowNodeCore2D<T,a>::dSdx[i]=FlowNodeCore2D<T,a>::dSdx[i]*m;
        FlowNodeCore2D<T,a>::dSdy[i]=FlowNodeCore2D<T,a>::dSdy[i]*m;
    }
    //FillNode2D(0,1);
    return *this;
}

//---------- MakeFlow2D ----------------------------------
template <class T, int a>
inline Flow2D* FlowNode2D<T,a>::MakeFlow2D(Flow2D* f) {
    f->kg(k);
    f->Rg(R);
    f->lam=lam;
    f->mu=mu;
    f->C  =  CP;
    f->P0(p+FlowNodeCore2D<T,a>::S[i2d_Rho]*(sqrt(U*U+V*V+1.e-10))/2);
    f->T0(f->P0()/FlowNodeCore2D<T,a>::S[i2d_Rho]/R);
    f->ROG(FlowNodeCore2D<T,a>::S[i2d_Rho]);
    f->U(U);
    f->V(V);
    return f;
}

template <class T, int a>
inline void FlowNode2D<T,a>::CopyFlowNodeCore2D(FlowNodeCore2D<T,a>& fnc)  {
 fnc = *this;
}

// Function for operations whis CondType bit flags without FlowNode2D class
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
inline ulong  isCond2D(ulong ct, CondType2D ctt) {
    return ((ct & ctt) == (ulong)ctt);
}
#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
inline ulong  isTurbulenceCond2D(ulong ct, TurbulenceCondType2D ctt) {
    return ((ct & ctt) == (ulong)ctt);
}

// <------------- 2D --------------->
#endif // _hyper_flow_node_hpp

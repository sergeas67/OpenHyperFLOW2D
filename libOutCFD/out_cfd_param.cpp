/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  2.0.1                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   last update: 14/04/2016                                                    *
*******************************************************************************/

#include "libOutCFD/out_cfd_param.hpp"

FP Re_Airfoil(FP chord, Flow2D* Flow) {
 return Flow->Wg()*chord*Flow->ROG()/Flow->mu;
}


FP Calc_I(FlowNode2D<FP,NUM_COMPONENTS>* node)
{
    return 0;
}

FP p_asterisk(FlowNode2D<FP,NUM_COMPONENTS>* node ) {
    FP A,Mach,WW,tmp_pow;
    A  = sqrt(node->k*node->R*node->Tg);
    WW = sqrt(node->U*node->U +
              node->V*node->V);
    Mach = WW/A;
    tmp_pow = node->k/(node->k-1.0);
    return  node->p * pow(1.0+(node->k-1.0)*0.5*Mach*Mach,tmp_pow);
}

FP Schliren(FlowNode2D<FP,NUM_COMPONENTS>* node ) {
   return sqrt(node->dSdx[i2d_Rho]*node->dSdx[i2d_Rho] + node->dSdy[i2d_Rho]*node->dSdy[i2d_Rho]);
}


FP T_asterisk(FlowNode2D<FP,NUM_COMPONENTS>* node ) {
    FP T_a;
    if (node->CP > 0.)
        T_a= (node->U*node->U+node->V*node->V)*0.5/node->CP;
    else
        T_a = 0.;
    return T_a;
}

FP CalcaveragePressure2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                             FP x0, // initial point of probed area
                             FP l,  // length of probed area
                             FP d   // diameter (or cross-section size) of probed area
                            ) {
    unsigned int n_pressure=0;
    FP p_mid = 0.;
    FP V_sum = 0;

#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:n_pressure,p_mid,V_sum)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:n_pressure,p_mid,V_sum)
#endif //_OPEN_MP 
#endif // _MPI_OPENMP 
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D) && 
                 i > (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i < (int)((l+x0)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j < (int)(d/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {
                FP V_i = 0;
                if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
                    V_i = 2 * M_PI * pJ->GetValue(i,j).y * FlowNode2D<FP,NUM_COMPONENTS>::dy * FlowNode2D<FP,NUM_COMPONENTS>::dx; 
                    V_sum += V_i;
                    p_mid += pJ->GetValue(i,j).p * V_i;
                } else {
                    p_mid += pJ->GetValue(i,j).p;
                }
                n_pressure++;
            }
        }

    if (n_pressure) {
        if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
            return p_mid/V_sum;
        } else {
            return p_mid/n_pressure;
        }
    } else {
        return 0.;
    }
}

FP CalcaverageTemperature2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                                FP x0, // initial point of probed area
                                FP l,  // length of probed area
                                FP d,  // diameter (or cross-section size) of probed area
                                int is_mid_ethalpy_T
                               ) {
    unsigned int n_temperature=0;
    FP T_mid=0.;
    FP V_sum = 0;

#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:n_temperature,T_mid, V_sum)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:n_temperature,T_mid, V_sum)
#endif //_OPEN_MP 
#endif // _MPI_OPENMP 
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D) && 
                 i > (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i < (int)((l+x0)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j < (int)(d/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {
                FP V_i = 0;
                if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
                    V_i = 2 * M_PI * pJ->GetValue(i,j).y * FlowNode2D<FP,NUM_COMPONENTS>::dy * FlowNode2D<FP,NUM_COMPONENTS>::dx; 

                    if (is_mid_ethalpy_T)
                        V_i = V_i*pJ->GetValue(i,j).CP;

                    V_sum += V_i;
                    T_mid += pJ->GetValue(i,j).Tg * V_i;
                } else {
                    T_mid += pJ->GetValue(i,j).Tg;
                }
                n_temperature++;    
            }
        }

    if (n_temperature) {
        if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
            return  T_mid/V_sum;
        } else {
            return T_mid/n_temperature;
        }
    } else {
        return 0.;
    }
}

FP CalcArea2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,             
              FP x0,  // initial X  point of probed area                  
              FP y0,  // initial Y  point of probed area                  
              FP dy   // radius (or cross-section size) of probed area 
              ) {                                                             
FP Sp=0.;
unsigned int i = (unsigned int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx);
unsigned int jj_start = (unsigned int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy);
unsigned int jj_end = (unsigned int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy);
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Sp)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Sp)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int j=jj_start;j<(int)jj_end;j++ ) {
        if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
            if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT)
                Sp+= FlowNode2D<FP,NUM_COMPONENTS>::dy;
            else
                Sp+= 2 * M_PI *FlowNode2D<FP,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).y;
        }
    }
    return(Sp);
}

FP CalcMassFlowRateX2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                       FP x0,  // initial X  point of probed area
                       FP y0,  // initial Y  point of probed area
                       FP dy   // diameter (or cross-section size) of probed area) 
                       ) {
FP Mp=0.;
unsigned int i = (unsigned int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx);
unsigned int jj_start = (unsigned int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy);
unsigned int jj_end = (unsigned int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy);
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Mp)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Mp)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int j=jj_start;j<(int)jj_end;j++ ) {
        if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
            if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT)
                Mp+= FlowNode2D<FP,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).S[i2d_RhoU];
            else
                Mp+= 2 * M_PI *FlowNode2D<FP,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).y*pJ->GetValue(i,j).S[i2d_RhoU];
        }
    }
    return(Mp);
}


FP CalcXForceYSym2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // initial point of probed area
                    FP l,  // length of probed area
                    FP d   // diameter (or cross-section size) of probed area
                    ) {
    FP Fp=0.; // pressure force
    FP Fd=0.; // drag force
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&
                 i >= (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i <= (int)((l+x0)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j <= (int)(d/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {

                FP Sp = 0.;
                FP Sd = 0.;

                if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<FP,NUM_COMPONENTS>::dy;
                    Sd = FlowNode2D<FP,NUM_COMPONENTS>::dx;
                } else {
                    Sp = 2 * M_PI * pJ->GetValue(i,j).y*FlowNode2D<FP,NUM_COMPONENTS>::dy;
                    Sd = 2 * M_PI * pJ->GetValue(i,j).y*FlowNode2D<FP,NUM_COMPONENTS>::dx;
                }
                // Pressure forces
                if (i > 0 && pJ->GetValue(i-1,j).isCond2D(CT_SOLID_2D) ) {
                    Fp+= Sp*pJ->GetValue(i,j).p;
                } else if ( i < (int)pJ->GetX()-1  && pJ->GetValue(i+1,j).isCond2D(CT_SOLID_2D)) {
                    Fp-= Sp*pJ->GetValue(i,j).p;
                }
                // Drag forces
                if (j < (int)pJ->GetY()-1 && !pJ->GetValue(i,j+1).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i,j+1).U > 0)
                        Fd+= - Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                    else
                        Fd-= - Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                } else if (j > 0 && !pJ->GetValue(i,j-1).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i,j-1).U > 0)
                        Fd+= - Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                    else
                        Fd-= - Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                }

            }
        }
    return(Fp+Fd);
}

FP CalcXForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                FP x0, // X initial point of probed area
                FP y0, // Y initial point of probed area
                FP dx, // X size of probed area
                FP dy  // Y of probed area
                ) {
    FP Fp=0.; // pressure force
    FP Fd=0.; // drag force
#ifdef _MPI_OPENMP
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP
#endif // _MPI_OPENMP
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&
                 i >= (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {

                FP Sp = 0.;
                FP Sd = 0.;

                if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<FP,NUM_COMPONENTS>::dy;
                    Sd = FlowNode2D<FP,NUM_COMPONENTS>::dx;
                } else {
                    Sp = 2 * M_PI * (j+0.5)*FlowNode2D<FP,NUM_COMPONENTS>::dy*FlowNode2D<FP,NUM_COMPONENTS>::dy;
                    Sd = 2 * M_PI * (j+0.5)*FlowNode2D<FP,NUM_COMPONENTS>::dy*FlowNode2D<FP,NUM_COMPONENTS>::dx;
                }
                // Pressure forces

                if (i > 0 && pJ->GetValue(i-1,j).isCond2D(CT_SOLID_2D)) {                         // [solid]<-[gas] -force
                    Fp -= Sp*pJ->GetValue(i,j).p;
                } else if ( i < (int)pJ->GetX()-1  && pJ->GetValue(i+1,j).isCond2D(CT_SOLID_2D)) { // [gas]->[solid] +force
                    Fp += Sp*pJ->GetValue(i,j).p;
                }

                // Drag forces
                if (j < (int)pJ->GetY()-1 && !pJ->GetValue(i,j+1).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i,j+1).U > 0)                                                 // flow --> +force
                        Fd+= Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                    else                                                                           // flow <-- -force
                        Fd-= Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                } else if (j > 0 && !pJ->GetValue(i,j-1).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i,j-1).U > 0)                                                 // flow --> +force
                        Fd+= Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                    else                                                                           // flow <-- -force
                        Fd-= Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dUdy);
                }

            }
        }

    //cout << "\nCalc Cx\nPressure force (X):" << Fp << " N" << endl;
    //cout << "Drag force (X):" << Fd << " N" << endl;

    return(Fp+Fd);
}

FP CalcYForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // X initial point of probed area
                    FP y0, // Y initial point of probed area
                    FP dx, // X size of probed area
                    FP dy  // Y of probed area
                   ) {
    FP Fp=0.; // pressure force
    FP Fd=0.; // drag force
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP 
#endif // _MPI_OPENMP 
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&    
                 i >= (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {

                FP Sp = 0.;
                FP Sd = 0.;

                if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<FP,NUM_COMPONENTS>::dx;
                    Sd = FlowNode2D<FP,NUM_COMPONENTS>::dy;
                } else {
                    Sp = 2 * M_PI * pJ->GetValue(i,j).y*FlowNode2D<FP,NUM_COMPONENTS>::dx;
                    Sd = 2 * M_PI * pJ->GetValue(i,j).y*FlowNode2D<FP,NUM_COMPONENTS>::dy;
                }
                // Pressure forces

                if (j > 0 && pJ->GetValue(i,j-1).isCond2D(CT_SOLID_2D) ) { //  [gas]
                    Fp-= Sp*pJ->GetValue(i,j).p;                           // [solid] down force (+)
                } else if ( j < (int)pJ->GetY()-1  && pJ->GetValue(i,j+1).isCond2D(CT_SOLID_2D)) {
                    Fp+= Sp*pJ->GetValue(i,j).p;                           //  [gas]  up force (-)
                }

                // Drag forces

                if (i < (int)pJ->GetX()-1 && !pJ->GetValue(i+1,j).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i+1,j).V > 0)
                        Fd+= -Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dVdx);
                    else
                        Fd-= -Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dVdx);
                } else if (i > 0 && !pJ->GetValue(i-1,j).isCond2D(CT_SOLID_2D)) {
                    if (pJ->GetValue(i-1,j).V > 0)
                        Fd+= -Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dVdx);
                    else
                        Fd-= -Sd*(pJ->GetValue(i,j).mu+pJ->GetValue(i,j).mu_t) * fabs(pJ->GetValue(i,j).dVdx);
                }

            }
        }
    
    //cout << "\nCalc Cy\nPressure force (Y):" << Fp << " N" << endl;
    //cout << "Drag force (Y):" << Fd << " N" << endl;
    return(Fp+Fd);
}

FP   Calc_Cp(FlowNode2D<FP,NUM_COMPONENTS>* CurrentNode, Flow2D* pF) {
    if (CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D) && pF->Wg() > 0.0)
        return(CurrentNode->p - pF->Pg())/(0.5*pF->ROG()*pF->Wg()*pF->Wg());
    else
        return 0;
}

FP GetFmid(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
           FP x0, 
           FP y0,
           FP dx, 
           FP dy) {
    
    FP Fmid=0;
    int    is_node;
#ifdef _MPI_OPENMP
#pragma omp parallel for reduction(+:Fmid)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fmid)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int j=0;j<(int)pJ->GetY();j++ ) {
        is_node = 0;
        for (int i=0;i<(int)pJ->GetX();i++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D) ) &&
                 i >= (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {
                 is_node = 1;
            }
        }
        if (is_node) {
            if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT) {
                Fmid += FlowNode2D<FP,NUM_COMPONENTS>::dy;
            } else {
                Fmid += 2 * M_PI * (j+0.5)*FlowNode2D<FP,NUM_COMPONENTS>::dy*FlowNode2D<FP,NUM_COMPONENTS>::dy;
            }
        }
    }
    //cout << "\nFmid=" << Fmid << endl;
    return Fmid;
}

FP GetS(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
        FP x0, 
        FP y0,
        FP dx, 
        FP dy) {
    
    FP S=0;
    int    is_node;
#ifdef _MPI_OPENMP
#pragma omp parallel for reduction(+:Fmid)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fmid)
#endif //_OPEN_MP
#endif // _MPI_OPENMP
    for (int i=0;i<(int)pJ->GetX();i++ ) {
        is_node = 0;
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D) ) &&
                 i >= (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<FP,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy)) {
                 is_node = 1;
            }
        }
        if (is_node) {
            S += FlowNode2D<FP,NUM_COMPONENTS>::dx;
        }
    }
    //cout << "\nS=" << S << endl;
    return S;
}

FP Calc_Cx_2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
              FP x0, 
              FP y0,
              FP dx, 
              FP dy,
              Flow2D* pF
             ) {

FP Pmax = pF->ROG()*pF->Wg()*pF->Wg()*0.5*GetS(pJ,x0,y0,dx,dy);
//cout << "\nPmax=" << Pmax << endl;
if (Pmax == 0.)
        return 0;
    else
        return CalcXForce2D(pJ,x0,y0,dx,dy)/Pmax;
}

FP   Calc_Cy_2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ, 
                    FP x0, 
                    FP y0,
                    FP dx, 
                    FP dy,
                    Flow2D* pF) {

FP Pmax = pF->ROG()*pF->Wg()*pF->Wg()*0.5*GetS(pJ,x0,y0,dx,dy);
//cout << "\nPmax=" << Pmax << endl;    
if (Pmax == 0.)
        return 0;
    else
        return CalcYForce2D(pJ,x0,y0,dx,dy)/Pmax;

    return 0;
}


void SmoothY(UMatrix2D<FP>* A) {
    for(int j=0;j<(int)A->GetY();j++) {
      for(int i=0;i<(int)A->GetX();i++) {
          if(j > 0 && j <  (int)A->GetY()-1 &&
             A->GetValue(i,j+1) >  0.   &&
             A->GetValue(i,j-1) >  0. ) {
             A->GetValue(i,j) = 0.5*(A->GetValue(i,j+1)+A->GetValue(i,j-1));
        }
      }
    }
}

void SmoothX(UMatrix2D<FP>* A) {
    for(int j=0;j<(int)A->GetY();j++) {
      for(int i=0;i<(int)A->GetX();i++) {
          if(i > 0 && i <  (int)A->GetX()-1 &&
             A->GetValue(i+1,j) >  0.   &&
             A->GetValue(i-1,j) >  0. ) {
             A->GetValue(i,j) = 0.5*(A->GetValue(i+1,j)+A->GetValue(i-1,j));
        }
      }
    }
}

void SaveXHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                     Flow2D* TestFlow, 
                     FP  Ts,
                     int     y_max,
                     int     y_min) {
    char  HeatFluxHeader[128];
    FP Cp;
    FP St;
    FP Trec  = (1 + 0.45 * (TestFlow->kg() - 1.0) * TestFlow->MACH() * TestFlow->MACH())*TestFlow->Tg();

#ifdef _REF_TEST_
    snprintf(HeatFluxHeader,128,"#VARIABLES = X, HeatFlux(X), Alpha(X), HeatFluxRef(X), AlphaRef(X), Re(X), Pr(X)");
    UArray<FP> HeatFluxRefArray(pJ->GetX());
    UArray<FP> HeatExcgCoeffRefArray(pJ->GetX());

    UArray<FP> ReArray(pJ->GetX());
    UArray<FP> PrArray(pJ->GetX());
    
#else
    snprintf(HeatFluxHeader,128,"#VARIABLES = X, HeatFlux(X),  Alpha(X), Cp(X), St(X)");
    
    UArray<FP> CpArray(pJ->GetX());
    UArray<FP> StArray(pJ->GetX());

#endif // _REF_TEST_
    
    *OutputData << HeatFluxHeader  << endl;

    UArray<FP> HeatFluxArray(pJ->GetX());
    UArray<FP> HeatExcgCoeffArray(pJ->GetX());

    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        FP Zero = 0.0;
        
        HeatFluxArray.SetElement(i,&Zero);
        HeatExcgCoeffArray.SetElement(i,&Zero);
        
        CpArray.SetElement(i,&Zero);
        StArray.SetElement(i,&Zero);

#ifdef _REF_TEST_
        ReArray.SetElement(i,&Zero);
        PrArray.SetElement(i,&Zero);
        
        HeatFluxRefArray.SetElement(i,&Zero);
        HeatExcgCoeffRefArray.SetElement(i,&Zero);
#endif // _REF_TEST_
    }

    for(int i=0; i < (int)pJ->GetX(); i++) {
        for(int j=max(0,y_min); j < min(y_max,(int)pJ->GetY()-1); j++) {

        if(pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {

            FlowNode2D<FP,NUM_COMPONENTS>* CurrentNode;
            FlowNode2D<FP,NUM_COMPONENTS>* UpNode;
            FlowNode2D<FP,NUM_COMPONENTS>* DownNode; 
            FlowNode2D<FP,NUM_COMPONENTS>* RightNode;
            FlowNode2D<FP,NUM_COMPONENTS>* LeftNode; 
            FP lam_eff;
            int    num_near_nodes = 1;

            lam_eff = pJ->GetValue(i,j).lam + pJ->GetValue(i,j).lam_t;

            int n1,n2,n3,n4,N1,N2,N3,N4;

            n1=pJ->GetValue(i,j).idXl;
            n2=pJ->GetValue(i,j).idXr;
            n3=pJ->GetValue(i,j).idYu;
            n4=pJ->GetValue(i,j).idYd;

            N1 = i - n1;
            N2 = i + n2;                                                                                                                                                                                                                                                                          N3 = j + n3;                                                                                                                   
            N3 = j + n3;
            N4 = j - n4;

            CurrentNode= &(pJ->GetValue(i,j));
            UpNode     = &(pJ->GetValue(i,N3));
            DownNode   = &(pJ->GetValue(i,N4));
            RightNode  = &(pJ->GetValue(N2,j));
            LeftNode   = &(pJ->GetValue(N1,j));

            if(LeftNode) {
              lam_eff += LeftNode->lam + LeftNode->lam_t;
              num_near_nodes++;
            }
            if (RightNode) {
                lam_eff += RightNode->lam + RightNode->lam_t;
                num_near_nodes++;
            }

            if (UpNode) {
                lam_eff +=UpNode->lam + UpNode->lam_t;
                num_near_nodes++;
            }
            if (DownNode) {
                lam_eff += DownNode->lam + DownNode->lam_t;
                num_near_nodes++;
            }

            lam_eff = lam_eff/num_near_nodes;
            FP Q     = lam_eff*(pJ->GetValue(i,j).Tg - Ts)/FlowNode2D<FP,NUM_COMPONENTS>::dy;
            FP alpha = lam_eff/FlowNode2D<FP,NUM_COMPONENTS>::dy;
            
            Cp = St = 0.0;

#ifdef _REF_TEST_
            FP Re    = (pJ->GetValue(i,pJ->GetY()-1).U * (i+0.5)*FlowNode2D<FP,NUM_COMPONENTS>::dx*pJ->GetValue(i,j).S[0])/pJ->GetValue(i,j).mu;
            FP Pr    =  pJ->GetValue(i,j).mu*pJ->GetValue(i,j).CP/pJ->GetValue(i,j).lam;
            FP Nu;
            
            ReArray.SetElement(i,&Re);
            PrArray.SetElement(i,&Pr);
            
            if(Re < 5.0e5 ) {
              Nu = 0.332*sqrt(Re)*pow(Pr,1.0/3.0);
            } else {
              Nu = 0.0296*pow(Re,0.8)*pow(Pr,1.0/3.0);
            }
            
            FP Alpha_Ref =  Nu*pJ->GetValue(i,j).lam/((i+0.5)*FlowNode2D<FP,NUM_COMPONENTS>::dx);
            FP Q_Ref     =  Alpha_Ref*(pJ->GetValue(i,j).Tg - Ts);
#endif // _REF_TEST_            
            St =     Q/(TestFlow->ROG()*TestFlow->Wg()*TestFlow->C*(Trec-Ts));
            Cp =     Calc_Cp(CurrentNode,TestFlow);
            
            if(HeatFluxArray.GetElement(i) != 0.) {
              Q  =     max(HeatFluxArray.GetElement(i),Q);
#ifdef _REF_TEST_
              Q_Ref = max(HeatFluxRefArray.GetElement(i),Q_Ref);
              HeatFluxRefArray.SetElement(i,&Q_Ref);
              Alpha_Ref = max(HeatExcgCoeffRefArray.GetElement(i),Alpha_Ref);
              HeatExcgCoeffRefArray.SetElement(i,&Alpha_Ref); 
#endif // _REF_TEST_            
              
              HeatFluxArray.SetElement(i,&Q);
              
              alpha = max(HeatExcgCoeffArray.GetElement(i),alpha);

              HeatExcgCoeffArray.SetElement(i,&alpha);

              CpArray.SetElement(i,&Cp);
              StArray.SetElement(i,&St);
            } else {
              HeatFluxArray.SetElement(i,&Q);
              HeatExcgCoeffArray.SetElement(i,&alpha);
              
              CpArray.SetElement(i,&Cp);
              StArray.SetElement(i,&St);
#ifdef _REF_TEST_
              HeatExcgCoeffRefArray.SetElement(i,&Alpha_Ref); 
              HeatFluxRefArray.SetElement(i,&Q_Ref);
#endif // _REF_TEST_            
            }
          }
        }
    }
    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        *OutputData << i*FlowNode2D<FP,NUM_COMPONENTS>::dx << " " << HeatFluxArray.GetElement(i)<< " " << HeatExcgCoeffArray.GetElement(i) << 
                                                                  " " << CpArray.GetElement(i) << " " << StArray.GetElement(i)  <<
#ifdef _REF_TEST_
                                                                  " " << HeatFluxRefArray.GetElement(i) << " " << HeatExcgCoeffRefArray.GetElement(i) << 
                                                                  " " << ReArray.GetElement(i) << " " << PrArray.GetElement(i)  <<
#endif // _REF_TEST_ 
                                                                  endl;           
    }
}

void SaveYHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                     FP  Ts) {
    char  HeatFluxHeader[128];
    snprintf(HeatFluxHeader,128,"#VARIABLES = Y, HeatFlux(Y)");
    *OutputData << HeatFluxHeader  << endl;
    UArray<FP> HeatFluxArray(pJ->GetY());
    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        FP Zero = 0.0;
        HeatFluxArray.SetElement(i,&Zero); 
    }
    for(int j=0; j < (int)pJ->GetY(); j++) {
         for(int i=0; i < (int)pJ->GetX()-1; i++){
        if(pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {
            FP lam_eff = pJ->GetValue(i,j).lam + pJ->GetValue(i,j).lam_t;
            int    num_near_nodes = 1;

            int n1,n2,n3,n4,N1,N2,N3,N4;

            n1=pJ->GetValue(i,j).idXl;
            n2=pJ->GetValue(i,j).idXr;
            n3=pJ->GetValue(i,j).idYu;
            n4=pJ->GetValue(i,j).idYd;

            N1 = i - n1;
            N2 = i + n2;                                                                                                                                                                                                                                                                          N3 = j + n3;                                                                                                                   
            N3 = j + n3;
            N4 = j - n4;

            FlowNode2D<FP,NUM_COMPONENTS>* UpNode     = &(pJ->GetValue(i,N3));
            FlowNode2D<FP,NUM_COMPONENTS>* DownNode   = &(pJ->GetValue(i,N4));
            FlowNode2D<FP,NUM_COMPONENTS>* RightNode  = &(pJ->GetValue(N2,j));
            FlowNode2D<FP,NUM_COMPONENTS>* LeftNode   = &(pJ->GetValue(N1,j));

            if(LeftNode) {
              lam_eff += LeftNode->lam + LeftNode->lam_t;
              num_near_nodes++;
            }
            if (RightNode) {
                lam_eff += RightNode->lam + RightNode->lam_t;
                num_near_nodes++;
            }

            if (UpNode) {
                lam_eff +=UpNode->lam + UpNode->lam_t;
                num_near_nodes++;
            }
            if (DownNode) {
                lam_eff += DownNode->lam + DownNode->lam_t;
                num_near_nodes++;
            }

            lam_eff = lam_eff/num_near_nodes;

            FP Q = lam_eff*(pJ->GetValue(i,j).Tg - Ts)/FlowNode2D<FP,NUM_COMPONENTS>::dx;
            if(HeatFluxArray.GetElement(j) != 0.) {
              Q =  max(HeatFluxArray.GetElement(j),Q);
              HeatFluxArray.SetElement(j,&Q); 
            } else {
              HeatFluxArray.SetElement(j,&Q); 
            }
           }
        }
    }
    for(int j=0; j < (int)HeatFluxArray.GetNumElements(); j++) {
        *OutputData << j*FlowNode2D<FP,NUM_COMPONENTS>::dy << " " << HeatFluxArray.GetElement(j) << endl;
    }
}

FP Calc_Cv(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,               
           FP x0, // initial point of probed area                        
           FP y0,                                                        
           FP dy, // diameter (or cross-section size) of probed area)
           FP p_amb,
           Flow2D* pF) {
  
FP Fv=0.;
FP Mp=0;
unsigned int i = (unsigned int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx);
unsigned int jj_start = (unsigned int)(y0/FlowNode2D<FP,NUM_COMPONENTS>::dy);
unsigned int jj_end = (unsigned int)((y0+dy)/FlowNode2D<FP,NUM_COMPONENTS>::dy);
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Fv)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fv)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int j=jj_start;j<(int)jj_end;j++ ) {
        if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
            if (FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_FLAT) {
                Fv+= FlowNode2D<FP,NUM_COMPONENTS>::dy*(pJ->GetValue(i,j).S[i2d_RhoU]*pJ->GetValue(i,j).U+(pJ->GetValue(i,j).p - p_amb));
            } else {
                Fv+= 2 * M_PI *FlowNode2D<FP,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).y*(pJ->GetValue(i,j).S[i2d_RhoU]*pJ->GetValue(i,j).U+(pJ->GetValue(i,j).p - p_amb));
            }
        }
    }
    
    Mp = CalcMassFlowRateX2D(pJ,x0,y0,dy);
    
    if(Mp > 0.0)
       return Fv/(pF->U()*Mp);
    else
       return 0;
}
                                                                                 

FP Calc_Cd(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,               
           FP x0, // initial point of probed area                        
           FP y0,                                                        
           FP dy, // radius (or cross-section size) of probed area   
           Flow2D* pF) {

    return CalcMassFlowRateX2D(pJ,x0,y0,dy)/pF->ROG()/pF->Wg()/CalcArea2D(pJ,x0,y0,dy);

}


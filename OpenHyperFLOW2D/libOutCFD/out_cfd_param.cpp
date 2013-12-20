/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 15/12/2013                                                    *
*******************************************************************************/

#include "libOutCFD/out_cfd_param.hpp"

double Calc_I(FlowNode2D<double,NUM_COMPONENTS>* node)
{
    return 0;
}

double p_asterisk(FlowNode2D<double,NUM_COMPONENTS>* node ) {
    double A,Mach,WW,tmp_pow;
    A  = sqrt(node->k*node->R*node->Tg);
    WW = sqrt(node->U*node->U +
              node->V*node->V);
    Mach = WW/A;
    tmp_pow = node->k/(node->k-1.0);
    return  node->p * pow(1.0+(node->k-1.0)*0.5*Mach*Mach,tmp_pow);
}

double Schliren(FlowNode2D<double,NUM_COMPONENTS>* node ) {
   return sqrt(node->dSdx[i2d_Ro]*node->dSdx[i2d_Ro] + node->dSdy[i2d_Ro]*node->dSdy[i2d_Ro]);
}


double T_asterisk(FlowNode2D<double,NUM_COMPONENTS>* node ) {
    double T_a;
    if (node->CP > 0.)
        T_a= (node->U*node->U+node->V*node->V)*0.5/node->CP;
    else
        T_a = 0.;
    return T_a;
}

double CalcaveragePressure2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                             double x0, // initial point of probed area
                             double l,  // length of probed area
                             double d   // diameter (or cross-section size) of probed area
                            ) {
    unsigned int n_pressure=0;
    double p_mid = 0.;
    double V_sum = 0;

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
                 i > (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i < (int)((l+x0)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j < (int)(d/FlowNode2D<double,NUM_COMPONENTS>::dy)) {
                double V_i = 0;
                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
                    V_i = 2 * M_PI * pJ->GetValue(i,j).r * FlowNode2D<double,NUM_COMPONENTS>::dy * FlowNode2D<double,NUM_COMPONENTS>::dx; 
                    V_sum += V_i;
                    p_mid += pJ->GetValue(i,j).p * V_i;
                } else {
                    p_mid += pJ->GetValue(i,j).p;
                }
                n_pressure++;
            }
        }

    if (n_pressure) {
        if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
            return p_mid/V_sum;
        } else {
            return p_mid/n_pressure;
        }
    } else {
        return 0.;
    }
}

double CalcaverageTemperature2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                                double x0, // initial point of probed area
                                double l,  // length of probed area
                                double d,  // diameter (or cross-section size) of probed area
                                int is_mid_ethalpy_T
                               ) {
    unsigned int n_temperature=0;
    double T_mid=0.;
    double V_sum = 0;

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
                 i > (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i < (int)((l+x0)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j < (int)(d/FlowNode2D<double,NUM_COMPONENTS>::dy)) {
                double V_i = 0;
                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
                    V_i = 2 * M_PI * pJ->GetValue(i,j).r * FlowNode2D<double,NUM_COMPONENTS>::dy * FlowNode2D<double,NUM_COMPONENTS>::dx; 

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
        if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC) {
            return  T_mid/V_sum;
        } else {
            return T_mid/n_temperature;
        }
    } else {
        return 0.;
    }
}

double CalcMassFlowRateX2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                           double x0,  // initial X  point of probed area
                           double y0,  // initial Y  point of probed area
                           double dy   // diameter (or cross-section size) of probed area) 
                          ){
    double Mp=0.;
    unsigned int i = (unsigned int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx);
    unsigned int jj_start = (unsigned int)(y0/FlowNode2D<double,NUM_COMPONENTS>::dy);
    unsigned int jj_end = (unsigned int)((y0+dy)/FlowNode2D<double,NUM_COMPONENTS>::dy);
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Mp)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Mp)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int j=jj_start;j<(int)jj_end;j++ ) {
        if ( !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
            if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT)
                Mp+= FlowNode2D<double,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).S[i2d_RoU];
            else
                Mp+= 2 * M_PI *FlowNode2D<double,NUM_COMPONENTS>::dy*pJ->GetValue(i,j).r*pJ->GetValue(i,j).S[i2d_RoU];
        }
    }
    return(Mp);
}


double CalcXForceYSym2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                        double x0, // initial point of probed area
                        double l,  // length of probed area
                        double d   // diameter (or cross-section size) of probed area
                       ) {
    double Fp=0.; // pressure force
    double Fd=0.; // drag force
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP
#endif // _MPI_OPENMP

    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&
                 i >= (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i <= (int)((l+x0)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j <= (int)(d/FlowNode2D<double,NUM_COMPONENTS>::dy)) {

                double Sp = 0.;
                double Sd = 0.;

                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<double,NUM_COMPONENTS>::dy;
                    Sd = FlowNode2D<double,NUM_COMPONENTS>::dx;
                } else {
                    Sp = 2 * M_PI * pJ->GetValue(i,j).r*FlowNode2D<double,NUM_COMPONENTS>::dy;
                    Sd = 2 * M_PI * pJ->GetValue(i,j).r*FlowNode2D<double,NUM_COMPONENTS>::dx;
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

double CalcXForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // X initial point of probed area
                    double y0, // Y initial point of probed area
                    double dx, // X size of probed area
                    double dy  // Y of probed area
                   ) {
    double Fp=0.; // pressure force
    double Fd=0.; // drag force
#ifdef _MPI_OPENMP
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP
#endif // _MPI_OPENMP
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&
                 i >= (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<double,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<double,NUM_COMPONENTS>::dy)) {

                double Sp = 0.;
                double Sd = 0.;

                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<double,NUM_COMPONENTS>::dy;
                    Sd = FlowNode2D<double,NUM_COMPONENTS>::dx;
                } else {
                    Sp = 2 * M_PI * (j+0.5)*FlowNode2D<double,NUM_COMPONENTS>::dy*FlowNode2D<double,NUM_COMPONENTS>::dy;
                    Sd = 2 * M_PI * (j+0.5)*FlowNode2D<double,NUM_COMPONENTS>::dy*FlowNode2D<double,NUM_COMPONENTS>::dx;
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

//    cout << "\nPressure force (X):" << Fp << " N" << endl;
//    cout << "Drag force (X):" << Fd << " N" << endl;

    return(Fp+Fd);
}

double CalcYForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // X initial point of probed area
                    double y0, // Y initial point of probed area
                    double dx, // X size of probed area
                    double dy  // Y of probed area
                   ) {
    double Fp=0.; // pressure force
    double Fd=0.; // drag force
#ifdef _MPI_OPENMP 
#pragma omp parallel for reduction(+:Fp,Fd)
#else
#ifdef _OPEN_MP
#pragma omp parallel for reduction(+:Fp,Fd)
#endif //_OPEN_MP 
#endif // _MPI_OPENMP 
    for (int i=0;i<(int)pJ->GetX();i++ )
        for (int j=0;j<(int)pJ->GetY();j++ ) {
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&    
                 i >= (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<double,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<double,NUM_COMPONENTS>::dy)) {

                double Sp = 0.;
                double Sd = 0.;

                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Sp = FlowNode2D<double,NUM_COMPONENTS>::dx;
                    Sd = FlowNode2D<double,NUM_COMPONENTS>::dy;
                } else {
                    Sp = 2 * M_PI * pJ->GetValue(i,j).r*FlowNode2D<double,NUM_COMPONENTS>::dx;
                    Sd = 2 * M_PI * pJ->GetValue(i,j).r*FlowNode2D<double,NUM_COMPONENTS>::dy;
                }
                // Pressure forces

                if (j > 0 && pJ->GetValue(i,j-1).isCond2D(CT_SOLID_2D) ) { //  [gas]
                    Fp+= Sp*pJ->GetValue(i,j).p;                           // [solid] down force (+)
                } else if ( j < (int)pJ->GetY()-1  && pJ->GetValue(i,j+1).isCond2D(CT_SOLID_2D)) {
                    Fp-= Sp*pJ->GetValue(i,j).p;                           //  [gas]  up force (-)
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
//    cout << "\nPressure force (Y):" << Fp << " N" << endl;
//    cout << "Drag force (Y):" << Fd << " N" << endl;
    return(Fp+Fd);
}

double   Calc_Cp(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, int i, int j, Flow2D* pF) {
    if (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) || pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D))
        return(pJ->GetValue(i,j).p - pF->Pg())/(0.5*pF->ROG()*pF->Wg()*pF->Wg());
    else
        return 0;
}


double Calc_Cx_2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                  double x0, 
                  double y0,
                  double dx, 
                  double dy,
                  Flow2D* pF
                 ) {
    double Fmid=0,Pmax=0;
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
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D) ) &&
                 i >= (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<double,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<double,NUM_COMPONENTS>::dy)) {
                 is_node = 1;
            }
        }
        if (is_node) {
            if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT) {
                Fmid += FlowNode2D<double,NUM_COMPONENTS>::dy;
            } else {
                Fmid += 2 * M_PI * (j+0.5)*FlowNode2D<double,NUM_COMPONENTS>::dy*FlowNode2D<double,NUM_COMPONENTS>::dy;
            }
        }
    }
    Pmax = pF->ROG()*pF->U()*pF->U()*0.5*Fmid;
    if (Pmax == 0.)
        return 0;
    else
        return CalcXForce2D(pJ,x0,y0,dx,dy)/Pmax;
}

double   Calc_Cy_2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, 
                    double x0, 
                    double y0,
                    double dx, 
                    double dy,
                    Flow2D* pF) {
    double Fmid=0,Pmax=0;
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
            if ( (pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) ||
                  pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) &&
                 i >= (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 i <= (int)((x0+dx)/FlowNode2D<double,NUM_COMPONENTS>::dx) &&
                 j >= (int)(y0/FlowNode2D<double,NUM_COMPONENTS>::dy) &&
                 j <= (int)((y0+dy)/FlowNode2D<double,NUM_COMPONENTS>::dy)) {
                is_node = 1;
            }
            if (is_node) {
                if (FlowNode2D<double,NUM_COMPONENTS>::FT == FT_FLAT) {
                    Fmid += FlowNode2D<double,NUM_COMPONENTS>::dy;
                } else {
                    Fmid += 2 * M_PI * pJ->GetValue(0,j).r*FlowNode2D<double,NUM_COMPONENTS>::dy;
                }
            }
        }
    }

    Pmax = pF->ROG()*pF->Wg()*pF->Wg()*0.5*Fmid;
    if (Pmax == 0.)
        return 0;
    else
        return CalcYForce2D(pJ,x0,y0,dx,dy)/Pmax;

    return 0;
}


void SmoothY(UMatrix2D<double>* A) {
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

void SmoothX(UMatrix2D<double>* A) {
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
                     UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                     double  Ts) {
    char  HeatFluxHeader[128];
    snprintf(HeatFluxHeader,128,"#VARIABLES = X, HeatFlux(X)");
    *OutputData << HeatFluxHeader  << endl;
    UArray<double> HeatFluxArray(pJ->GetX());
    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        double Zero = 0.0;
        HeatFluxArray.SetElement(i,&Zero);
    }

    for(int i=0; i < (int)pJ->GetX(); i++) {
        for(int j=0; j < (int)pJ->GetY()-1; j++) {

        if(pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {

            FlowNode2D<double,NUM_COMPONENTS>* UpNode;
            FlowNode2D<double,NUM_COMPONENTS>* DownNode; 
            FlowNode2D<double,NUM_COMPONENTS>* RightNode;
            FlowNode2D<double,NUM_COMPONENTS>* LeftNode; 
            double lam_eff;
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
            double Q = lam_eff*(pJ->GetValue(i,j).Tg - Ts)/FlowNode2D<double,NUM_COMPONENTS>::dx;

            if(HeatFluxArray.GetElement(i) != 0.) {
              Q =  max(HeatFluxArray.GetElement(i),Q);
              HeatFluxArray.SetElement(i,&Q);
            } else {
              HeatFluxArray.SetElement(i,&Q);
            }
          }
        }
    }
    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        *OutputData << i*FlowNode2D<double,NUM_COMPONENTS>::dx << " " << HeatFluxArray.GetElement(i) << endl;
    }
}

void SaveYHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                     double  Ts) {
    char  HeatFluxHeader[128];
    snprintf(HeatFluxHeader,128,"#VARIABLES = Y, HeatFlux(Y)");
    *OutputData << HeatFluxHeader  << endl;
    UArray<double> HeatFluxArray(pJ->GetY());
    for(int i=0; i < (int)HeatFluxArray.GetNumElements(); i++) {
        double Zero = 0.0;
        HeatFluxArray.SetElement(i,&Zero); 
    }
    for(int j=0; j < (int)pJ->GetY(); j++) {
         for(int i=0; i < (int)pJ->GetX()-1; i++){
        if(pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {
            double lam_eff = pJ->GetValue(i,j).lam + pJ->GetValue(i,j).lam_t;
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

            FlowNode2D<double,NUM_COMPONENTS>* UpNode     = &(pJ->GetValue(i,N3));
            FlowNode2D<double,NUM_COMPONENTS>* DownNode   = &(pJ->GetValue(i,N4));
            FlowNode2D<double,NUM_COMPONENTS>* RightNode  = &(pJ->GetValue(N2,j));
            FlowNode2D<double,NUM_COMPONENTS>* LeftNode   = &(pJ->GetValue(N1,j));

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

            double Q = lam_eff*(pJ->GetValue(i,j).Tg - Ts)/FlowNode2D<double,NUM_COMPONENTS>::dy;
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
        *OutputData << j*FlowNode2D<double,NUM_COMPONENTS>::dy << " " << HeatFluxArray.GetElement(j) << endl;
    }
}


/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2014 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*************************************************/
#ifndef _flow_hpp
#define _flow_hpp

#include <math.h>
enum Func {
    TAUf,Pf,EPSf,Qf,Yf,Ff,Rf
};

enum FixedValue {
     FV_VELOCITY,
     FV_MACH
};

class Flow {

    double  r;                  
    double  t0;                 
    double  p0;                 
    double  Lambda;             
    double  k;                  

    double LMax;
    double LMin;
    double TestLam;
    int    iter;

    void MakeVar();
    void InitVar(double Other_k, double OtherT0,double OtherP0,double OtherR);
    void InitVar();

    inline double  Test(Func F, double Value);
    inline double  Test(Func F, double Value,int Area);
    inline double  TestFunc(Func F, double Value);

    double  TestFF(double AnyLambda);
    double  TestQF(double AnyLambda);
    double  TestTAU(double AnyLambda);
    double  TestPF(double AnyLambda) {
        return pow(TestTAU(AnyLambda),k/(k-1));
    }
    double  TestEPS(double AnyLambda) {
        return pow(TestTAU(AnyLambda),1/(k-1));
    }
    double  TestYF(double AnyLambda) {
        return TestQF(AnyLambda)/TestPF(AnyLambda);
    }
    double  TestZF(double AnyLambda) {
        return(AnyLambda + 1/AnyLambda);
    }
    double  TestRF(double AnyLambda) {
        return TestPF(AnyLambda)/TestFF(AnyLambda) ;
    }

public:

    double  C;
    double  lam;
    double  mu;

    double LMAX() {
        return sqrt((k+1)/(k-1));
    }
    double Tg() {
        return t0*TAU();
    }
    double EPS() {
        return TestEPS(Lambda);
    }
    double ZF(double NewZ, int Field=1);
    double Pr() {
        return C*mu/lam;
    }
    double ROG() {
        return EPS()*P0()/Rg()/T0();
    }
    double ROG(double newRo) {
        return EPS(newRo/P0()*Rg()*T0())*P0()/Rg()/T0();
    }
    double Pg() {
        return p0*PF();
    }
    double Pg(double NewPg);
    double Akr();
    double Asound();
    double MACH(double NewMach);
    double kg(double New_k);
    double TAU(double NewT);
    double PF(double NewPI);
    double EPS(double NewEPS);
    double YF(double NewY);
    double FF(double NewF, int Field= 1);
    double RF(double NewR);
    double YF();
    double FF();
    double RF();
    double MACH();
    double LAM();
    double LAM(double NewLambda);
    double kg();
    double Rg();
    double Rg(double NewR);
    double Tg(double NewT);
    double Wg();
    double Wg(double NewW);
    double T0();
    double T0(double NewT0);
    double P0();
    double P0(double NewP0);
    double TAU();
    double PF();
    double QF();
    double QF(double NewQ, int Field= 1);
    double ZF();
    double BF();
    double AF();

    void   CorrectFlow(double T, double p, double ref_val, FixedValue fv = FV_MACH );

    Flow();                                 
    Flow(double Other_Cp,double OtherT0, double OtherP0, double OtherR,  double Other_lam=0.01, double Other_mu=5.e-5 );
    Flow(Flow &OtherFlow);
    virtual ~Flow();                        

    Flow& operator = (Flow &NewFlow);

};

#endif // _flow_hpp






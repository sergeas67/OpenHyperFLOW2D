/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2015 by Serge A. Suchkov  *
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

    FP  r;                  
    FP  t0;                 
    FP  p0;                 
    FP  Lambda;             
    FP  k;                  

    FP LMax;
    FP LMin;
    FP TestLam;
    int    iter;

    void MakeVar();
    void InitVar(FP Other_k, FP OtherT0,FP OtherP0,FP OtherR);
    void InitVar();

    inline FP  Test(Func F, FP Value);
    inline FP  Test(Func F, FP Value,int Area);
    inline FP  TestFunc(Func F, FP Value);

    FP  TestFF(FP AnyLambda);
    FP  TestQF(FP AnyLambda);
    FP  TestTAU(FP AnyLambda);
    FP  TestPF(FP AnyLambda) {
        return pow(TestTAU(AnyLambda),k/(k-1));
    }
    FP  TestEPS(FP AnyLambda) {
        return pow(TestTAU(AnyLambda),1/(k-1));
    }
    FP  TestYF(FP AnyLambda) {
        return TestQF(AnyLambda)/TestPF(AnyLambda);
    }
    FP  TestZF(FP AnyLambda) {
        return(AnyLambda + 1/AnyLambda);
    }
    FP  TestRF(FP AnyLambda) {
        return TestPF(AnyLambda)/TestFF(AnyLambda) ;
    }

public:

    FP  C;
    FP  lam;
    FP  mu;

    FP LMAX() {
        return sqrt((k+1)/(k-1));
    }
    FP Tg() {
        return t0*TAU();
    }
    FP EPS() {
        return TestEPS(Lambda);
    }
    FP ZF(FP NewZ, int Field=1);
    FP Pr() {
        return C*mu/lam;
    }
    FP ROG() {
        return EPS()*P0()/Rg()/T0();
    }
    FP ROG(FP newRo) {
        return EPS(newRo/P0()*Rg()*T0())*P0()/Rg()/T0();
    }
    FP Pg() {
        return p0*PF();
    }
    FP Pg(FP NewPg);
    FP Akr();
    FP Asound();
    FP MACH(FP NewMach);
    FP kg(FP New_k);
    FP TAU(FP NewT);
    FP PF(FP NewPI);
    FP EPS(FP NewEPS);
    FP YF(FP NewY);
    FP FF(FP NewF, int Field= 1);
    FP RF(FP NewR);
    FP YF();
    FP FF();
    FP RF();
    FP MACH();
    FP LAM();
    FP LAM(FP NewLambda);
    FP kg();
    FP Rg();
    FP Rg(FP NewR);
    FP Tg(FP NewT);
    FP Wg();
    FP Wg(FP NewW);
    FP T0();
    FP T0(FP NewT0);
    FP P0();
    FP P0(FP NewP0);
    FP TAU();
    FP PF();
    FP QF();
    FP QF(FP NewQ, int Field= 1);
    FP ZF();
    FP BF();
    FP AF();

    void   CorrectFlow(FP T, FP p, FP ref_val, FixedValue fv = FV_MACH );

    Flow();                                 
    Flow(FP Other_Cp,FP OtherT0, FP OtherP0, FP OtherR,  FP Other_lam=0.01, FP Other_mu=5.e-5 );
    Flow(Flow &OtherFlow);
    virtual ~Flow();                        

    Flow& operator = (Flow &NewFlow);

};

#endif // _flow_hpp






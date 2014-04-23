/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2014 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*************************************************/
#include "libFlow/flow.hpp"

Flow::Flow () {
    Lambda = 0.01;
    InitVar();

    C = k*r/(k-1);
    lam = 0.01;
    mu  = 5.e-5;
}


Flow::Flow(double Other_Cp,double OtherT0, double OtherP0, double OtherR,  double Other_lam, double Other_mu ) {
    Lambda = 0.01;
    C = Other_Cp;
    InitVar(C/(C-OtherR), OtherT0, OtherP0, OtherR);
    lam = Other_lam;
    mu  = Other_mu;
}

Flow::Flow (Flow &OtherFlow) {
    Lambda = OtherFlow.Lambda;
    k      = OtherFlow.k;
    p0     = OtherFlow.p0;
    t0     = OtherFlow.t0;
    r      = OtherFlow.r;
    C      = OtherFlow.C;
    lam    = OtherFlow.lam;
    mu     = OtherFlow.mu;
}

Flow::~Flow () {
}

void   Flow::MakeVar() {
    return;
}


void   Flow::InitVar(double Other_k, double OtherT0, double OtherP0, double OtherR) {
    k  = Other_k;
    t0 = OtherT0;
    p0 = OtherP0;
    r  = OtherR;
    C  = k*r/(k-1);
}
void   Flow::InitVar() {
    InitVar(1.4, 300., 1.e5, 300.);
}

double Flow::TestTAU(double AnyLambda) {
    return(1-(k-1)/(k+1)*AnyLambda*AnyLambda);
}


double Flow::TestQF(double AnyLambda) {
    double rValue;
    rValue = pow((k+1)/2,1/(k-1))*AnyLambda;
    return rValue*pow(1-(k-1)/(k+1)*AnyLambda*AnyLambda,1/(k-1));
}


double Flow::LAM() {
    return Lambda;
}

double Flow::LAM(double NewLambda) {
    if(LMAX() > NewLambda && NewLambda > 0.) {
        Lambda = NewLambda;
        return Lambda;
    } else {
        return -1;
    }
}



double Flow::kg() {
    return k;
}
double Flow::kg(double New_k) {
    if(New_k > 0.) {
        k = New_k;
        return k;
    } else {
        return -1;
    }
}


double Flow::Rg() {
    return r;
}

double Flow::Rg(double NewR) {
    if(NewR > 0.) {
        r = NewR;
        return r;
    } else {
        return -1;
    }
}

double Flow::Tg(double NewT) {
    double TestT;
    if(t0 > NewT &&  NewT > 0.) {
        TestT = NewT/t0;
        Lambda = Test(TAUf, TestT);
        return Tg();
    } else {
        return -1;
    }
}


double Flow::BF() {

    return(sqrt(1.-1./k/k));

}

double Flow::AF() {

    return k*pow(2/(k+1),k/(k-1))*sqrt((k+1)/(k-1));

}

double Flow::QF() {
    return TestQF(Lambda);
}

double Flow::QF(double  NewQ, int Field) {
    if(Test(Qf, NewQ, Field) > 0.) {
        Lambda = Test(Qf, NewQ, Field);
        return TestQF(Lambda);
    } else {
        return -1;
    }
}

double Flow::T0() {
    return t0;
}

double Flow::T0(double NewT0) {
    if(NewT0 > 0.) {
        t0 = NewT0;
        return t0;
    } else {
        return t0;
    }
}

double Flow::P0() {
    return p0;
}

double Flow::P0(double NewP0) {
    if(NewP0 > 0.) {
        p0 = NewP0;
        return p0;
    } else {
        return p0;
    }
}

double Flow::Wg() {
    return Lambda*Akr();
}

double Flow::Wg(double NewW) {
    if(NewW > 0.) {
        if(NewW < Akr()*LMAX())
            Lambda = NewW/Akr();
        else return -1;

        return NewW;
    } else {
        return -1;
    }
}

double Flow::Akr() {
    return sqrt(2*k/(k+1)*r*t0);
}

double Flow::TAU() {
    return TestTAU(Lambda);
}


double Flow::Asound() {
    return sqrt(k*r*t0*TAU());
}

double Flow::TAU(double NewTAU) {
    if(1. > NewTAU && NewTAU > 0.) {
        Lambda = Test(TAUf, NewTAU);
        return TestTAU(Lambda);
    } else {
        return -1;
    }
}

double Flow::PF() {
    return TestPF(Lambda);
}

double Flow::PF(double NewPI) {
    if(1. > NewPI && NewPI > 0.) {
        Lambda = Test(Pf, NewPI);
        return TestPF(Lambda);
    } else {
        return -1;
    }
}

double Flow::ZF() {
    return TestZF(Lambda);
}

double Flow::ZF(double NewZ, int Field) {
    if(NewZ*NewZ < 4.) return -1.;
    if(Field < 0) {
        Lambda = (NewZ-sqrt(NewZ*NewZ-3.999999))/2;
    } else {
        Lambda = (NewZ+sqrt(NewZ*NewZ-3.999999))/2;
    }

    return TestZF(Lambda);
}
double Flow::MACH() {
    return Wg()/Asound();
}

double Flow::MACH(double NewMach) {
    if(NewMach < 0.) return -1;
    Lambda = sqrt((k+1)/2*NewMach*NewMach/(1+((k-1)/2*NewMach*NewMach)));
    return NewMach;
}

double Flow::EPS(double NewEPS)

{
    Lambda = Test(EPSf, NewEPS);
    return TestEPS(Lambda);
}

double Flow::YF() {
    return TestYF(Lambda);
}

double Flow::YF(double NewY)

{
    Lambda = Test(Yf, NewY);
    return TestYF(Lambda);
}

double Flow::TestFF(double AnyLambda) {
    return(AnyLambda*AnyLambda+1)*pow(1-(k-1)/(k+1)*
                                      AnyLambda*AnyLambda,1/(k-1));
}

double Flow::FF() {
    return TestFF(Lambda);
}

double Flow::FF(double  NewF, int Field)

{

    Lambda = Test(Ff, NewF, Field);
    return TestFF(Lambda);

}
double Flow::RF()

{
    return TestRF(Lambda);
}

double Flow::RF(double NewR)

{
    Lambda = Test(Rf, NewR);
    return TestRF(Lambda);
}

double Flow::Pg(double NewPg) {
    double TmpPI;
    if(NewPg >= P0()) return Pg();
    TmpPI = NewPg/P0();
    PF(TmpPI);
    return Pg();
}

double Flow::TestFunc(Func F, double Val)

{
    switch(F) {
    case TAUf:  return TestTAU(Val);
    case Pf:    return TestPF(Val);
    case EPSf:  return TestEPS(Val);
    case Qf:    return TestQF(Val);
    case Yf:    return TestYF(Val);
    case Ff:    return TestFF(Val);
    case Rf:    return TestRF(Val);
    default:    return -1.;
    }

}

double Flow::Test(Func F, double Val) {
    LMax = LMAX();
    LMin = 0.01;
    iter = 0;

    do {
        iter++;
        TestLam = (LMax + LMin)/2;
        if(TestFunc(F,TestLam) < Val) {
            LMax = TestLam;
        } else {
            LMin = TestLam;
        }

        if(iter > 100) return -1;
    }while(fabs((Val - TestFunc(F,TestLam))/Val)>0.01);

    return TestLam;
}

double Flow::Test(Func F, double Val , int Area) {
    if(Area < 0) {
        LMax = 0.01;
        LMin = 1.;
    } else {
        LMax = LMAX();;
        LMin = 1.;
    }
    iter = 0;

    do {
        iter++;
        TestLam = (LMax + LMin)/2;
        if(TestFunc(F, TestLam) < Val) {
            LMax = TestLam;
        } else {
            LMin = TestLam;
        }
        if(iter > 100) return -1;
    }while(fabs((Val - TestFunc(F, TestLam))/Val)>0.01);

    return TestLam;
}

Flow& Flow::operator = (Flow& NewFlow) {
    Lambda = NewFlow.Lambda;
    k      = NewFlow.k;
    p0     = NewFlow.p0;
    t0     = NewFlow.t0;
    r      = NewFlow.r;
    C      = NewFlow.C;
    lam    = NewFlow.lam;
    mu     = NewFlow.mu;

    return *this;
}

void   Flow::CorrectFlow(double T, double p, double ref_val, FixedValue fv) {

    if(fv == FV_MACH) {
       MACH(ref_val);
       t0 = T/TAU();
       p0 = p/PF();
    } else if(fv == FV_VELOCITY) {
       double res_p = 1.0, res_t = 1.0;
       int iter=0;
       do {
           MACH(ref_val/Asound());
           t0 = T/TAU();
           p0 = p/PF();
           res_p = fabs((p0-p/PF())/p0);
           res_t = fabs((t0-T/TAU())/t0);
           Wg(ref_val);
           iter++;
         } while ((res_p > 0.0001 || res_t > 0.0001) && iter < 1000);
    }
}
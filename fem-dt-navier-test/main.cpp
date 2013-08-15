#include <dolfin.h>
#include <boost/make_shared.hpp>
#include "navier_velocity.h"
#include "navier_pressure.h"
#include "navier_psi.h"
#include "navier_d1.h"
#include "navier_d2.h"
#include "navier_d3.h"
#include "navier_n1.h"
#include "navier_n2.h"
#include "navier_n3.h"
#include "navier_s1.h"
#include "navier_s2.h"
#include "navier_s3.h"
#include "navier_s4.h"
#include "navier_s5.h"


using namespace dolfin;

class DomainWalls : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return on_boundary && x[1] < 1 - DOLFIN_EPS;
  }
};

class DomainAllWalls : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return on_boundary;
  }
};

class DomainTop : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return on_boundary && x[1] > 1 - DOLFIN_EPS;
  }
};


int main(int argc, char *argv[])
{
    for(int i=0; i<argc; i++){
        info("param %d = %s", i, argv[i]);
    }
    if (argc < 6){
        if (MPI::process_number() == 0) {
            info("Usage: fem-dt-navier-test N DT T CONV RE");
        }
        exit(-1);
    }
    int N = 30;
    double DT = 0.5;
    double T = 10;
    std::string CONV = "d1";
    double REYNOLDS = 100;

    N = atoi(argv[1]);
    DT = atof(argv[2]);
    T = atof(argv[3]);
    CONV = argv[4];
    REYNOLDS = atof(argv[5]);

    if (MPI::process_number() == 0) {
        info("Using: %d %f %f %s %f", N, DT, T, CONV.c_str(), REYNOLDS);
    }

    UnitSquareMesh mesh(N, N);
    navier_velocity::FunctionSpace U(mesh);
    navier_pressure::FunctionSpace P(mesh);
    navier_psi::FunctionSpace PSI(mesh);
    Constant noslip(0, 0);
    Constant moveright(1, 0);
    Constant zero(0);
    DomainWalls walls;
    DomainTop top;
    DomainAllWalls all;
    DirichletBC bc0(U, noslip, walls);
    DirichletBC bc1(U, moveright, top);
    DirichletBC bc_psi(PSI, zero, all);
    std::vector<const BoundaryCondition*> bcs;
    bcs.push_back(&bc0);
    bcs.push_back(&bc1);

    Function u1(U);
    Function u0(U);
    Function u12(U);
    Function p0(P);
    Function p1(P);
    Function psi1(PSI);
    Constant tau(DT);
    Constant Re(REYNOLDS);

    boost::shared_ptr<Form> bf_pred;
    boost::shared_ptr<Form> lf_pred;
    if (CONV == "d1"){
        bf_pred = boost::shared_ptr<Form>(new navier_d1::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_d1::LinearForm(U));
    } else if (CONV == "d2") {
        bf_pred = boost::shared_ptr<Form>(new navier_d2::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_d2::LinearForm(U));
    } else if (CONV == "d3") {
        bf_pred = boost::shared_ptr<Form>(new navier_d3::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_d3::LinearForm(U));
    } else if (CONV == "n1") {
        bf_pred = boost::shared_ptr<Form>(new navier_n1::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_n1::LinearForm(U));
    } else if (CONV == "n2") {
        bf_pred = boost::shared_ptr<Form>(new navier_n2::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_n2::LinearForm(U));
    } else if (CONV == "n3") {
        bf_pred = boost::shared_ptr<Form>(new navier_n3::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_n3::LinearForm(U));
    } else if (CONV == "s1") {
        bf_pred = boost::shared_ptr<Form>(new navier_s1::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_s1::LinearForm(U));
    } else if (CONV == "s2") {
        bf_pred = boost::shared_ptr<Form>(new navier_s2::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_s2::LinearForm(U));
    } else if (CONV == "s3") {
        bf_pred = boost::shared_ptr<Form>(new navier_s3::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_s3::LinearForm(U));
    } else if (CONV == "s4") {
        bf_pred = boost::shared_ptr<Form>(new navier_s4::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_s4::LinearForm(U));
    } else if (CONV == "s5") {
        bf_pred = boost::shared_ptr<Form>(new navier_s5::BilinearForm(U, U));
        lf_pred = boost::shared_ptr<Form>(new navier_s5::LinearForm(U));
    } else {
        info("Unknown CONV: %s", CONV.c_str());
        exit(-1);
    }
    bf_pred->set_coefficient("tau", reference_to_no_delete_pointer(tau));
    bf_pred->set_coefficient("Re", reference_to_no_delete_pointer(Re));
    bf_pred->set_coefficient("u0", reference_to_no_delete_pointer(u0));
    lf_pred->set_coefficient("tau", reference_to_no_delete_pointer(tau));
    lf_pred->set_coefficient("u0", reference_to_no_delete_pointer(u0));
    lf_pred->set_coefficient("p0", reference_to_no_delete_pointer(p0));

    navier_pressure::BilinearForm bf_pressure(P, P);
    navier_pressure::LinearForm lf_pressure(P);
    lf_pressure.tau = tau;
    lf_pressure.u12 = u12;
    lf_pressure.p0 = p0;

    navier_velocity::BilinearForm bf_velocity(U, U);
    navier_velocity::LinearForm lf_velocity(U);
    bf_velocity.tau = tau;
    lf_velocity.tau = tau;
    lf_velocity.u12 = u12;
    lf_velocity.p0 = p0;
    lf_velocity.p1 = p1;

    navier_psi::BilinearForm bf_psi(PSI, PSI);
    navier_psi::LinearForm lf_psi(PSI);
    lf_psi.u1 = u1;

    char s[1024];
    sprintf(s, "result_%s_%d_%f_%f/velocity.pvd", CONV.c_str(), N, DT, REYNOLDS);
    File fvelo(s);
    sprintf(s, "result_%s_%d_%f_%f/pressure.pvd", CONV.c_str(), N, DT, REYNOLDS);
    File fpress(s);
    sprintf(s, "result_%s_%d_%f_%f/psi.pvd", CONV.c_str(), N, DT, REYNOLDS);
    File fpsi(s);
    double t = 0;
    while(t <= T) {
        if (MPI::process_number() == 0){
            info("t = %lf", t);
        }
        solve(*bf_pred == *lf_pred, u12, bcs);
        solve(bf_pressure == lf_pressure, p1);
        solve(bf_velocity == lf_velocity, u1, bcs);
        solve(bf_psi == lf_psi, psi1, bc_psi);

        double psi_max = psi1.vector()->max();
        double psi_min = psi1.vector()->min();
        psi_max = MPI::max(psi_max);
        psi_min = MPI::min(psi_min);
        if (MPI::process_number() == 0){
            info("%lf %lf %lf", t, psi_max, psi_min);
        }
        fvelo << u1;
        fpress << p1;
        fpsi << psi1;
        u0 = u1;
        p0 = p1;
        t += DT;
    }
    return 0;
}

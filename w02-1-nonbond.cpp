#include <iostream>
#include <cmath>
using namespace std;

class ATOM {
 public:
  double x[3];  // position vector of atom in Angstroms
  double chg;   // partial charges in electrons
  double R0;    // LJ (Lennard-Jones) R0 parameter in Angstroms
  double D0;    // LJ D0 parameter in kcal/mol

  ATOM() { x[0] = x[1] = x[2] = chg = R0 = D0 = 0; };
};

class NONBOND {
 public:
  double r;                        // separation distance
  double r2;                       // r*r
  double vdw;                      // van der waals energy in kcal/mol
  double cou;                      // coulomb energy in kcal/mol
  double cal_r(ATOM &, ATOM &);    // calculate r and r2
  double cal_eng(ATOM &, ATOM &);  // calculate nonbond energy

  NONBOND() { r = r2 = vdw = cou = 0; };
};

int main(int argc, char **argv) {
  ATOM a[2];
  
  // data for Oxygen
  a[0].x[0] = 0;
  a[0].x[1] = 0;
  a[0].x[2] = 0;
  a[0].chg = -0.82;
  a[0].R0 = 3.5532;
  a[0].D0 = 0.1848;

  // data for Hydrogen
  a[1].x[0] = 0;
  a[1].x[1] = 0.7;
  a[1].x[2] = -0.7;
  a[1].chg = 0.41;
  a[1].R0 = 0.9;
  a[1].D0 = 0.01;

  NONBOND nb;
  nb.cal_eng(a[0], a[1]);
  cout << "vdw: " << nb.vdw << " kcal/mol, cou: " << nb.cou << " kcal/mol"
       << endl;

  return 0;
}

double NONBOND::cal_r(ATOM &a, ATOM &b) {
  r2 = (a.x[0] - b.x[0]) * (a.x[0] - b.x[0]) +
       (a.x[1] - b.x[1]) * (a.x[1] - b.x[1]) +
       (a.x[2] - b.x[2]) * (a.x[2] - b.x[2]);
  r = sqrt(r2);
  return r;
}

double NONBOND::cal_eng(ATOM &a, ATOM &b) {
  double R0, D0;

  R0 = 0.5 * (a.R0 + b.R0);
  D0 = sqrt(a.D0 * b.D0);
  cout << "R0: " << R0 << endl;
  cout << "D0: " << D0 << endl;

  vdw = D0 * (pow(R0 / cal_r(a, b), 12) - 2 * pow(R0 / cal_r(a, b), 6));
  cou = 332.0637 * a.chg * b.chg / cal_r(a, b);
  return 0;
}
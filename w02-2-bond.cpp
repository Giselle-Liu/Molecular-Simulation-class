#include <iostream>
#include <cmath>
using namespace std;


class ATOM {
 public:
  double x[3];  // position vector of atom in Angstroms
  double chg;   // partial charges in electrons
  double R0;    // LJ R0 parameter in Angstroms
  double D0;    // LJ D0 parameter in kcal/mol

  ATOM() { x[0] = x[1] = x[2] = chg = R0 = D0 = 0; };
};

class BOND {
 public:
  ATOM* atm[2];
  double len;        // in Angstroms
  double len0;       // equilib bond length in Angstroms
  double Kb;         // force constant in kcal/mol A2
  double eng;        // bond energy in kcal/mol
  double cal_len();  // calculate bond length from two vectors
  double cal_eng();  // calculate bond energy

  BOND() { len = len0 = Kb = eng = 0; };
};

int main(int argc, char** argv) {
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

  BOND bond;
  bond.atm[0] = &a[0];
  bond.atm[1] = &a[1];
  bond.Kb = 500;
  bond.len0 = 1;
  bond.cal_eng();
  cout << "bond length is " << bond.len << " A, energy is " << bond.eng
       << " kcal/mol" << endl;

  return 0;
}

double BOND::cal_len() {
  double r2;
  r2 = (atm[0]->x[0] - atm[1]->x[0]) * (atm[0]->x[0] - atm[1]->x[0]) +
       (atm[0]->x[1] - atm[1]->x[1]) * (atm[0]->x[1] - atm[1]->x[1]) +
       (atm[0]->x[2] - atm[1]->x[2]) * (atm[0]->x[2] - atm[1]->x[2]);

  len = sqrt(r2);

  return 0;
}

double BOND::cal_eng() {
  cal_len();
  //	eng=0.5*Kb*pow((cal_len()-len0),2);
  eng = 0.5 * Kb * pow((len - len0), 2);

  return 0;
}
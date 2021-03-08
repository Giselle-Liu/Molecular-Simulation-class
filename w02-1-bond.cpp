#include <iostream>
#include <cmath>
using namespace std;

/* run this program using the console pauser or add your own getch,
 * system("pause") or input loop */

class ATOM {
 public:
  ATOM() { x[0] = x[1] = x[2] = chg = Ro = Do = 0; };

 public:
  double x[3];  // position vector of atom in angstroms
  double chg;   // partial charges in electrons
  double Ro;    // LJ Ro parameter in Angstroms
  double Do;    // LJ Do parameter in kcal/mol
};

class BOND {
 public:
  BOND() { len = len0 = Kb = eng = 0; };
  ATOM* atm[2];
  double len;        // in angstroms
  double len0;       // equilib bond length in angstroms
  double Kb;         // force constant in kcal/mol A2
  double eng;        // bond energy in kcal/mol
  double cal_len();  // calculate bond length from two vectors
  double cal_eng();  // calculate bond energy
};

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

int main(int argc, char** argv) {
  ATOM a[2];
  // data for Oxygen
  a[0].x[0] = 0;
  a[0].x[1] = 0;
  a[0].x[2] = 0;
  a[0].chg = -0.82;
  a[0].Ro = 3.5532;
  a[0].Do = 0.1848;

  // data for Hydrogen
  a[1].x[0] = 0;
  a[1].x[1] = 0.7;
  a[1].x[2] = -0.7;
  a[1].chg = 0.41;
  a[1].Ro = 0.9;
  a[1].Do = 0.01;

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

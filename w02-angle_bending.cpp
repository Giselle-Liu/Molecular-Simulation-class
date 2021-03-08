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

class ANGLE {
 public:
  ATOM* atm[3];
  ANGLE() { ang0 = K = eng = 0; };
  double ang0;
  double ang;
  double K;
  double eng;
  double cal_angle();
  double cal_eng();
};

double ANGLE::cal_angle() {
  double costheta, b12, b22, ip;
  double p[2][2];

  p[0][0] = atm[1]->x[0] - atm[0]->x[0];
  p[0][1] = atm[1]->x[1] - atm[0]->x[1];
  p[0][2] = atm[1]->x[2] - atm[0]->x[2];

  p[1][0] = atm[2]->x[0] - atm[0]->x[0];
  p[1][1] = atm[2]->x[1] - atm[0]->x[1];
  p[1][2] = atm[2]->x[2] - atm[0]->x[2];

  b12 = p[0][0] * p[0][0] + p[0][1] * p[0][1] + p[0][2] * p[0][2];
  b22 = p[1][0] * p[1][0] + p[1][1] * p[1][1] + p[1][2] * p[1][2];
  ip = p[0][0] * p[1][0] + p[0][1] * p[1][1] + p[0][2] * p[1][2];

  costheta = ip / (b12 * b22);

  return acos(costheta);
};
double ANGLE::cal_eng() {
  double ang = cal_angle();
  eng = 0.5 * K / pow(sin(ang0 * acos(-1) / 180), 2) *
        pow(cos(ang) - cos(ang0 * acos(-1) / 180), 2);

  return 0;
}

int main(int argc, char** argv) {
  ATOM a[3];
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

  // data for Hydrogen
  a[2].x[0] = 1;
  a[2].x[1] = 0;
  a[2].x[2] = 0;
  a[2].chg = 0.41;
  a[2].Ro = 0.9;
  a[2].Do = 0.01;

  ANGLE ang;
  ang.atm[0] = &a[0];
  ang.atm[1] = &a[1];
  ang.atm[2] = &a[2];

  ang.ang0 = 109.5;
  ang.K = 100;

  ang.cal_eng();

  cout << "bond angle is " << acos(ang.ang) / acos(-1) * 180 << " A, energy is "
       << ang.eng << " kcal/mol" << endl;

  return 0;
}

#include <star.h>
#include <stdlib.h>
#include <iostream>

double M(star2d &s, unsigned int ell);
double Mcgs(star2d &s, unsigned int ell);
double Mgeom(star2d &s, unsigned int ell);

double S(star2d &s, unsigned int ell);
double Scgs(star2d &s, unsigned int ell);
double Sgeom(star2d &s, unsigned int ell);

double Iint(star2d &s);
double Icgs(star2d &s);
double Igeom(star2d &s);

double chi(star2d &s);

int main(int argc,char *argv[]) {

  if(argc!=2) {
    printf("Usage: %s model_file \n", argv[0]);
    return 1;
  };

  star2d s;
  s.read(argv[1]);

  const double chis = chi(s);

  const double M0g = Mgeom(s,0), S1g = Sgeom(s,1);

  for ( int ell = 0; ell < s.nth && ell <= 10; ell++)
    if (0 == ell % 2) { // even
      const double Mlg = Mgeom(s, ell);
      const double Mlbar = Mlg / M0g / pow(M0g * chis , ell);
      // printf("M_{%d} = %.6e cm^%d\n", ell, Mlg, ell+1);
      // printf("\\bar{M}_{%d} = %.6e\n", ell, Mlbar);
      printf("%d,%.6e,%.6e\n", ell, Mlg, Mlbar);
    } else { // odd
      const double Slg = Sgeom(s,ell);
      const double Slbar = Slg / M0g / pow(M0g * chis , ell);
      // printf("S_{%d} = %.6e cm^%d\n", ell, Slg, ell+1);
      // printf("\\bar{S}_{%d} = %.6e\n", ell, Slbar);
      printf("%d,%.6e,%.6e\n", ell, Slg, Slbar);
    };

  // I computed from integral
  double Ival = Igeom(s);
  // Equatorial omega
  double we = s.map.leg.eval_00(s.w.row(-1),PI/2)(0);
  printf("%.6f,%.6f,%.6e\n",
         Ival,
         chis,
         we * s.units.Omega / C_LIGHT); // in rad/cm

  return 0;

};

double M(star2d &s, unsigned int ell)
{
  if (1 == ell%2)
    return 0;

  matrix Pl = s.map.leg.P1_00.row(ell/2) / sqrt(2*ell+1);
  double Ml = 2.0*PI * (s.map.gl.I,
                      s.rho * Pl * pow(s.r,ell+2.0) * s.map.rz,
                      s.map.leg.I_00)(0);

  return Ml;
};

double Mcgs(star2d &s, unsigned int ell)
{

  double Ml = M(s, ell);
  Ml *= s.units.rho * pow(s.units.r, ell+3.0);

  return Ml;
};

double Mgeom(star2d &s, unsigned int ell)
{

  double Ml = Mcgs(s, ell);
  Ml *= GRAV / (C_LIGHT * C_LIGHT);

  return Ml;
};

double S(star2d &s, unsigned int ell)
{
  if (0 == ell%2)
    return 0;

  matrix dPl = s.map.leg.P1_10.row((ell-1)/2) / sqrt(2.0*ell+1.0);
  double Sl = - 4.0*PI / (ell + 1.0) *
               (s.map.gl.I,
                s.rho * s.w * sin(s.th)
                * dPl * pow(s.r,ell+3.0) * s.map.rz,
                s.map.leg.I_00)(0);

  return Sl;
};

double Scgs(star2d &s, unsigned int ell)
{

  double Sl = S(s, ell);
  Sl *= s.units.rho * s.units.Omega * pow(s.units.r, ell+4.0);

  return Sl;
};

double Sgeom(star2d &s, unsigned int ell)
{

  double Sl = Scgs(s, ell);
  Sl *= GRAV / pow(C_LIGHT, 3.0);

  return Sl;
};

double Iint(star2d &s)
{
   return 4./3.*PI * (s.map.gl.I,
                      s.rho * pow(s.r,4.0) * s.map.rz,
                      s.map.leg.I_00)(0);
};

double Icgs(star2d &s)
{
  double Ival = Iint(s);
  Ival *= s.units.rho * pow(s.units.r, 5.0);

  return Ival;
};

double Igeom(star2d &s)
{
  double Ival = Icgs(s);
  Ival *= GRAV / (C_LIGHT * C_LIGHT);

  return Ival;
};

double chi(star2d &s)
{
  const double M0 = M(s,0);
  const double S1 = S(s,1);

  return C_LIGHT * s.units.Omega * S1 /
    (GRAV * s.units.rho * s.units.r * M0*M0);
};

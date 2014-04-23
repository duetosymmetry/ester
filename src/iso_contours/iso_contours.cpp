#include <star.h>
#include <stdlib.h>
#include <iostream>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

struct contour_root_params
{
  star2d *s;
  matrix *field_theta;
  double level;
};

double contour_gsl_f( double r, void * vparams);

int main(int argc,char *argv[]) {

  if(argc!=3) {
    printf("Usage: %s model_file ncontours \n", argv[0]);
    return 1;
  };

  star2d s;
  s.read(argv[1]);

  matrix &field = s.rho;

  const int ncontours = atoi(argv[2]);
  if (ncontours < 1 || ncontours >= s.nr) {
    printf("ERROR: Choose saner number of contours.\n");
    return -1;
  };

  // The last contour is simply the surface.
  // Do that one separately.
  matrix equatorial_zeta = vector(s.z(s.nr-1, 0) / ncontours,
                                   s.z(s.nr-1, 0) * (1. - 1. / ncontours),
                                   ncontours-1);

  matrix levels = s.map.eval_z( field, equatorial_zeta, zeros( 1, ncontours-1 ) );

  matrix field_theta, r_theta;

  matrix contours(ncontours, s.nth);

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *sol;
  T = gsl_root_fsolver_brent;
  sol = gsl_root_fsolver_alloc (T);

  contour_root_params p;
  p.s = &s;
  p.field_theta = &field_theta;
  gsl_function F;
  F.function = &contour_gsl_f;
  F.params   = &p;

  for (int ith = 0; ith < s.nth; ith++) {
    field_theta = s.map.leg.eval_00(field, s.map.th(ith));
    r_theta     = s.map.leg.eval_00(s.map.r, s.map.th(ith));

    for (int icontour = 0; icontour < ncontours-1; icontour++) {
      p.level = levels(icontour);

      gsl_root_fsolver_set (sol, &F, s.z(0), s.z(s.nr-1) );

      int iter = 0, max_iter=100;
      int status;
      double zeta;
      do {
        iter++;
        status = gsl_root_fsolver_iterate (sol);
        zeta = gsl_root_fsolver_root (sol);
        double z_lo = gsl_root_fsolver_x_lower (sol);
        double z_hi = gsl_root_fsolver_x_upper (sol);
        status = gsl_root_test_interval (z_lo, z_hi,
                                         0, 1.e-6);

      } while (status == GSL_CONTINUE && iter < max_iter);

      if (status != GSL_SUCCESS)
        printf ("Warning: didn't converge!\n");

      // printf(" contours( %d, %d ) = %.4f\n", icontour, ith, zeta);
      contours( icontour, ith ) = s.map.gl.eval( r_theta, zeta)(0);

    };
  };

  // Surface contour
  for (int ith = 0; ith < s.nth; ith++) {
    contours( ncontours - 1, ith) = s.map.R(s.ndomains, ith);
  };

  // Output
  // Can't just use .write() because I want to output a rectangular csv file
  for (int ith = 0; ith < s.nth; ith++) {
    printf("%.14e",s.th(ith));
    for (int icontour = 0; icontour < ncontours; icontour++)
      printf(",%.14e",contours(icontour, ith));
    printf("\n");
  };

  // s.th.write();
  // contours.write();

  gsl_root_fsolver_free (sol);

  return 0;

};

double contour_gsl_f( double zeta, void * vparams)
{
  contour_root_params * p = (contour_root_params *)vparams;
  const double field_val = p->s->map.gl.eval( *(p->field_theta), zeta )(0);
  //printf("field(zeta=%.3f,th=%.3f) = %.4f (want %.4f)\n", zeta, p->theta(0), field_val, p->level );
  return field_val - p->level;
};

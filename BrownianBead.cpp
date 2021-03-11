/*
*    Plug-in program for plot treatement in Xvin.
 *
 *    V. Croquette
  */
#ifndef _BROWNIANBEAD_C_
#define _BROWNIANBEAD_C_

#include <allegro.h>
# include "xvin.h"

/* If you include other regular header do it here*/

#include "HydroCorrections.hh"
#include "BrSimulator.hh"


/* But not below this define */
# define BUILDING_PLUGINS_DLL
 # include "BrownianBead.hh"


int do_BrownianBead_simulate_Brownianmotion(void)
{
  O_p *op = NULL;
  d_s *dsx = NULL, *dsy = NULL, *dsz=NULL;
  pltreg *pr = NULL;
  int i = 0;
  static double Ra = 500;
  static double time = 10;
  static double dt = 1;
  static double force = 13;
  static double angle = M_PI/2;
  static double anisotropy = 0.005;
  static double dsdna = 120;
  static int two_D = 0;
  static int no_coupling = 0;
  int n = 0;

  float *x = NULL, *y = NULL, *z=NULL;

  if(updating_menu_state != 0)	return D_O_K;


  if (ac_grep(cur_ac_reg,"%pr",&pr) != 1)
    return win_printf_OK("cannot find data");



  i = win_scanf("Indicate the radius of the bead %lf (nm) \n"
  "Indicate the simulation time %lf (s) \n"
  "Indicate the simulation time step %lf (us) \n"
  "Indicate the magnetic force %lf (pN)\n"
  "Indicate the angle between anchoring point and magnetic pole %lf (rad)\n"
  "Indicate the magnetic anisotropy %lf (percent)\n"
  "Indicate the number of double stranded bases %lf \n"
  "Keep in 2D (no x translation and rotation only around x) %b \n"
  "Remove translation/rotation coupling %b\n",&Ra,&time,&dt,&force,&angle,&anisotropy,&dsdna,&two_D,&no_coupling);
  if (i == WIN_CANCEL)	return OFF;

  dt = dt * 1e-6;

  BrSimulator B(Ra,force,angle,time,dt,anisotropy,dsdna,1-two_D,no_coupling);
  while (B.iterate() == 0)
  {
  ;
  }
  n = B.export_simulation(&x,&y,&z);

  op = create_and_attach_one_plot(pr,16,16,0);
  dsx = op->dat[0];
  float *t;
  t = static_cast<float*>(calloc(n,sizeof(float)));
  for (i=0;i<n;i++)
  {
    t[i] = dt*i;
  }
  dsx->xd=t;
  dsx->yd = x;
  dsx->nx = dsx->mx = dsx->ny = dsx->my = n;
  dsy = create_and_attach_one_ds(op, 16, 16, 0);
  float *ty;
  ty = static_cast<float*>(calloc(n,sizeof(float)));
  for (i=0;i<n;i++)
  {
    ty[i] = dt*i;
  }
  dsy->xd=ty;
  dsy->yd=y;
  dsy->nx = dsy->mx =   dsy->ny = dsy->my = n;
  dsz = create_and_attach_one_ds(op, 16, 16, 0);
  float *tz;
  tz = static_cast<float*>(calloc(n,sizeof(float)));

  for (i=0;i<n;i++)
  {
    tz[i] = dt*i;
  }
  dsz->nx = dsz->mx =   dsz->ny = dsz->my = n;
  dsz->xd = tz;
  dsz->yd = z;
  dsz->nx = dsz->mx = n;

  dt = dt * 1e6;


  create_attach_select_y_un_to_op(op, IS_METER, 0 , (float)0.001, -6, 0, "\\mu m");
  create_attach_select_x_un_to_op(op, IS_SECOND, 0 , (float)1, 0, 0, "s");

  refresh_plot(pr, UNCHANGED);
  return D_O_K;
}



MENU *BrownianBead_plot_menu(void)
{
  static MENU mn[32];

  if (mn[0].text != NULL)	return mn;
  add_item_to_menu(mn,"Simulate Brownian Bead", do_BrownianBead_simulate_Brownianmotion,NULL,0,NULL);

  return mn;
}

int	BrownianBead_main(int argc, char **argv)
{
  (void)argc;  (void)argv;  add_plot_treat_menu_item ( "BrownianBead", NULL, BrownianBead_plot_menu(), 0, NULL);
  return D_O_K;
}

int	BrownianBead_unload(int argc, char **argv)
{
  (void)argc;  (void)argv;  remove_item_to_menu(plot_treat_menu, "BrownianBead", NULL, NULL);
  return D_O_K;
}
#endif

/*
Fuchsian discrete hyperbolic geometry example for thesis example in chapter 3.

compile: gcc fuchsian.c -lGL -lglut -lGLU -lm -o  fuch
./fuch

toorandom@gmai.com

*/

#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#else

#include <GL/glut.h>
#include <math.h>
#endif


#define MAX_W 5001


double maxnormdelta = 0.2;
int firsttime = 1;
complex ga, gb, gc, gd;

int globalrand;
typedef struct mat2_s
{
  complex a, b, c, d;
} mat2;


typedef struct gens_s
{
  mat2 A, B, Ai, Bi, AB;
} gens;

typedef struct word_s
{
  unsigned char l[MAX_W];
  int siz;
} word;


mat2
mult (mat2 M, mat2 N)
{
  mat2 R;
  R.a = M.a * N.a + M.b * N.c;
  R.b = M.a * N.b + M.b * N.d;
  R.c = M.c * N.a + M.d * N.c;
  R.d = M.c * N.b + M.d * N.d;
  return R;
}



double
determinant (mat2 M)
{

  return M.a * M.d - M.b * M.c;

}


void
get_circle_center (double x1, double y1, double x2, double y2, double x3,
		   double y3, double *center)
{
  double mr = (y2 - y1) / (x2 - x1);
  double mt = (y3 - y2) / (x3 - x2);
  double X, Y;

  X = (mr * mt * (y3 - y1) + mr * (x2 + x3) - mt * (x1 + x2)) / (2 * (mr - mt));
  Y = (-1 / mr) * (X - ((x1 + x2) / 2)) + ((y1 + y2) / 2);
  center[0] = X;
  center[1] = Y;
  return;
}



complex
mobius_eval (mat2 M, complex z)
{
  return (M.a * z + M.b) / (M.c * z + M.d);

}

void
print_matrix (mat2 M)
{

  printf ("%.02f+i%.02f\t%.02f+i%.02f\n", creal (M.a), cimag (M.a),
	  creal (M.b), cimag (M.b));
  printf ("%.02f+i%.02f\t%.02f+i%.02f\n\n", creal (M.c), cimag (M.c),
	  creal (M.d), cimag (M.d));
  return;

}

mat2
gen_kl_elem (gens G, word W)
{

  mat2 Gens[5];
  mat2 R, T1, T2, Id;
  int i;

  Id.a = 1;
  Id.b = 0;
  Id.c = 0;
  Id.d = 1;

  T1 = Id;

  Gens[0] = Id;
  Gens[1] = G.A;
  Gens[2] = G.B;
  Gens[3] = G.Ai;
  Gens[4] = G.Bi;


  for (i = 0; i < W.siz; i++)
    {
      // print_matrix(Gens[W.l[i]]); 
      T2 = mult (T1, Gens[W.l[i]]);
      T1 = T2;
    }

  R = T1;
  return R;

}

word
random_n_word (int n)
{

  int i;
  word W;
  W.siz = n;
  for (i = 0; i < n; i++)
    {
      W.l[i] = rand () % 5;
//printf("%d ",W.l[i]);
    }
  W.l[n] = 0xff;
//printf("\n");
  return W;
}




gens
get_parabolic_comm (complex Ta, complex Tb)
{

  gens G;
  complex z0, Tab;
  complex detA, detB;

  Tab = (Ta * Tb + csqrt (Ta * Tb * Ta * Tb - 4 * (Ta * Ta + Tb * Tb))) / 2.0;

  z0 = ((Tab - 2.0) * Tb) / (Tb * Tab - 2.0 * Ta + 2.0 * I * Tab);

  G.A.a = Ta / 2.0;
  G.A.b = (Ta * Tab - 2.0 * Tb + 4.0 * I) / ((2.0 * Tab + 4) * z0);
  G.A.c = ((Ta * Tab - 2.0 * Tb - 4 * I) * z0) / (2.0 * Tab - 4.0);
  G.A.d = Ta / 2.0;

  G.B.a = (Tb - 2.0 * I) / 2.0;
  G.B.b = Tb / 2.0;
  G.B.c = Tb / 2.0;
  G.B.d = (Tb + 2.0 * I) / 2.0;

  G.AB.a = Tab / 2.0;
  G.AB.b = (Tab - 2.0) / (2.0 * z0);
  G.AB.c = ((Tab + 2.0) * z0) / 2.0;
  G.AB.d = Tab / 2.0;

  detA = G.A.a * G.A.d - G.A.b * G.A.c;
  detB = G.B.a * G.B.d - G.B.b * G.B.c;

  G.Ai.a = G.A.d / detA;
  G.Ai.b = -G.A.b / detA;
  G.Ai.c = -G.A.c / detA;
  G.Ai.d = G.A.a / detA;

  G.Bi.a = G.B.d / detB;
  G.Bi.b = -G.B.b / detB;
  G.Bi.c = -G.B.c / detB;
  G.Bi.d = G.B.a / detB;

  return G;


}

void
mobius_fixed_points (mat2 M, complex * fix)
{


  if (M.c != 0)
    {
      fix[0] =
	((M.a - M.d) +
	 csqrt ((M.a + M.d) * (M.a + M.d) -
		4 * (M.a * M.d - M.b * M.c))) / 2 * M.c;
      fix[1] =
	((M.a - M.d) -
	 csqrt ((M.a + M.d) * (M.a + M.d) -
		4 * (M.a * M.d - M.b * M.c))) / 2 * M.c;
    }

  if ((M.c == 0) && (M.a == M.d))
    {
      fix[0] = 0;
      fix[1] = 0;
    }


  if (M.c == 0)
    {
      fix[0] = -M.b / (M.a - M.d);
      fix[1] = -M.b / (M.a - M.d);
    }



}

void plot_circle(double *c, double r) {

double t;
double ss = 100;

if(r > 0.01)
return; 


glBegin(GL_LINE_LOOP);
for(t=0;t<2*M_PI;t += (2*M_PI)/ss) {

glVertex2f(c[0]+r*cos(t), c[1]+r*sin(t));


}

glEnd();


}

void
kleinian_action_s1 (complex ta, complex tb)
{

  complex x, z, w;
  double t, x1, y1, x2, y2, x3, y3, center[2];
  int i, j;
  gens G;
  mat2 fuchsianG[100000];

  //ta = -(1.87+I*0.1);
  // tb = -(1.87-I*0.1);

  //ta = -2;
  //tb = -2;


//  get_circle_center(cos(M_PI/5),sin(M_PI/5),cos(M_PI/7),sin(M_PI/7),cos(M_PI/4),sin(M_PI/4),(double *)&center);
//printf("%f %f\n",center[0],center[1]);
//exit(1);


  G = get_parabolic_comm (ta, tb);

//printf("NORM Cdelta = %f\n",fabs(cabs(G.A.c) - cabs(gc)));



  if (!firsttime)
    {
      // printf("Checking VALUE!\n");
      if (fabs (cabs (G.A.c) - cabs (gc)) > maxnormdelta)
	{
	  // glClearColor (1, 1, 1, 1);
	  // glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  // glutSwapBuffers ();
	  printf ("Sleeping, no draw will be done for:\n");
	  print_matrix (G.A);
	  //sleep(1);
	  return;
	}
    }
  else
    {
      firsttime = 0;
    }
//printf("det B=%f\nB =:\n",determinant(G.B));
//  print_matrix (G.B);

  glClearColor (1, 1, 1, 1);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glColor3f (0, 0, 0);
//printf("Plotting for:\n");
//print_matrix(G.A);
//  printf ("Tr z = %f+I%f\t Tr w = %f + I%f\n", creal (ta), cimag (ta),
	//  creal (tb), cimag (tb));
  //glBegin (GL_POINTS);

glBegin(GL_POINTS);

  for (i = 0; i < 50000; i++)
    {
      fuchsianG[i] = gen_kl_elem (G, random_n_word (i % 20));


/*
      for (t = 0; t < 2 * M_PI; t += 0.05)
	{
	  z = cos (t) + I * sin (t);
	  w = mobius_eval (fuchsianG[i], z);
	  glVertex2f (creal (w), cimag (w));
	}
*/

      z = mobius_eval (fuchsianG[i], cos (M_PI / 4) + I * sin (M_PI / 4));
      x1 = creal (z);
      y1 = cimag (z);
      z = mobius_eval (fuchsianG[i], cos (M_PI / 5) + I * sin (M_PI / 5));
      x2 = creal (z);
      y2 = cimag (z);
      z = mobius_eval (fuchsianG[i], cos (M_PI / 7) + I * sin (M_PI / 7));
      x3 = creal (z);
      y3 = cimag (z);


      get_circle_center (x1, y1, x2, y2, x3, y3, (double *) &center);

      glVertex2f (center[0], center[1]);

  //   plot_circle(center, sqrt( pow(x1-center[0],2) + pow(y1-center[1],2)));

    }
glEnd();
  //glEnd ();
  glutSwapBuffers ();
  ga = G.A.a;
  gb = G.A.b;
  gc = G.A.c;
  gd = G.A.d;

  return;
}




void
winInit ()
{
  gluOrtho2D (-1.5, 1.5, -1.5, 1.5);
  // glEnable (GL_LINE_SMOOTH);   
  glEnable (GL_POINT_SMOOTH);
  glHint (GL_POINT_SMOOTH_HINT, GL_NICEST);
  // glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
  // glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
  // glEnable(GL_BLEND);
  // glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}

void
display (void)
{
  double x, y;
  complex z, w;
  // glClearColor(1,1,1,1);
  //glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);


  for (x = -1.87; x < 3.0; x += 0.01)
    {
      firsttime = 1;
      for (y = -0.4; y < 0.9; y += 0.01)
	{
	  z = x + I * y;
	  w = x - I * y;
	  //printf ("Tr z = %f+I%f\t Tr w = %f + I%f\n", creal (z), cimag (z),  creal (w), cimag (w));
	  kleinian_action_s1 (z, w);
	}

    }
}


int
main (int argc, char **argv)
{
  srand (getpid ());

  


  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize (1024, 1024);
  glutCreateWindow ("Kleinianos");
  glClearColor (1, 1, 1, 1);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  winInit ();
  glutDisplayFunc (display);
  glutMainLoop ();
  return 0;
}

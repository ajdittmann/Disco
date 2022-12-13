//
// This code was created by Jeff Molofee '99 (ported to Linux/GLUT by Richard Campbell '99)
// If you've found this code useful, please let me know.
//
// Visit me at www.demonews.com/hosted/nehe 
// (email Richard Campbell at ulmont@bellsouth.net)
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <fstream>

#include <hdf5.h>

#include "gifer.h"

#ifdef OSX
#include <GLUT/glut.h>    // Header File For The GLUT Library 
#include <OpenGL/gl.h>	// Header File For The OpenGL32 Library
#include <OpenGL/glu.h>	// Header File For The GLu32 Library
#else
#include <GL/glut.h>    // Header File For The GLUT Library 
#include <GL/gl.h>	// Header File For The OpenGL32 Library
#include <GL/glu.h>	// Header File For The GLu32 Library
#endif
#include <unistd.h>     // needed to sleep

/* ASCII code for the escape key. */
#define ESCAPE 27

#define VAL_FLOOR 0.0 // 1.0 - 1.0e-10 //1 //0.95//0 //-8e-3//-3e-2 //(-HUGE_VAL) //.96
#define VAL_CEIL  2.0 //+ 1.0e-10 //-1 //4.5e-3 //1.05//5.25e-21 //5.25e-9 //8e-3//3e-2 //5.24e-5 //(HUGE_VAL)  //1.04
#define FIXMAXMIN 1
#define COLORMAX 13
#define CAM_BACKUP  1.5
#define ZRORDER 1 // 1: checkpoints organized with faster index r (default,new)
                  // 0: checkpoints organized with slower index r (old)
#define DIAGMODE 1 // 1: checkpoints store RZ diagnostics (default,new)
                   // 0: checkpoints store R diagnostics (old)

static int WindowWidth  = 800;
static int WindowHeight = 600;
 
int CommandMode;
int FullScreenMode=0;

int dim3d = 0;
int t_off = 0;
int p_off = 0;
int flipCM = 0;
int cmap = 6;
int draw_1d = 0;
int draw_bar = 0;
int draw_t   = 0;
int draw_spiral = 0;
int draw_jet = 0;
int draw_planet = 0;
int draw_scale = 0;
int valq=0;
int draw_border = 0;
int logscale = 0;
int floors=0;
int help_screen=0;
int print_vals=0;
int fix_zero=0;
int geometry = 0;

int first_frame = 1;

double val_floor = VAL_FLOOR;
double val_ceil = VAL_CEIL;

double rotate_angle = M_PI/2.;

double offx, offy, rescale;
double maxval = -HUGE_VAL;
double minval = HUGE_VAL;

double t;
int Nr,Nz,Nq,Npl,NpDat,N1d,Nc,KK, layer;
int * Np = NULL;
double * r_jph = NULL;
double * z_kph = NULL;
double ** p_iph = NULL;
double phimax;
double *** theZones = NULL;
double ** rzZones = NULL;
double ** thePlanets = NULL;
double ** theRadialData = NULL;
int *Np_All = NULL;
int *Tindex_All = NULL;
int *Id_phi0 = NULL;

double max_1d = 1.0;

char filename[1024];
char **filenameList = NULL;
int nfiles = 0;
int currentFile = 0;

void get_rgb( double , float * , float * , float * , int, int );
void loadSliceZ(char *filename, int k);
void loadSlicePhi(char *filename);
void loadDiagnostics(char *filename, int k);
void loadTime(char *filename);
void loadFile(int fileIndex, int zslice);
void DrawGLScene();

double getval( double * thisZone , int q ){
   if( q!=-1 ) return( thisZone[q] );
   double rho = thisZone[0];
//   double X   = thisZone[5];
   double P   = thisZone[1];
   double vr = thisZone[2];
   double om = thisZone[3];

/*
   double gam = 1.333333;
   double cs = sqrt(gam*P/rho);
    return fabs(vr)/cs;
    */

   double Br = thisZone[5];
   double Bp = thisZone[6];
   double Bz = thisZone[7];
   return  sqrt(P/rho); //return( sqrt(5*P/(3*rho)) ); //fabs(P/pow(rho,5./3.)-1.) );// fabs(thisZone[1]/pow(thisZone[0],5./3.)-1.) );
//   return( rho*P );
}

void getMaxMin(void){
   
   maxval = -HUGE_VAL;
   minval = HUGE_VAL;
   int q = valq;
   double val;
   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      for( i=0 ; i<Np[j] ; ++i ){
         val = getval( theZones[j][i], q );
         if( logscale ) val = log(val)/log(10.);
         if( maxval < val ) maxval = val;
         if( minval > val ) minval = val;
      }
      if( dim3d ){
         for( k=0 ; k<Nz ; ++k ){
            int jk;
            if(ZRORDER)
                jk = k*Nr+j;
            else
                jk = j*Nz+k;
            val = getval( rzZones[jk], q );
            if( logscale ) val = log(val)/log(10.);
            if( maxval < val ) maxval = val;
            if( minval > val ) minval = val;
         }
      }
   }
   if( floors == 1){
      if( maxval > VAL_CEIL  ) maxval = VAL_CEIL;
      if( minval < VAL_FLOOR ) minval = VAL_FLOOR;
   }
   if( FIXMAXMIN && floors==1 ){
      maxval = VAL_CEIL;
      minval = VAL_FLOOR;
   }
   if(floors == 2)
    {
       maxval = val_ceil; 
       minval = val_floor;
    }
   if( fix_zero ){
      //maxval = .5*(maxval-minval);
      if(-minval > maxval) maxval = -minval;
      minval = -maxval;
   }
   //if( floors ) minval = maxval-5.0;
   //if( floors ) minval = log( getval( theZones[0][Np[0]-1] , q ) )/log(10.);
}

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readString(char *file, char *group, char *dset, char *buf, int len)
{
   hid_t strtype = H5Tcopy(H5T_C_S1);
   H5Tset_size(strtype, H5T_VARIABLE);

   hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );

   if(!H5Lexists(h5grp, dset, H5P_DEFAULT))
   {
       buf[0] = '\0';
       return;
   }
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hid_t space = H5Dget_space(h5dst);
   hsize_t dims[1];
   H5Sget_simple_extent_dims(space, dims, NULL);

   char **rdata = (char **)malloc(dims[0] * sizeof(char *));
   H5Dread(h5dst, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
   strncpy(buf, rdata[0], len-1);
   buf[len-1] = '\0';

   H5Dvlen_reclaim(strtype, space, H5P_DEFAULT, rdata);
   free(rdata);
   H5Dclose(h5dst);
   H5Sclose(space);
   H5Tclose(strtype);
   H5Gclose(h5grp);
   H5Fclose(h5fil);
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

/* The number of our GLUT window */
int window; 

// Here are the fonts: 
#ifdef OSX
void * glutFonts[7] = { 
#else
void * glutFonts[7] = { 
#endif
    GLUT_BITMAP_9_BY_15, 
    GLUT_BITMAP_8_BY_13, 
    GLUT_BITMAP_TIMES_ROMAN_10, 
    GLUT_BITMAP_TIMES_ROMAN_24, 
    GLUT_BITMAP_HELVETICA_10, 
    GLUT_BITMAP_HELVETICA_12, 
    GLUT_BITMAP_HELVETICA_18 
}; 

// This function prints some text wherever you want it. 
#ifdef OSX
void glutPrint(float x, float y, float z , void ** font, char* text, float r, float g, float b, float a) 
#else
void glutPrint(float x, float y, float z , void * font, char* text, float r, float g, float b, float a) 
#endif

{ 
    if(!text || !strlen(text)) return; 
    int blending = 0; 
    if(glIsEnabled(GL_BLEND)) blending = 1; 
    glEnable(GL_BLEND); 
    glColor4f(r,g,b,a); 
    glRasterPos3f(x,y,z); 
    while (*text) { 
        glutBitmapCharacter(font, *text); 
        text++; 
    } 
    if(!blending) glDisable(GL_BLEND); 
}  
 
// A general OpenGL initialization function.  Sets all of the initial parameters. 
void InitGL(int Width, int Height)	        // We call this right after our OpenGL window is created.
{
  glClearColor(0.8f, 0.8f, 0.8f, 0.0f);		// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			        // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);		        // Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void ReSizeGLScene(int Width, int Height)
{
  if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
    Height=1;

#ifdef OSX
  //glutReshapeWindow(Width, Height);
#endif

  glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);

  WindowHeight = Height;
  WindowWidth = Width;
}

/*
void TakeScreenshot(const char *useless){
   int dimx = WindowWidth;
   int dimy = WindowHeight;
   printf("Writing out.ppm... ");
   float *pixels = new float[3*dimx*dimy];
   glReadBuffer(GL_BACK);
   glPixelStorei(GL_PACK_ALIGNMENT,1);
   glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_FLOAT, pixels);

   char *filename = (char *)"out.ppm";
   std::ofstream fp(filename);
   fp << "P3\n" << dimx << " " << dimy << std::endl << "255\n";
   int i,j;
   for( i=dimy-1; i>=0; i--){    
      for( j=0; j<dimx; ++j ){    
         int pixelPos = (i*dimx+j)*3;
         int r = (int)(pixels[pixelPos+0]*255.0);
         int g = (int)(pixels[pixelPos+1]*255.0);
         int b = (int)(pixels[pixelPos+2]*255.0);
         fp << r << " " << g << " " << b << std::endl;
      }    
   }    
   fp.close();
   printf("done!\n");
}
*/

void makePalette(int pal[])
{
    pal[0] = 0;
    pal[1] = 0;
    pal[2] = 0;
    pal[3] = 220;
    pal[4] = 220;
    pal[5] = 220;

    int i;
    for(i=2; i<256; i++)
    {
        float val = (i-2)/253.0;
        float rrr, ggg, bbb;
        get_rgb(val, &rrr, &ggg, &bbb, cmap, flipCM);
        pal[3*i+0] = 255*rrr;
        pal[3*i+1] = 255*ggg;
        pal[3*i+2] = 255*bbb;
    }
}

void calcGIFpixels(int dimx, int dimy, int *gifstream)
{
    float *pixels = (float *)malloc(3*dimx*dimy * sizeof(float));
    float *pixels_bw = (float *)malloc(dimx*dimy * sizeof(float));

    glReadBuffer(GL_BACK);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_FLOAT, pixels);

    float bgcolor[4];
    glGetFloatv(GL_COLOR_CLEAR_VALUE, bgcolor);

    int old_cmap = cmap;
    cmap = 5;  // black-to-white, rrr = ggg = bbb = val
    DrawGLScene();
    DrawGLScene();

    glReadBuffer(GL_BACK);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels(0, 0, dimx, dimy, GL_RED, GL_FLOAT, pixels_bw);
   
    cmap = old_cmap;
    DrawGLScene();

    int i, j;

    for(j=0; j<dimy; j++)
        for(i=0; i<dimx; i++)
        {
            int pixIdx = dimx*j + i;
            int gifIdx = dimx*(dimy-j-1) + i;
            float rrr = pixels[3*pixIdx];
            float ggg = pixels[3*pixIdx+1];
            float bbb = pixels[3*pixIdx+2];

            if(rrr == 0.0f && ggg == 0.0f && bbb == 0.0f)
                gifstream[gifIdx] = 0;
            else if(fabs(rrr - bgcolor[0]) < 1.0e-7
                    && fabs(ggg - bgcolor[1]) < 1.0e-7
                    && fabs(bbb - bgcolor[2]) < 1.0e-7)
                gifstream[gifIdx] = 1;
            else
                gifstream[gifIdx] = (int)(253 * pixels_bw[pixIdx]) + 2;
        }

    free(pixels);
    free(pixels_bw);
}

void gifify()
{
    int palette[3*256];
    makePalette(palette);

    int dimx = glutGet(GLUT_WINDOW_WIDTH); //WindowWidth;
    int dimy = glutGet(GLUT_WINDOW_HEIGHT); //WindowHeight;

    char gifname[1024 + 5];
    strcpy(gifname, filename);

    char *dot_ptr = strrchr(gifname, '.');
    if(dot_ptr != NULL)
        *dot_ptr = '\0';
    strcat(gifname, ".gif" );

    int *gifstream = (int *)malloc(dimx*dimy * sizeof(int));;

    calcGIFpixels(dimx, dimy, gifstream);
   
    printf("Saving %s\n", gifname);
    makeGIF(dimx, dimy, 1, palette, 256, 1, -1, 1.0, gifstream, gifname);

    free(gifstream);

}

void gifify_all()
{
    int palette[3*256];
    makePalette(palette);

    int dimx = glutGet(GLUT_WINDOW_WIDTH); //WindowWidth;
    int dimy = glutGet(GLUT_WINDOW_HEIGHT); //WindowHeight;

    char gifname[] = "disco_run.gif";

    int *gifstream = (int *)malloc(nfiles * dimx*dimy * sizeof(int));
    int i;

    int old_file = currentFile;
    for(i=0; i<nfiles; i++)
    {
        currentFile = i;
        loadFile(currentFile, KK);
        DrawGLScene();
        calcGIFpixels(dimx, dimy, gifstream + dimx*dimy*i);
    }
    currentFile = old_file;
    loadFile(currentFile, KK);
   
    printf("Saving %s\n", gifname);
    makeGIF(dimx, dimy, nfiles, palette, 256, 1, -1, 10.0, gifstream, gifname);

    free(gifstream);

}

void draw1dBackground(double RotationAngleX, double RotationAngleY,
                        double RotationAngleZ, double camdist, double xoff,
                        double yoff, double zoff, int q)
{
      glColor3f(1.0,1.0,1.0);
      glBegin(GL_QUADS);
      glVertex3f(-1./rescale-xoff,-1./rescale-yoff, camdist + zoff );
      glVertex3f( 1./rescale-xoff,-1./rescale-yoff, camdist + zoff );
      glVertex3f( 1./rescale-xoff, 1./rescale-yoff, camdist + zoff );
      glVertex3f(-1./rescale-xoff, 1./rescale-yoff, camdist + zoff );
      glEnd();
      glLineWidth( 3.0f );
      glColor3f(0.0,0.0,0.0);
}
void draw1dRadialData(double RotationAngleX, double RotationAngleY, 
                        double RotationAngleZ, double camdist, double xoff, 
                        double yoff, double zoff, int q)
{
    glBegin( GL_LINE_STRIP );
    int j;
    //      for( j=0 ; j<Nr ; ++j ){
    for( j=0 ; j<Nr ; ++j )
    {
        //if(logscale) val = log(getval(theZones[i],q))/log(10.);
        double rm = r_jph[j-1];
        double rp = r_jph[j]; 
        double val = theRadialData[j][q]/max_1d;//*pow(.5*(rp+rm),1.5);//2.*(theRadialData[j][q]-minval)/(maxval-minval) - 1.;//getval(theZones[i],q);
        if( q==1 ) val *= 1e2;
        val = 2.*val - 1.;
        double c0 = 2.*((.5*(rp+rm)-r_jph[-1])/(r_jph[Nr-1]-r_jph[-1]) - .5)/rescale;
        double c1 = val/rescale;
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+.001 );
    }    
    glEnd();
    /*
    glColor3f(0.5,0.5,0.5);
    glBegin( GL_LINE_STRIP );
    for( j=0 ; j<Nr ; ++j )
    {
        //if(logscale) val = log(getval(theZones[i],q))/log(10.);
        double rm = r_jph[j-1];
        double rp = r_jph[j];
        double r = .5*(rp+rm); 
        double val = pow(r,-1.5)/max_1d;// *pow(.5*(rp+rm),1.5);//2.*(theRadialData[j][q]-minval)/(maxval-minval) - 1.;//getval(theZones[i],q);         if( q==1 ) val *= 1e2;
        val = 2.*val - 1.;
        double c0 = 2.*((.5*(rp+rm)-r_jph[-1])/(r_jph[Nr-1]-r_jph[-1]) - .5)/rescale;                 double c1 = val/rescale;
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+.001 );
    }       
    glEnd();
    */
}

void draw1dPlanet(double RotationAngleX, double RotationAngleY,
                double RotationAngleZ, double camdist, double xoff,
                double yoff, double zoff, int q)
{
    double eps = 0.01;
    int p;
    for( p=0 ; p<Npl ; ++p )
    {
        double r   = thePlanets[p][0];
        double x = 2.*(( r - r_jph[-1] )/( r_jph[Nr-1]-r_jph[-1] ) - .5)/rescale;
        glColor3f(0.0,0.0,0.0);
        glLineWidth(2.0);
        glBegin(GL_LINE_LOOP);
        glVertex3f( x-xoff+eps , 1./rescale-yoff+eps , camdist+.001 );
        glVertex3f( x-xoff-eps , 1./rescale-yoff+eps , camdist+.001 );
        glVertex3f( x-xoff-eps , 1./rescale-yoff-eps , camdist+.001 );
        glVertex3f( x-xoff+eps , 1./rescale-yoff-eps , camdist+.001 );
        glEnd();
    }
}

void drawZCell(double x1p, double x1m, double x2p, double x2m, double x3,
                double camdist, double xoff, double yoff, double zoff)
{
    if(geometry == 1)
    {
        double rp = x1p/rescale;
        double rm = x1m/rescale;
        double phip = x2p;
        double phim = x2m;
        double dp = phip-phim;
        double z0 = x3/rescale;

        double c0 = rm*cos(phim);
        double c1 = rm*sin(phim);
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

        int np = (int) (rp*(phip-phim) / (rp-rm)) + 1;
        int nm = (int) (rm*(phip-phim) / (rp-rm)) + 1;

        c0 = rp*cos(phim);
        c1 = rp*sin(phim);
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

        int p;
        for(p=1; p<np; p++)
        {
            double ph = phim + (p*(phip-phim))/np;
            c0 = rp*cos(ph)/cos(dp/np);
            c1 = rp*sin(ph)/cos(dp/np);
            glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );
        }

        c0 = rp*cos(phip);
        c1 = rp*sin(phip);
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

        c0 = rm*cos(phip);
        c1 = rm*sin(phip);
        glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

        for(p=nm-1; p>0; p--)
        {
            double ph = phim + (p*(phip-phim))/nm;
            c0 = rm*cos(ph);
            c1 = rm*sin(ph);
            glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );
        }
    }
    else if(geometry == 2)
    {
        double rp = x1p/rescale;
        double rm = x1m/rescale;
        double phip = x2p;
        double phim = x2m;
        double sinth = sin(x3);
        double costh = cos(x3);
        double dp = phip-phim;

        double c0 = rm*sinth*cos(phim);
        double c1 = rm*sinth*sin(phim);
        double c2 = rm*costh;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int np = (int) (rp*(phip-phim) / (rp-rm)) + 1;
        int nm = (int) (rm*(phip-phim) / (rp-rm)) + 1;

        c0 = rp*sinth*cos(phim);
        c1 = rp*sinth*sin(phim);
        c2 = rp*costh;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int p;
        for(p=1; p<np; p++)
        {
            double ph = phim + (p*(phip-phim))/np;
            c0 = rp*sinth*cos(ph)/cos(dp/np);
            c1 = rp*sinth*sin(ph)/cos(dp/np);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }

        c0 = rp*sinth*cos(phip);
        c1 = rp*sinth*sin(phip);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        c0 = rm*sinth*cos(phip);
        c1 = rm*sinth*sin(phip);
        c2 = rm*costh;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        for(p=nm-1; p>0; p--)
        {
            double ph = phim + (p*(phip-phim))/nm;
            c0 = rm*sinth*cos(ph);
            c1 = rm*sinth*sin(ph);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }
    }
    else if(geometry == 3)
    {
        double mup = x1p;
        double mum = x1m;
        double phip = x2p;
        double phim = x2m;
        double s_nu = sin(x3);
        double c_nu = cos(x3);
        double dp = phip-phim;

        double rm = sinh(mum) * s_nu / rescale;
        double zm = cosh(mum) * c_nu / rescale; 
        double rp = sinh(mup) * s_nu / rescale;
        double zp = cosh(mup) * c_nu / rescale; 

        double c0 = rm*cos(phim);
        double c1 = rm*sin(phim);
        double c2 = zm;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int np = (int) (rp*(phip-phim) / (rp-rm)) + 1;
        int nm = (int) (rm*(phip-phim) / (rp-rm)) + 1;

        c0 = rp*cos(phim);
        c1 = rp*sin(phim);
        c2 = zp;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int p;
        for(p=1; p<np; p++)
        {
            double ph = phim + (p*(phip-phim))/np;
            c0 = rp*cos(ph)/cos(dp/np);
            c1 = rp*sin(ph)/cos(dp/np);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }

        c0 = rp*cos(phip);
        c1 = rp*sin(phip);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        c0 = rm*cos(phip);
        c1 = rm*sin(phip);
        c2 = zm;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        for(p=nm-1; p>0; p--)
        {
            double ph = phim + (p*(phip-phim))/nm;
            c0 = rm*cos(ph);
            c1 = rm*sin(ph);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }
    }
    else
    {
        double xm = x1m/rescale;
        double xp = x1p/rescale;
        double ym = x2m/rescale;
        double yp = x2p/rescale;
        double z = x3/rescale;
        glVertex3f( xm-xoff, ym-yoff, camdist-zoff+z );
        glVertex3f( xp-xoff, ym-yoff, camdist-zoff+z );
        glVertex3f( xp-xoff, yp-yoff, camdist-zoff+z );
        glVertex3f( xm-xoff, yp-yoff, camdist-zoff+z );
    }
}

void drawPhiCell(double x1p, double x1m, double x3p, double x3m, double x2,
                    double camdist, double xoff, double yoff, double zoff)
{
    if(geometry == 1)
    {
        double cp = cos(x2);
        double sp = sin(x2);
        double rm = x1m/rescale;
        double rp = x1p/rescale;
        double zm = x3m/rescale;
        double zp = x3p/rescale;

        glVertex3f( rm*cp-xoff, rm*sp-yoff, camdist-zoff+zm );
        glVertex3f( rp*cp-xoff, rp*sp-yoff, camdist-zoff+zm );
        glVertex3f( rp*cp-xoff, rp*sp-yoff, camdist-zoff+zp );
        glVertex3f( rm*cp-xoff, rm*sp-yoff, camdist-zoff+zp );
    }
    else if(geometry == 2)
    {
        double rp = x1p/rescale;
        double rm = x1m/rescale;
        double thp = x3p;
        double thm = x3m;
        double sinp = sin(x2);
        double cosp = cos(x2);
        double dth = thp-thm;

        double c0 = rm*sin(thm)*cosp;
        double c1 = rm*sin(thm)*sinp;
        double c2 = rm*cos(thm);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int np = (int) (rp*(thp-thm) / (rp-rm)) + 1;
        int nm = (int) (rm*(thp-thm) / (rp-rm)) + 1;

        c0 = rp*sin(thm)*cosp;
        c1 = rp*sin(thm)*sinp;
        c2 = rp*cos(thm);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        int p;
        for(p=1; p<np; p++)
        {
            double th = thm + (p*(thp-thm))/np;
            c0 = rp*sin(th)*cosp/cos(dth/np);
            c1 = rp*sin(th)*sinp/cos(dth/np);
            c2 = rp*cos(th);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }

        c0 = rp*sin(thp)*cosp;
        c1 = rp*sin(thp)*sinp;
        c2 = rp*cos(thp);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        c0 = rm*sin(thp)*cosp;
        c1 = rm*sin(thp)*sinp;
        c2 = rm*cos(thp);
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);

        for(p=nm-1; p>0; p--)
        {
            double th = thm + (p*(thp-thm))/nm;
            c0 = rm*sin(th)*cosp;
            c1 = rm*sin(th)*sinp;
            c2 = rp*cos(th);
            glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff);
        }
    }
    else if(geometry == 3)
    {
        double mup = x1p;
        double mum = x1m;
        double nup = x3p;
        double num = x3m;
        double sinp = sin(x2);
        double cosp = cos(x2);

        double sh_mum = sinh(mum);
        double ch_mum = cosh(mum);
        double sh_mup = sinh(mup);
        double ch_mup = cosh(mup);
        double s_num = sin(num);
        double c_num = cos(num);
        double s_nup = sin(nup);
        double c_nup = cos(nup);
        
        double rmm = sh_mum * s_num / rescale;
        double rmp = sh_mum * s_nup / rescale;
        double rpm = sh_mup * s_num / rescale;
        double rpp = sh_mup * s_nup / rescale;

        double zmm = ch_mum * c_num / rescale;
        double zmp = ch_mum * c_nup / rescale;
        double zpm = ch_mup * c_num / rescale;
        double zpp = ch_mup * c_nup / rescale;

        double c0 = rmm*cosp;
        double c1 = rmm*sinp;
        double c2 = zmm;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        c0 = rpm*cosp;
        c1 = rpm*sinp;
        c2 = zpm;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        c0 = rpp*cosp;
        c1 = rpp*sinp;
        c2 = zpp;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );

        c0 = rmp*cosp;
        c1 = rmp*sinp;
        c2 = zmp;
        glVertex3f( c0-xoff, c1-yoff, c2+camdist-zoff );
    }
    else
    {
        double xm = x1m/rescale;
        double xp = x1p/rescale;
        double zm = x3m/rescale;
        double zp = x3p/rescale;
        double y = x2/rescale;
        glVertex3f( xm-xoff, y-yoff, camdist-zoff+zm );
        glVertex3f( xp-xoff, y-yoff, camdist-zoff+zm );
        glVertex3f( xp-xoff, y-yoff, camdist-zoff+zp );
        glVertex3f( xm-xoff, y-yoff, camdist-zoff+zp );
    }
}

void drawZSlice(double camdist, double xoff, double yoff, double zoff, int q,
                 int draw_border_now)
{
    if( draw_border_now )
        zoff -= .0001;
    int i,j;

    for( j=0 ; j<Nr ; ++j )
    {
        double rm = r_jph[j-1];
        double rp = r_jph[j];
        double rc = .5*(rp+rm);
        double dr = rp-rm;

        if( !draw_border )
        {
            rm = rc-.55*dr;
            rp = rc+.55*dr;
        }
        for( i=0 ; i<Np[j] ; ++i )
        {
            double phip = p_iph[j][i];
            double phim;
            if( i==0 ) 
                phim = p_iph[j][Np[j]-1]; 
            else 
                phim = p_iph[j][i-1];

            if( t_off && !p_off )
            {
                phip -= t;
                phim -= t;
            }
            if( p_off )
            { 
                phip -= thePlanets[1][1]; 
                phim -= thePlanets[1][1];
            }

            double dp   = phip-phim;
            while( dp < 0.0 ) dp += phimax;
            while( dp > phimax ) dp -= phimax;
            double phi = phip - 0.5*dp;
            phim = phip - dp;

            if( !draw_border )
            {
                phip = phi+.55*dp;
                phim = phi-.55*dp;
            }

            double val = (getval(theZones[j][i],q)-minval)/(maxval-minval);
            if(logscale) 
                val = (log(getval(theZones[j][i],q))/log(10.)-minval)
                        /(maxval-minval);
            if( val > 1.0 ) val = 1.0;
            if( val < 0.0 ) val = 0.0;
            //double u = getval( theZones[j][i] , 2 );
            //if( uMax < u ) uMax = u;

            float rrr,ggg,bbb;
            get_rgb( val , &rrr , &ggg , &bbb , cmap, flipCM );

            if( (!dim3d || (sin(phi)>0 || cos(phi+.25)<0.0)) && dim3d !=2 )
            {
                //if( dp < 1.5 && (!dim3d || (sin(phi)>0 && cos(phi)>0.0)) && dim3d !=2 ){
                if( !draw_border_now )
                { 
                    glColor3f( rrr , ggg , bbb );
                    glBegin(GL_POLYGON);
                }
                else
                {
                    glLineWidth(3.0f);
                    glColor3f(0.0,0.0,0.0);
                    glBegin(GL_LINE_LOOP);
                }

                double z0 = 0.0;
                //if( dim3d ) 
                    z0 = 0.5*(z_kph[KK-1] + z_kph[KK]);

                drawZCell(rp, rm, phip, phim, z0, camdist, xoff, yoff, zoff);

                glEnd();
            }
        }
    }
    
    if( draw_border_now )
        zoff += .0001;
}

void drawPhiSlice(double camdist, double xoff, double yoff, double zoff, int q,
                    int draw_border_now)
{
    if( draw_border_now )
        zoff -= .0001;

    int j, k;
    for( j=0 ; j<Nr ; ++j )
    {
        for( k=0 ; k<Nz ; ++k )
        {
            int jk;
            if(ZRORDER)
                jk = k*Nr+j;
            else
                jk = j*Nz+k;
            double rp = r_jph[j];
            double rm = r_jph[j-1];
            double zp = z_kph[k];
            double zm = z_kph[k-1];

            double phi = 0.0;//rzZones[jk][Nq];

            double val = (getval(rzZones[jk],q)-minval)/(maxval-minval);
            if(logscale)
                val = (log(getval(rzZones[jk],q))/log(10.)-minval)
                            /(maxval-minval);
            if( val > 1.0 ) 
                val = 1.0;
            if( val < 0.0 ) 
                val = 0.0;
            float rrr,ggg,bbb;
            get_rgb( val , &rrr , &ggg , &bbb , cmap, flipCM );
            if( !draw_border_now )
            { 
                glColor3f( rrr , ggg , bbb );
                glBegin(GL_POLYGON);
            }
            else
            {
                glLineWidth(2.0f);
                glColor3f(0.0,0.0,0.0);
                glBegin(GL_LINE_LOOP);
            } 
                  
            drawPhiCell(rp, rm, zp, zm, phi, camdist, xoff, yoff, zoff);
            glEnd();
        }

        if( draw_border_now || dim3d == 1 )
        {
            glLineWidth(2.0f);
            glColor3f(0.0,0.0,0.0);
            glBegin(GL_LINE_LOOP);

            double rp = r_jph[Nr-1];
            //double rm =-r_jph[Nr-1];
            double rm = r_jph[-1];
            double zp = z_kph[Nz-1];
            double zm = z_kph[  -1];
            
            double phi = 0.0;
            drawPhiCell(rp, rm, zp, zm, phi, 0.0, xoff, yoff, zoff);
            glEnd();
        }
    }
    
    if( draw_border_now )
        zoff += .0001;
}

void drawData(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff, int q)
{
    if(dim3d < 2)
    {
        drawZSlice(camdist, xoff, yoff, zoff, q, 0);
        if(draw_border)
            drawZSlice(camdist, xoff, yoff, zoff, q, 1);
    }
    if(dim3d)
    {
        drawPhiSlice(camdist, xoff, yoff, zoff, q, 0);
        if(draw_border)
            drawPhiSlice(camdist, xoff, yoff, zoff, q, 1);
    }
}

// Draw a color bar.
void drawColorBar(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff, int q)
{
    glRotatef(-RotationAngleZ, 0, 0, 1);
    glRotatef(-RotationAngleY, 0, 1, 0);
    glRotatef(-RotationAngleX, 1, 0, 0);
    double xb = 0.6;
    double hb = 1.0;
    double wb = 0.02;
    int Nb = 1000;
    double dy = hb/(double)Nb;
    int k;
    for( k=0 ; k<Nb ; ++k )
    {
        double y = (double)k*dy - .5*hb;
        double val = (double)k/(double)Nb;
        float rrr,ggg,bbb;
        get_rgb( val , &rrr , &ggg , &bbb , cmap, flipCM );
        glLineWidth(0.0f);
        glColor3f( rrr , ggg , bbb );
        glBegin(GL_POLYGON);
        glVertex3f( xb    , y    , camdist + .001 ); 
        glVertex3f( xb+wb , y    , camdist + .001 ); 
        glVertex3f( xb+wb , y+dy , camdist + .001 ); 
        glVertex3f( xb    , y+dy , camdist + .001 ); 
        glEnd();
    }

    int Nv = 8;

    for( k=0 ; k<Nv ; ++k )
    {
        double y = (double)k*hb/(double)(Nv-1) - .5*hb;
        double val = ((double)k/(double)(Nv-1))*(maxval-minval) + minval;
        char valname[256];

        sprintf(valname,"%+.16e",val);

        glLineWidth(1.0f);
        glColor3f(0.0,0.0,0.0);
        glBegin(GL_LINE_LOOP);
        glVertex3f( xb    , y , camdist + .0011 );
        glVertex3f( xb+wb , y , camdist + .0011 );
        glEnd();
        glutPrint( xb+1.5*wb , y-.007 , camdist + .001 ,glutFonts[6] , valname , 0.0f, 0.0f , 0.0f , 0.5f );
    }

    glRotatef(RotationAngleX, 1, 0, 0);
    glRotatef(RotationAngleY, 0, 1, 0);
    glRotatef(RotationAngleZ, 0, 0, 1);
}

void drawPlanet(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff)
{
    glRotatef(-RotationAngleZ, 0, 0, 1);
    glRotatef(-RotationAngleY, 0, 1, 0);
    glRotatef(-RotationAngleX, 1, 0, 0);
    double eps = 0.01;
    int p;
    for( p=0 ; p<Npl ; ++p )
    {
        double r   = thePlanets[p][0]/rescale;
        double phi = thePlanets[p][1];
        if( t_off && !p_off ) phi -= t;
        if( p_off ) phi -= thePlanets[1][1];

        glColor3f(0.0,0.0,0.0);
        glLineWidth(2.0);

        glBegin(GL_LINE_LOOP);
        glVertex3f( r*cos(phi)-xoff+eps , r*sin(phi)-yoff+eps , camdist+.001 );
        glVertex3f( r*cos(phi)-xoff-eps , r*sin(phi)-yoff+eps , camdist+.001 );
        glVertex3f( r*cos(phi)-xoff-eps , r*sin(phi)-yoff-eps , camdist+.001 );
        glVertex3f( r*cos(phi)-xoff+eps , r*sin(phi)-yoff-eps , camdist+.001 );
        glEnd();
    }
    glRotatef(RotationAngleX, 1, 0, 0);
    glRotatef(RotationAngleY, 0, 1, 0);
    glRotatef(RotationAngleZ, 0, 0, 1);
}

void drawSpiral(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff)
{

    //double rp   = 5.0;//thePlanets[1][0];
    double Mach = 3.65;
    double p0 = 1.2;
    //double dr = .07;
    int Nr = 200;
    //double Rmin = 0.5;
    //double Rmax = 1.1;
    int k;
    
    glLineWidth(3.0f);
    glColor3f(0.0,0.0,0.0);
    //glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINE_LOOP);
    
    for( k=0 ; k<Nr ; ++k )
    {
        //double phi0 = ((double)k+.5)/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
        double x = ((double)k+.5)/(double)Nr;
        double r = .001*pow(.5/.001,x);
        //double r = .5*pow(1.5/.5,x);
        double phi0 = p0-log(r/0.001)*Mach;//(3.-2.*sqrt(1./r)-r)*5.;
        //double phi0 = (3.-2.*sqrt(1./r)-r)*20.;
        //if( r<1. ) phi0 = -phi0;

        //         double x0 = rp + dr*cos(phi0);
        //         double y0 = dr*sin(phi0);

        double phi = phi0;
        //         double r = rp;       

        //         double phi = atan2(y0,x0);
        //         double phi = ((double)k+.5)/(double)Nr*2.*M_PI*9.32 + 4.*M_PI;
        //         double r   = sqrt(x0*x0+y0*y0);
        //         double r = pow(phi/10.,-2.);
        //         double r   = 2./(1.+sin(phi));//1.0;//((double)k+0.5)/(double)Nr*(Rmax-Rmin) + Rmin;
        //         if( r<1. ) phi = -phi;
        r /= rescale;

        glVertex3f( r*cos(phi)-xoff , r*sin(phi)-yoff , camdist + .0011 );
        if( k%2==1 )
        {
            glEnd();
            glBegin(GL_LINE_LOOP);
        }
    }
    glEnd();

    /*
    e += 0.01;
    glBegin(GL_LINE_LOOP);
    for( k=0 ; k<Nr ; ++k ){
    double phi0 = ((double)k-.5)/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
    double x0 = 1.+e*cos(phi0);
    double y0 = e*sin(phi0);
    double phi = atan2(y0,x0);
    double r   = sqrt(x0*x0+y0*y0);
    //         double phi = (double)k/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
    //         double r   = 2./(1.+sin(phi));//1.0;//((double)k+0.5)/(double)Nr*(Rmax-Rmin) + Rmin;
    //         if( r<1. ) phi = -phi;
    r /= rescale;
    glVertex3f( r*cos(phi)-xoff , r*sin(phi)-yoff , camdist + .0011 );
    if( k%2==1 ){
    glEnd();
    glBegin(GL_LINE_LOOP);
    }
    }
    glEnd();
    */
}

void drawText(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff)
{
    char tprint[256];
    sprintf(tprint,"t = %.2e",t);
    glutPrint( -.6 , .5 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
    //    sprintf(tprint,"uMax = %.1f",uMax);
    //    glutPrint( -.8 , .4 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
}

void drawHelp(double RotationAngleX, double RotationAngleY,
                    double RotationAngleZ, double camdist, double xoff,
                    double yoff, double zoff)
{
    int NLines = 10;
    double dy = -.04;
    char help[NLines][256];
    sprintf(help[0],"Help Display:");
    sprintf(help[1],"b - Toggle Colorbar");
    sprintf(help[2],"B - Flip Colormap");
    sprintf(help[3],"[ / ] - Change Colormap");
    sprintf(help[4],"f - Toggle Max/Min Floors");
    sprintf(help[5],"p - Toggle Planet Data");
    sprintf(help[6],"1-9 - Choose Primitive Variable to Display");
    sprintf(help[7],"arrow keys - Move Camera");
    sprintf(help[8],"z/x - Zoom in/out");
    sprintf(help[9],"h - Toggle Help Screen");
    int i;
    for( i=0 ; i<NLines ; ++i )
    {
        glutPrint( -.8 , i*dy , camdist + .001 , glutFonts[6] , help[i] , 
                    0.0f, 0.0f , 0.0f , 0.5f );
    }
} 

/* 
 * The main drawing function. 
 */
void DrawGLScene(){

   glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
   glLoadIdentity();

   if(first_frame)
   {
       ReSizeGLScene(WindowWidth, WindowHeight);
       first_frame = 0;
   }


   double RotationAngleX = 0.0;
   double RotationAngleY = 0.0;
   double RotationAngleZ = 0.0;

   if( dim3d==1 ) RotationAngleX = -60.0;
   if( dim3d==2 ) RotationAngleX = -90.0;
 
   glTranslatef( 0.0 , 0.0 , -CAM_BACKUP );

   glRotatef(RotationAngleX, 1, 0, 0);
   glRotatef(RotationAngleY, 0, 1, 0);
   glRotatef(RotationAngleZ, 0, 0, 1);

   int q = valq;

   getMaxMin();
   if( print_vals==0 ){
      printf("Max = %e Min = %e\n",maxval,minval);
      print_vals=1;
   }

   double camdist = 0.0;//-CAM_BACKUP;
   double xoff    = offx;
   double yoff    = offy;
   double zoff    = 0.0;

   if( draw_1d )
   {
      draw1dBackground(RotationAngleX, RotationAngleY, RotationAngleZ,
                        camdist, xoff, yoff, zoff, q);
      if(theRadialData != NULL)
         draw1dRadialData(RotationAngleX, RotationAngleY, RotationAngleZ,
                            camdist, xoff, yoff, zoff, q);
      if( draw_planet )
         draw1dPlanet(RotationAngleX, RotationAngleY, RotationAngleZ,
                            camdist, xoff, yoff, zoff, q);
   }
   else
   {
       drawData(RotationAngleX, RotationAngleY, RotationAngleZ,
                        camdist, xoff, yoff, zoff, q);

       if( draw_bar )
          drawColorBar(RotationAngleX, RotationAngleY, RotationAngleZ, camdist,
                    xoff, yoff, zoff, q);

       if( draw_planet )
          drawPlanet(RotationAngleX, RotationAngleY, RotationAngleZ, camdist,
                        xoff, yoff, zoff);

       if( draw_spiral )
          drawSpiral(RotationAngleX, RotationAngleY, RotationAngleZ, camdist,
                        xoff, yoff, zoff);
   }

   if( draw_t )
       drawText(RotationAngleX, RotationAngleY, RotationAngleZ, camdist,
                        xoff, yoff, zoff);

   if( help_screen )
       drawHelp(RotationAngleX, RotationAngleY, RotationAngleZ, camdist,
                        xoff, yoff, zoff);

   if( CommandMode )
   {
      //TakeScreenshot("out.ppm");
      gifify();
      exit(0);
   }

   glutSwapBuffers();

}

/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
   usleep(100);
   if (key == ESCAPE){
      glutDestroyWindow(window); 
      exit(0);                   
   }
   if (key == '[' ){
       --cmap;
       cmap = (cmap+COLORMAX) % COLORMAX;
   }  
   if (key == ']' ){
       ++cmap;
       cmap = cmap % COLORMAX;
   }
   if( key >= (int)'0' && key < (int)'1'+Nq ) valq = (int)key-(int)'1';
   if( key == 'b' ) draw_bar = !draw_bar;
   if( key == 'B' ) flipCM = !flipCM;
   if( key == 'd' ) {++dim3d; if(dim3d==3) dim3d=0;}
   if( key == 'f' ){
       floors = (floors+1)%3;
       print_vals = 0;
   }
   if( key == 'g' ) draw_border = !draw_border;
   if( key == 'h' ) help_screen = !help_screen;
   if( key == 'l' ) logscale = !logscale;
   if( key == 's' ) draw_scale = !draw_scale;
   if( key == 't' ) draw_t = !draw_t;
   if( key == 'u' ) draw_1d = !draw_1d;
   if( key == 'x' ){
      rescale *= 1.3;
      offx /= 1.3;
      offy /= 1.3;
   }
   if( key == 'z' ){
      rescale /= 1.3;
      offx *= 1.3;
      offy *= 1.3;
   }
   if( key == 'j' ){
      max_1d *= 1.3;
   }
   if( key == 'J' ){
      max_1d /= 1.3;
   }
   if (key == 'F') {
      if (FullScreenMode) {
         glutReshapeWindow(WindowWidth, WindowHeight);
         glutPositionWindow(0,0);
         FullScreenMode = 0; 
      }else{
         glutFullScreen();
         FullScreenMode = 1; 
      }
   }
   if( key == 'p' ) draw_planet = !draw_planet;
   if( key == 's' ) draw_spiral = !draw_spiral;
   if( key == 'Z' ) fix_zero = !fix_zero;
   if( key == 'q' )
   {
       layer++;
       if(layer > 8)
          layer = 0;

       if(layer == 0)
           KK = 0;
       else if(layer == 1)
           KK = Nz/8;
       else if(layer == 2)
           KK = Nz/4;
       else if(layer == 3)
           KK = (3*Nz)/8;
       else if(layer == 4)
           KK = Nz/2;
       else if(layer == 5)
           KK = (5*Nz)/8;
       else if(layer == 6)
           KK = (3*Nz)/4;
       else if(layer == 7)
           KK = (7*Nz)/8;
       else
           KK = Nz - 1;
       loadSliceZ(filename, KK);
       loadDiagnostics(filename, KK);
   }
   if( key == ',')
   {
       currentFile--;
       while(currentFile < 0)
           currentFile += nfiles;
       loadFile(currentFile, KK);
   }
   if( key == '.')
   {
       currentFile++;
       while(currentFile >= nfiles)
           currentFile -= nfiles;
       loadFile(currentFile, KK);
   }
   if( key == '<')
   {
       currentFile -= 5;
       while(currentFile < 0)
           currentFile += nfiles;
       loadFile(currentFile, KK);
   }
   if( key == '>')
   {
       currentFile += 5;
       while(currentFile >= nfiles)
           currentFile -= nfiles;
       loadFile(currentFile, KK);
   }
   if(key == '/')
   {
       currentFile = 0;
       loadFile(currentFile, KK);
   }
   if(key == 'm')
   {
       val_ceil = maxval;
       val_floor = minval;
   }
   if(key == 'o')
   {
       gifify();
   }
   if(key == 'O')
   {
       gifify_all();
   }

   glutPostRedisplay();
}

void specialKeyPressed(int key, int x, int y){
   usleep(100);
   if (key == GLUT_KEY_LEFT ) offx -= .1;
   if (key == GLUT_KEY_RIGHT) offx += .1;

   if (key == GLUT_KEY_UP   ) offy += .1;
   if (key == GLUT_KEY_DOWN ) offy -= .1;
}

void freeData()
{
    int p, i, j;

   for( p=0 ; p<Npl ; ++p ){
      free(thePlanets[p]);
   }
   free(thePlanets);

   for( j=0 ; j<Nr ; ++j ){
      for( i=0 ; i<Np[j] ; ++i ){
         free( theZones[j][i] );
      }
      free( theZones[j] );
   }
   for( j=0 ; j<Nr*Nz ; ++j ){
      free( rzZones[j] );
   }
   free( rzZones );

   free( theZones );
   for( j=0 ; j<Nr ; ++j ){
      free( p_iph[j] );
      free( theRadialData[j] );
   }
   free( p_iph );
   if(theRadialData != NULL)
   {
       for( j=0 ; j<Nr ; ++j ){
          free( theRadialData[j] );
       }
       free( theRadialData );
   }
   free( Np );
   r_jph--;
   free( r_jph );
   z_kph--;
   free( z_kph );

   free(Np_All);
   free(Tindex_All);
   free(Id_phi0);
}

void loadTime(char *filename)
{
   char group1[256];
   strcpy( group1 , "Grid" );
   readSimple( filename , group1 , (char *)"T" , &t , H5T_NATIVE_DOUBLE );
}

void loadGrid(char *filename)
{
   char group1[256];
   char group2[256];
   char group3[256];
   char group4[256];
   strcpy( group1 , "Grid" );
   strcpy( group2 , "Data" );
   strcpy( group3 , "Pars" );
   strcpy( group4 , "Opts" );

   hsize_t dims[3];
   
   getH5dims( filename , group1 , (char *)"r_jph" , dims );
   Nr = dims[0]-1;
   getH5dims( filename , group1 , (char *)"z_kph" , dims );
   Nz = dims[0]-1;

   Np = (int *) malloc( Nr*sizeof(int) );
   r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
   z_kph = (double *) malloc( (Nz+1)*sizeof(double) );
   
   printf("t = %.2f, Nr = %d Nz = %d\n",t,Nr,Nz);

   readSimple( filename , group1 , (char *)"r_jph" , r_jph , H5T_NATIVE_DOUBLE );
   readSimple( filename , group1 , (char *)"z_kph" , z_kph , H5T_NATIVE_DOUBLE );
   r_jph++;
   z_kph++;

   readSimple( filename , group3 , (char *)"Phi_Max", &phimax, 
                    H5T_NATIVE_DOUBLE);

   char buf[256];

   readString(filename, group4, (char *)"GEOMETRY", buf, 256);
   if(strlen(buf) == 0)
       strcpy(buf, "cylindrical");

   printf("Geometry is: %s\n", buf);
   printf("%lu %d %d\n", strlen(buf), strcmp(buf,"cartesian"), 
                            strcmp(buf,"cylindrical"));

   if(strcmp(buf, "cylindrical") == 0)
      geometry = 1;
   else if(strcmp(buf, "spherical") == 0)
      geometry = 2;
   else if(strcmp(buf, "prolateSpheroidal") == 0)
      geometry = 3;
   else
      geometry = 0;

   int start[2]    = {0,0};
   int loc_size[2] = {Nz,Nr};
   int glo_size[2] = {Nz,Nr};
   if(!ZRORDER)
   {
       loc_size[0] = Nr; loc_size[1] = Nz;
       glo_size[0] = Nr; glo_size[1] = Nz;
   }

   Np_All = (int *) malloc(Nr*Nz * sizeof(int));
   Tindex_All = (int *) malloc(Nr*Nz * sizeof(int));
   Id_phi0 = (int *) malloc(Nr*Nz * sizeof(int));

   readPatch( filename , group1 , (char *)"Np"      , Np_All     , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( filename , group1 , (char *)"Index"   , Tindex_All , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( filename , group1 , (char *)"Id_phi0" , Id_phi0    , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);


   getH5dims( filename , group2 , (char *)"Cells" , dims );
   Nc = dims[0];
   Nq = dims[1]-1;
   printf("Nc = %d Nr = %d Nq=%d\n",Nc,Nr,Nq);
}

void loadPlanets(char *filename)
{
   char group2[256];
   strcpy( group2 , "Data" );

   hsize_t dims[3];
   int p;

   if(thePlanets == NULL)
   {
       getH5dims( filename , group2 , (char *)"Planets" , dims );
       Npl = dims[0];
       NpDat = dims[1];
       printf("Npl=%d NpDat=%d\n",Npl, NpDat);

       thePlanets = (double **) malloc( Npl*sizeof(double *) );
       for( p=0 ; p<Npl ; ++p ){ 
          thePlanets[p] = (double *) malloc( 2*sizeof(double) );
       }
   }
   
   int start[2]    = {0,0};
   int loc_size[2] = {1,NpDat};
   int glo_size[2] = {Npl,NpDat};
   for( p=0 ; p<Npl ; ++p ){
      start[0] = p;
      double thisPlanet[NpDat];
      readPatch( filename , group2 , (char *)"Planets" , thisPlanet , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
      int dp;
      printf("Planet = ");
      for( dp=0 ; dp<NpDat ; ++dp ) printf("%e ",thisPlanet[dp]);
      printf("\n");
      thePlanets[p][0] = thisPlanet[3];
      thePlanets[p][1] = thisPlanet[4];
   }
}

void loadFile(int fileIndex, int zslice)
{
    printf("File %d of %d\n", fileIndex+1, nfiles);
    strcpy(filename, filenameList[fileIndex]);
    loadTime(filenameList[fileIndex]);
    loadPlanets(filenameList[fileIndex]);
    loadSliceZ(filenameList[fileIndex], zslice);
    loadSlicePhi(filenameList[fileIndex]);
    loadDiagnostics(filenameList[fileIndex], zslice);
}

void loadSliceZ(char *filename, int k)
{
   char group1[256];
   char group2[256];
   strcpy( group1 , "Grid" );
   strcpy( group2 , "Data" );

   int Tindex[Nr];

   printf("Loading Z slice...\n");

   printf("Zones:");

   //Free arrays if already allocated
   int i,j;
   if(theZones != NULL)
   {
       for( j=0 ; j<Nr ; ++j )
       {
           for( i=0 ; i<Np[j] ; ++i )
               free(theZones[j][i]);
           free(theZones[j]);
       }
       free(theZones);
   }
   if(p_iph != NULL)
   {
       for( j=0 ; j<Nr ; ++j )
           free(p_iph[j]);
       free(p_iph);
   }

   //Grab zone indices
   for( j=0 ; j<Nr ; ++j ){
      int jk = k*Nr+j;
      if(!ZRORDER)
          jk = j*Nz+k;
      Np[j]     = Np_All[jk];
      Tindex[j] = Tindex_All[jk];
   }

   //Allocate arrays
   theZones = (double ***) malloc( Nr*sizeof(double **) );
   int q;
   for( j=0 ; j<Nr ; ++j ){
      theZones[j] = (double **) malloc( Np[j]*sizeof( double * ) );
      for( i=0 ; i<Np[j] ; ++i ){
         theZones[j][i] = (double *) malloc( Nq*sizeof( double ) );
      }
   }

   p_iph = (double **) malloc( Nr*sizeof( double * ) );
   for( j=0 ; j<Nr ; ++j ){
      p_iph[j] = (double *) malloc( Np[j]*sizeof( double ) );
   }

   printf("  Allocated,");

   //Load data

   int start[2]    = {0,0};
   int loc_size[2] = {0,Nq+1};
   int glo_size[2] = {Nc,Nq+1};

   for( j=0 ; j<Nr ; ++j ){
      loc_size[0] = Np[j];
      start[0] = Tindex[j];

      double TrackData[Np[j]*(Nq+1)];
      readPatch( filename , group2 , (char *)"Cells" , TrackData , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size);
      for( i=0 ; i<Np[j] ; ++i ){
         p_iph[j][i] = TrackData[i*(Nq+1) + Nq];
         for( q=0 ; q<Nq ; ++q ){
            theZones[j][i][q] = TrackData[i*(Nq+1) + q];
         }
      }
   }

   printf("  Initialized.\n");
}

void loadSlicePhi(char *filename)
{
    printf("Loading Phi slice...\n");
    
    int j,k;
    char group2[256];
    strcpy( group2 , "Data" );

    printf("rzZones:");

    if(rzZones != NULL)
    {
        for(j=0; j<Nr*Nz; j++)
            free(rzZones[j]);
        free(rzZones);
    }

    rzZones = (double **) malloc( Nr*Nz*sizeof(double *) );
    for( j=0 ; j<Nr*Nz ; ++j )
        rzZones[j] = (double *) malloc( (Nq+1)*sizeof( double ) );

    printf("  Allocated,");

    int glo_size[2] = {Nc, Nq+1};
    int loc_size[2] = {1, Nq+1};
    int start[2] = {0, 0};
   
   if(ZRORDER)
   {
       for( k=0 ; k<Nz ; ++k ){
          for( j=0 ; j<Nr ; ++j ){
             int jk = k*Nr+j;
             start[0] = Id_phi0[jk];
             readPatch( filename , group2 , (char *)"Cells" , rzZones[jk] , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
          }
       }
   }
   else
   {
       for( j=0 ; j<Nr ; ++j ){
          for( k=0 ; k<Nz ; ++k ){
             int jk = j*Nz+k;
             start[0] = Id_phi0[jk];
             readPatch( filename , group2 , (char *)"Cells" , rzZones[jk] , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
          }
       }
   }

   printf("  Initialized.\n");
}

void loadDiagnostics(char *filename, int k)
{
    int j;
    hsize_t dims[3];
    char group2[256];
    strcpy( group2 , "Data" );

    printf("Diagnostics:");

   if(DIAGMODE == 0)
   {
       getH5dims( filename , group2 , (char *)"Radial_Diagnostics" , dims );
       N1d = dims[1];
   }
   else if(DIAGMODE == 1)
   {
       //Diagnostics are stored in an array Nz x Nr x Ndiag
       getH5dims( filename , group2 , (char *)"Diagnostics" , dims );
       N1d = dims[2];
   }

   theRadialData = (double **) malloc( Nr*sizeof(double *) );
   for( j=0 ; j<Nr ; ++j )
      theRadialData[j] = (double *) malloc( N1d*sizeof(double) );
   
    printf("  Allocated,");

   if(DIAGMODE == 0)
   {
       int loc_size[2] = {1, N1d};
       int glo_size[2] = {Nr, N1d};
       int start[2] = {0, 0};
       for( j=0 ; j<Nr ; ++j ){
          start[0] = j;
          double thisQ[N1d];
          readPatch( filename , group2 , (char *)"Radial_Diagnostics" , thisQ , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
          int nq;
          for( nq = 0 ; nq < N1d ; ++nq ) theRadialData[j][nq] = thisQ[nq];
       }
   }
   else if(DIAGMODE == 1)
   {
       int loc_size3[3] = {1,1,N1d};
       int glo_size3[3] = {Nz,Nr,N1d};
       int start3[3] = {k,0,0};
       for( j=0 ; j<Nr ; ++j ){
          start3[1] = j;
          double thisQ[N1d];
          readPatch( filename , group2 , (char *)"Diagnostics" , thisQ , H5T_NATIVE_DOUBLE , 3 , start3 , loc_size3 , glo_size3 );
          int nq;
          for( nq = 0 ; nq < N1d ; ++nq ) theRadialData[j][nq] = thisQ[nq];
       }
   }

    printf("  Initialized.\n");
}

int main(int argc, char *argv[]) 
{
   if( argc < 2 ){
      printf("Please specify the input file.\n");
      exit(1);
   }

    filenameList = &(argv[1]);
    nfiles = argc - 1;
    currentFile = 0;

    strcpy( filename , filenameList[0] );
   CommandMode=0;
   //if( argc>2 ){ CommandMode=1; FullScreenMode=1; }

   loadTime(filename);
   loadGrid(filename);
   loadPlanets(filename);

   KK = Nz/2;
   layer = 4;
   loadSliceZ(filename, KK);
   loadSlicePhi(filename);
   loadDiagnostics(filename, KK);
   

   //double RR = r_iph[0][Nr[0]-1];
   //double r_max = .98*RR;
   //double r_min = 0.0;//r_iph[0][0];
   //double r_min = .82*RR;

   //printf("Rmin = %.2e Rmax = %.2e\n",r_min,r_max);
   //double thalf = .5*t_jph[Nt-1];
   rescale = 1.7*r_jph[Nr-1];  //(r_max-r_min);
   //offx = .5*(r_min+r_max)/rescale;
   //offy = 0.5;//.5*(r_min+r_max)*sin(thalf)/rescale;
   //offy = .5*(r_min+r_max)/rescale;
   //rescale *= .5;
   //rescale *= .4;
   rescale *= 1.2;
   //offy = 0.75/rescale;
   //offx = 0.6/rescale;//0.0;//1.0/rescale;
   offy = 0.0;
   offx = 0.0;//2.5/rescale;//0.0;//1.0/rescale;

//////////////////////////////
   glutInit(&argc, argv);  
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);  
   glutInitWindowSize(WindowWidth, WindowHeight);
   glutInitWindowPosition(0, 0);
   window = glutCreateWindow("DISCO Viewer");
   glutDisplayFunc(&DrawGLScene);  
   if(FullScreenMode) glutFullScreen();
   glutIdleFunc(&DrawGLScene);
   glutReshapeFunc(&ReSizeGLScene);
   glutKeyboardFunc(&keyPressed);
   glutSpecialFunc(&specialKeyPressed);
   InitGL(WindowWidth, WindowHeight);
   glutMainLoop();  
//////////////////////////////

   freeData();

   return 0;
}


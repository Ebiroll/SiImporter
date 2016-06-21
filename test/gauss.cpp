/*****************************************
 *
 *   Module: GAUSS.C SI 36011 40106
 *  
 *   Author:  Copy of GAUSHANN.PAS developed by Bo Johansson LFV Sturup
 *            converted to C by Kristian Tr|nnberg
 * 
 *   Description: Contains Interface functions for transformations
 *                GAUSS HANOVER but with same interface as 
 *                Slant Stereographical Projection model see slant.c
 * 
 *   Revision History:
 *   $Log$
 *   Revision 1.1  2015/07/06 08:03:15  olas
 *   Added aftsnim to sunnet
 *
 *   Revision 1.1  2013/09/13 08:50:30  olas
 *   Corrected projection method
 *
 *   Revision 1.1.1.1  2013/07/15 08:59:55  olas
 *   Import of ufs into ares
 *
 *   Revision 1.14  2013/06/14 15:14:40  alfa
 *   Adde fields for billing + log-file handling for it
 *
 *   Revision 1.13  2010/08/26 10:35:23  hegu
 *   Corrected "a_hat"
 *
 *   Revision 1.12  2010/08/25 14:33:53  hegu
 *   Corrected "a_hat"
 *
 *   Revision 1.11  2008/10/24 12:52:32  hegu
 *   corrected rounding error in SlGeoToLl2
 *
 *   Revision 1.10  2008/05/29 07:41:07  hegu
 *   added functions for converting lat/long with 100th sec precision
 *
 *   Revision 1.9  2008/05/28 17:18:50  hegu
 *   handling of cat19 and 20
 *
 *   Revision 1.8  2008/05/06 17:53:46  hegu
 *   increased order projection formulas series expansions to fourth order
 *
 *   Revision 1.7  2008/05/06 17:41:07  hegu
 *   increased order projection formulas series expansions to fourth order
 *
 *   Revision 1.6  2008/05/06 14:02:51  hegu
 *   increased order projection formulas series expansions to fourth order
 *
 *   Revision 1.5  2008/05/06 12:12:54  hegu
 *   increased order projection formulas series expansions to fourth order
 *
 *   Revision 1.4  2006/01/13 09:29:53  hegu
 *   Corrected rounding of lat/long seconds (59.5+ was incorrectly rounded to 60) in SlGeoToLl.
 *
 *   Revision 1.3  2005/04/18 14:18:56  hegu
 *   Fixed rounding error in SlGeoToLL
 *
 *   Revision 1.2  2005/04/07 08:56:15  hegu
 *   Prototypes should be in math.h
 *
 *   Revision 1.1  2003/11/21 16:11:19  hegu
 *   sl.h now called gauss.h and moved
 *   to common/util
 *
 *   Revision 1.6  2003/07/04 13:47:26  hegu
 *   removed hardcoded parameters. WGS84 earthmodel used (same as in MRT)
 *
 *   Revision 1.5  2003/07/03 16:27:23  hegu
 *   no message
 *
 *   Revision 1.2  2003/05/20 16:13:50  lafo
 *   Removed unused objects
 *
 *   Revision 1.1.1.1  2000/04/17 13:57:19  olas
 *   Build 14.3 from Riga
 *
 *   Revision 1.3  1999/02/17 10:56:42  alfa
 *   Build 9.0 Clean Up
 *
 *   Revision 1.2  1998/02/24 16:06:01  alfa
 *   Test on compiler variable __RIGA removed in GaussInit
 *   Caused wrong RefMeridian if not correct set in Makefile
 *
 *   Revision 1.1  1997/04/10 07:17:45  moor
 *   Initial revision
 *
 *   Revision 1.5  1996/06/10 18:17:46  krtr
 *   Ref meridian set at init of module
 *
 * Revision 1.4  95/12/19  11:44:56  krtr
 * Removed static decl on MeriaConv
 * 
 * Revision 1.3  94/04/22  13:54:53  krtr
 * Made function MeridianConv public
 * 
 * Revision 1.2  94/03/15  17:00:40  krtr
 * New SI no
 * 
 * Revision 1.1  94/01/16  13:55:42  krtr
 * Initial revision
 * 
 * 
 *   date     action
 *   ------   -------
 *   940112   First version
 ********************************************************************/
static char rcsid[]="$Header$";

 /* include constants */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gauss.h"

/*
* Modules vars
 */
  
static double  Rad2Deg = 180.0 / 3.141592653589793238; 
static double  Deg2Rad = 3.141592653589793238 / 180.0; 

static double  RefMeridian;
static double  Mean_Radius;
static double  origo_x ; 
static double  origo_y ;
static double  origo_lat;
static double  origo_long;
static double  beta1, beta2, beta3, beta4;
static double  delta1, delta2, delta3, delta4;
static double  ecc;
static double  ecc2, ecc4, ecc6,  ecc8;
static double  A0, B0, C0, D0; 
static double  A1, B1, C1, D1; 

static double earth_A; // earth equatorial radius
static double earth_B; // earth polar radius


/* forward decl */
void  ll2xy  ( double lat,double lng, double *x_ptr,double *y_ptr);
static void  xy2ll  (double *lat_ptr, double *long_ptr,
                      double x, double y);

static void  Origo_ll2xy  (   double  lat,double lng,double *x_ptr, double *y_ptr);
static void  Origo_xy2ll  (double *lat,double *lng, double  x,double y);

static void     SetOrigo (double lat,double lng);
static void     GetOrigo (double *lat, double* lng,double *x_ptr,double *y_ptr);
static void     SetEarthRadie ( double meter);

/****************************************************************************
* ll2xy

 Description: lat     latitude in radians
              long    longitude in radians
              x,y     calculated grid coordinates


    Function: Calculates grid coordinates x,y from the geographical
              coordinates lat, long in radians.


 ****************************************************************************/

void ll2xy (double lat,double lng, double *x_ptr,double *y_ptr)
{

    double       xi, eta;
    double       conf_lat, slat;

    slat      =    sin (lat);

    conf_lat  =    lat - slat * cos(lat) * (A0 + B0 * slat * slat + C0 * pow(slat, 4) 
                                            + D0 * pow(slat, 6));
   
    xi        =    atan (tan (conf_lat) / cos (lng - RefMeridian));
    
    eta       =    cos (conf_lat) * sin (lng - RefMeridian);
    eta       =    0.5 * log ((1 + eta) / (1 - eta));
    
    *y_ptr    =    Mean_Radius * (xi + beta1 * sin(2 * xi) * cosh(2 * eta) 
                               + beta2 * sin(4 * xi) * cosh(4 * eta) 
                               + beta3 * sin(6 * xi) * cosh(6 * eta) 
                               + beta4 * sin(8 * xi) * cosh(8 * eta));

    *x_ptr    =    Mean_Radius * (eta + beta1 * cos(2 * xi) * sinh(2 * eta)
                               + beta2 * cos(4 * xi) * sinh(4 * eta)
                               + beta3 * cos(6 * xi) * sinh(6 * eta)
                               + beta4 * cos(8 * xi) * sinh(8 * eta));

}

/****************************************************************************
* xy2ll

 Description: lat     calculated latitude in radians
              long    calculated longitude in radians
              x,y     grid coordinates

    Function: Calculates geographical coordinates lat, long in
              radians from grid coordinates x,y


 ****************************************************************************/

static void  xy2ll  (double *lat_ptr, double *long_ptr, double x, double y)
{
    double      xi, eta;
    double      conf_lat, slat;

    xi      =    y / Mean_Radius;
    eta     =    x / Mean_Radius;

    xi = xi - delta1 * sin(2 * xi) * cosh(2 * eta)
              - delta2 * sin(4 * xi) * cosh(4 * eta)
              - delta3 * sin(6 * xi) * cosh(6 * eta)
              - delta4 * sin(8 * xi) * cosh(8 * eta);

    eta = eta - delta1 * cos(2 * xi) * sinh(2 * eta)
                - delta2 * cos(4 * xi) * sinh(4 * eta)
                - delta3 * cos(6 * xi) * sinh(6 * eta)
                - delta4 * cos(8 * xi) * sinh(8 * eta);

    conf_lat  =    asin (sin (xi) / cosh (eta));

    *long_ptr =    atan (sinh (eta) / cos (xi)) + RefMeridian;

    slat = sin (conf_lat);

    *lat_ptr  =    conf_lat + slat * cos (conf_lat) * (A1 + B1 * slat * slat 
                              + C1 * pow(slat, 4) + D1 * pow(slat, 6));

}


/*****************************************************************************
* Origo_ll2xy

 Description: lat     latitude in radians
              long    longitude in radians
              x,y     calculated grid coordinates


    Function: 

 ****************************************************************************/

static void  Origo_ll2xy  (   double  lat,double lng,double *x_ptr, double *y_ptr)
{
  ll2xy(lat,lng,x_ptr,y_ptr);       /*  convert to x,y */
  *x_ptr = *x_ptr - origo_x;          /*  shift relative origin */
  *y_ptr = *y_ptr - origo_y;
}
/*****************************************************************************
* Origo_xy2ll

 Description: lat     calculated latitude in radians
              long    calculated longitude in radians
              x,y     grid coordinates

    Function: 
        
 ****************************************************************************/

static void  Origo_xy2ll  (double *lat,double *lng, double  x,double y)
{
  x = x + origo_x; /*  shift relative origin  */
  y = y + origo_y;
  xy2ll(lat,lng,x,y);   /*  convert to lat,long */
}


/*****************************************************************************
* SetRefMeridian

 Description: rad        Reference meridian in radians


    Function: 

 ****************************************************************************/

void  SetRefMeridian (double rad)
{
  RefMeridian = rad; 
  SetOrigo(origo_lat,origo_long);       /*  calc. new origin */
}


/*{****************************************************************************
* GetRefMeridian
 Description: rad        Referensmeridian angiven i radianer


    Function: Anv�ndes f�r att h�mta aktuell referensmeridian f�r koordinat-
              systemet.

 ****************************************************************************/

double  GetRefMeridian()
{
  return(RefMeridian);
}


/* ****************************************************************************
* SetOrigo
 Description: lat     latitude of origin in radians
              long    longitude of origin in radians


    Function:

 ****************************************************************************/

static void  SetOrigo (double lat,double lng)
{
  origo_lat = lat;
  origo_long = lng;
  ll2xy(lat,lng,&origo_x,&origo_y);
}

/*****************************************************************************
* GetOrigo
 Description: lat     latitude of origin in radians
              long    longitude of origin in radians
              x       x-coordinate for origin
              y       y-coordinate for origin


    Function: Get lat,long and x,y values for def. origin

 ****************************************************************************/

static void GetOrigo (double *lat, double* lng,double *x_ptr,double *y_ptr)
{
  *lat  = origo_lat;
  *lng = origo_long;
  *x_ptr    = origo_x;
  *y_ptr    = origo_y;
}


/*****************************************************************************
* SetEarthRadie
 Description: meter        Earth mean radius in meters


    Function: 

 ****************************************************************************/
static void  SetEarthRadie ( double meter)
{
  Mean_Radius = meter;
}
      

/*
* MeridianConv in degrees
 */
float MeridianConv(char *ll)
{
  SGeo geo;
  double diff_ref;
  double meco;

  if(!SlLlToGeo(ll,&geo)) printf("Merid: illegal\n");
  diff_ref =  RefMeridian - geo.lng;
  meco = Rad2Deg * diff_ref * sin (geo.lat);
  return(meco);    
  
}

/*{****************************************************************************
* GaussInit 
 *                           INITIATION ROUTINE                             *
 *                                                                          *
 ****************************************************************************/
static int GaussInit(SGeo geo) /* system origin */
{

    double  f, n, a_hat, k_0, scale;
    
    origo_x     = 0;
    origo_y     = 0;
    origo_lat   = 0;
    origo_long  = 0;
    
    //  Reference : "Gauss Conformal Projection (Transverse Mercator), 
    //                           Krugers formulae. Bo-Gunnar Reit", 2001-12-10,
    //                           National Land Survey of Sweden, Geodetic research
    //                           division, SE-80182 Gavle Sweden. 
    //                           http://www.lm.se/geodesi/
 
    //  WGS 84 Earth model used. MRT uses these 3 parameters: 

    earth_A   =   6378137.0;   //  equatorial radius (WGS 84)
    earth_B   =   6356752.3;   //  polar radius (WGS 84)  
    k_0       =   1.0;         //  map global scale factor

    //  Derived quantities

    f  =  (earth_A - earth_B) / earth_A;  // flattening

    n  =  f / (2.0 - f);

    a_hat  =  earth_A * (1 + n * n / 4.0 + n * n * n * n / 64.0 ) / (1 + n); // Earth mean radius 

    //  Calculate beta (to 4th order)

    beta1  =  n / 2.0 - (2.0 / 3.0) * n * n + (5.0 / 16.0) * n * n * n + (41.0 / 180.0) * n * n * n * n; // + h.o. 
    beta2  =  (13.0 / 48.0) * n * n - (3.0 / 5.0) * n * n * n + (557.0 / 1440.0) * n * n * n * n; // + h.o.
    beta3  =  (61.0 / 240) * n * n * n - (103.0 / 140.0) * n * n * n * n; // + h.o.
    beta4  =  (49561.0 / 161280.0) * n * n * n * n; // + h.o. 

    //  Calculate delta (to 4th order)

    delta1  =  n / 2.0 - (2.0 / 3.0) * n * n  + (37.0 / 96.0) * n * n * n - n * n * n * n / 360.0; // + h.o. 
    delta2  =  n * n / 48.0 + n * n * n / 15.0 - (437.0 / 1440.0) * n * n * n * n; // + h.o.
    delta3  =  (17.0 / 480.0) * n * n * n - (37.0 / 840.0) * n * n * n * n; // + h.o.
    delta4  =  (4397.0 / 161280.0) * n * n * n * n; // + h.o.

    //  Define map scale

    scale   =  k_0 * a_hat; 

    //  Calculate Eccentricity squared, and powers of it

    ecc2  =  f * (2.0 - f); 

    ecc   =  sqrt (ecc2);
    ecc4  =  ecc * ecc * ecc * ecc;
    ecc6  =  ecc * ecc * ecc * ecc * ecc * ecc;
    ecc8  =  ecc * ecc * ecc * ecc * ecc * ecc * ecc * ecc;

    A0 = ecc2; 
    B0 = (5 * ecc4 - ecc6) / 6.0;
    C0 = (104 * ecc6 - 45 * ecc8) / 120.0;
    D0 = (1237 * ecc8) / 1260.0;

    A1 = ecc2 + ecc4 + ecc6 + ecc8;
    B1 = - (7 * ecc4 + 17 * ecc6 + 30 * ecc8) / 6.0;
    C1 = (224 * ecc6 + 889 * ecc8) / 120.0;
    D1 = - (4279 * ecc8) / 1260.0; 

    RefMeridian  = geo.lng;

    Mean_Radius  = scale;

    SetOrigo(geo.lat,geo.lng);

    return (1);

}

/*********************************************************************
 *
* Function: SlGeoToXy
 *
 * Description: This function calculates lat long in radians to xy in
 *              world coordinates (m)
 *        
 ********************************************************************/
void SlGeoToXy(SGeo *geo_pos, SPos *poswc,int n)
{ 
  double x,y;
    

  if (n<=0) return;
  while (n--)
  {
    Origo_ll2xy(geo_pos->lat,geo_pos->lng,&x,&y);
    
    poswc->x = x;
    poswc->y = y;

    poswc++;
    geo_pos++;
  }
}

/*********************************************************************
 *
* Function: SlXyToGeo
 *
 * Description: This function calculates lat/long in radians for 
 *              the given xy.
 *
 * 
 ********************************************************************/
void SlXyToGeo(SPos *poswc,SGeo *geo_pos,int no)
{
  double lat,lng;
  
  if (no<=0) return;
  while (no--)
  {
    Origo_xy2ll(&lat,&lng,poswc->x,poswc->y);
    geo_pos->lat = lat;
    geo_pos->lng = lng;
    geo_pos++;
    poswc++;
  }
}


/*********************************************************************
 *
* Function: SlLlToGeo
 *
 * Description: Convert lat/long string to radians 
 *              format e.g 'ddmmssNdddmmssE'
 ********************************************************************/
int SlLlToGeo(char *posll,SGeo *geo)
{
  int dx,mx,sx,dy,my,sy,i;
  char ns,we;
  double lat_degr,long_degr;
  char buff[40],*p;

  
  sscanf(posll,"%20s",buff);  /* skip leading space and cr */
  p = buff;
  i = 0;
  while (isdigit(p[i])) i++;
  if ((i != 6) || (p[i] !='N' && p[i] != 'S')) return(0);
  i++;
  while (isdigit(p[i])) i++;
  if ((i != 14) || (p[i] !='W' && p[i] != 'E')) return(0);

  if((sscanf(posll,
    "%2d%2d%2d%c%3d%2d%2d%c",&dy,&my,&sy,&ns,&dx,&mx,&sx,&we))!=8)
  return(0);
  if ((ns == 'N' ||  ns == 'S') &&
        (we == 'E' ||  we == 'W'))
  {
      if (!( dy < 90 && my<60  && sy <60 && dx < 180 && mx < 60 &&
          sx < 60))
      return(0);   /* False */
  }  else return(0);   /* False */
    
  lat_degr = dy + ((double)my)/60.0 + ((double)sy)/3600.0; /* make degrees */
  if (ns == 'S') lat_degr = -lat_degr;
  geo->lat = lat_degr * Deg2Rad;                                 /* make radians */

  long_degr = dx + ((double)mx)/60.0 + ((double)sx)/3600.0; /* make degrees */
  if (we == 'W') long_degr = -long_degr;
  geo->lng = long_degr * Deg2Rad;                                /* make radians */
 
  
  return(1);   /* True */
}

/*********************************************************************
 *
* Function: SlLlToGeo2
 *
 * Description: Convert lat/long string to radians  (100th secs)
 *              format e.g 'ddmmssssNdddmmssssE'
 ********************************************************************/
int SlLlToGeo2(char *posll,SGeo *geo)
{
  int dx,mx,sx,dy,my,sy, ssx, ssy, i;
  char ns,we;
  double lat_degr,long_degr;
  char buff[40],*p;

  sscanf(posll,"%20s",buff);  /* skip leading space and cr */
  p = buff;
  i = 0;
  while (isdigit(p[i])) i++;
  if ((i != 8) || (p[i] !='N' && p[i] != 'S')) return(0);
  i++;
  while (isdigit(p[i])) i++;
  if ((i != 18) || (p[i] !='W' && p[i] != 'E')) return(0);

  if((sscanf(posll,
    "%2d%2d%2d%2d%c%3d%2d%2d%2d%c",&dy,&my,&sy, &ssy, &ns,&dx,&mx,&sx,&ssx,&we))!= 10)
  return(0);
  if ((ns == 'N' ||  ns == 'S') &&
        (we == 'E' ||  we == 'W'))
  {
      if (!( dy < 90 && my<60  && sy <60 && dx < 180 && mx < 60 &&
          sx < 60))
      return(0);   /* False */
  }  else return(0);   /* False */
    
  lat_degr = dy + ((double)my)/60.0 + ((double)sy)/3600.0 + ((double)ssy)/360000.0; /* make degrees */
  if (ns == 'S') lat_degr = -lat_degr;
  geo->lat = lat_degr * Deg2Rad;                                 /* make radians */

  long_degr = dx + ((double)mx)/60.0 + ((double)sx)/3600.0 + ((double)ssx)/360000.0; /* make degrees */
  if (we == 'W') long_degr = -long_degr;
  geo->lng = long_degr * Deg2Rad;                                /* make radians */
 
  
  return(1);   /* True */
}

/*********************************************************************
 *
* Function: SlGeoToLl
 *
 * Description: Convert lat/long in radians to lat/long  string
 * 
 ********************************************************************/
int SlGeoToLl(SGeo *geo, char *posll)
{
  double x,y;
  int secs, dx,mx,sx,dy,my,sy;
  char ns='N',we='E';

  x = geo->lng;
  y = geo->lat;
  if (x<0) { x = -x; we = 'W';}
  if (y<0) { y = -y; ns = 'S';}
 
  x *= Rad2Deg;
  y *= Rad2Deg;
  
  secs = (int)(y*3600.0 + 0.5); /* total y in seconds */
  dy = secs/3600;
  my = secs/60 - dy*60; 
  sy = secs - my*60 - dy*3600;

  secs = (int)(x*3600.0 + 0.5); /* total x in seconds */
  dx = secs/3600;
  mx = secs/60 - dx*60; 
  sx = secs - mx*60 - dx*3600;
  
  sprintf(posll,"%02d%02d%02d%c%03d%02d%02d%c",dy,my,sy,ns,dx,mx,sx,we);

  return (1); /* True */
}

/*********************************************************************
 *
* Function: SlGeoToLl2
 *
 * Description: Convert lat/long in radians to lat/long  string (100th secs)
 * 
 ********************************************************************/
int SlGeoToLl2(SGeo *geo, char *posll)
{
  double x,y;
  int fsecs, secs, dx,mx,sx,dy,my,sy, ssx, ssy;
  char ns='N',we='E';
  double  pi = 4*atan((double) 1.0);  
  x = geo->lng;
  y = geo->lat;
  if (x<0) { x = -x; we = 'W';}
  if (y<0) { y = -y; ns = 'S';}
 
  x *= 180.0/pi;
  y *= 180.0/pi;
  
  ssy = (int)(y*360000.0 + 0.5);
  secs = ssy/100;
  ssy = (int)(100 * (y*3600.0 - secs) );   
  dy = secs/3600;
  my = secs/60 - dy*60; 
  sy = secs - my*60 - dy*3600;

  ssx = (int)(x*360000.0 + 0.5);
  secs = ssx/100;
  ssx = (int)(100 * (x*3600.0 - secs) );   
  dx = secs/3600;
  mx = secs/60 - dx*60; 
  sx = secs - mx*60 - dx*3600;
  
  sprintf(posll,"%02d%02d%02d%02d%c%03d%02d%02d%02d%c",dy,my,sy,ssy,ns,dx,mx,sx,ssx,we);

  return (1); /* True */
}

/*********************************************************************
 *
* Function: SlInitiation
 *
 * Description: This function initiates the slant projection
 *              constants
 *
 ********************************************************************/
int SlInitiation(char *posll)
{ 
  SGeo geo;
 
  if(!SlLlToGeo(posll,&geo))
  {
    printf("illegal lat/long");
    return(0);
  }
  GaussInit(geo);

  return(1);
}

/*********************************************************************
 *
* Function: SlXyToATCASMeco
 *
 * Description: This function calculates the meridian convergence in radians
 *              for given xy for the standard meridian
 *
 ********************************************************************/
void SlXyToATCASMeco(SPos  poswc,float *meco_ptr)
{
  double diff_ref;
  double meco;
  double lat,lng;

  Origo_xy2ll(&lat,&lng,poswc.x,poswc.y);
  diff_ref =  RefMeridian - lng;
  meco = diff_ref * sin (lat);
  *meco_ptr = meco;
  
}
/*********************************************************************
 *
* Function: SlXyToLatLongText
 *
 * Description: This function calculates the meridian convergence in radians
 *              for given xy for the standard meridian
 *
 ********************************************************************/

char *SlXyToLatLongText   (SPos                         xyPos,
                           char                        *pointLatLongText)

{

    SGeo            geoPos;
                    
    //  Coordinates to radians

    SlXyToGeo (&xyPos,
               &geoPos,
               1);

    //  Radians to lat/long string

    SlGeoToLl (&geoPos,
               pointLatLongText);

    return (pointLatLongText);
    
}

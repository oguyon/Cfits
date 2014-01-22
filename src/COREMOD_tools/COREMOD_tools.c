#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "../Cfits.h"
#include "../00CORE/00CORE.h"

#define SBUFFERSIZE 1000

char errormessage[SBUFFERSIZE];



int create_counter_file(char *fname, long NBpts)
{
  long i;
  FILE *fp;
  
  if((fp=fopen(fname,"w"))==NULL)
    {
      sprintf(errormessage,"cannot create file \"%s\"",fname);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }

  for(i=0;i<NBpts;i++)
    fprintf(fp,"%ld %f\n",i,(PRECISION) (1.0*i/NBpts));
  
  fclose(fp);

  return(0);
}

int bubble_sort(PRECISION *array, long count)
{
  register long a,b;
  register PRECISION t;

  for(a=1; a<count; a++)
    for(b=count-1; b>=a; b--)
      if(array[b-1]>array[b])
	{
	  t = array[b-1];
	  array[b-1] = array[b];
	  array[b] = t;
	}
  return(0);
}

void qs(PRECISION *array, long left, long right)
{
  register long i,j;
  PRECISION x,y;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;
      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs(array,left,j);
  if(i<right) qs(array,i,right);
}

void qs_long(long *array, long left, long right)
{
  register long i,j;
  long x,y;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;
      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs_long(array,left,j);
  if(i<right) qs_long(array,i,right);
}

void qs_double(double *array, long left, long right)
{
  register long i,j;
  double x,y;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;
      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs_double(array,left,j);
  if(i<right) qs_double(array,i,right);
}

void qs3(PRECISION *array, PRECISION *array1, PRECISION *array2, long left, long right)
{
  register long i,j;
  PRECISION x,y;
  PRECISION y1,y2;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      y1 = array1[i];
      array1[i] = array1[j];
      array1[j] = y1;

      y2 = array2[i];
      array2[i] = array2[j];
      array2[j] = y2;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs3(array,array1,array2,left,j);
  if(i<right) qs3(array,array1,array2,i,right);
}



void qs3_float(float *array, float *array1, float *array2, long left, long right)
{
  register long i,j;
  float x,y;
  float y1,y2;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      y1 = array1[i];
      array1[i] = array1[j];
      array1[j] = y1;

      y2 = array2[i];
      array2[i] = array2[j];
      array2[j] = y2;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs3_float(array,array1,array2,left,j);
  if(i<right) qs3_float(array,array1,array2,i,right);
}

void qs3_double(double *array, double *array1, double *array2, long left, long right)
{
  register long i,j;
  double x,y;
  double y1,y2;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      y1 = array1[i];
      array1[i] = array1[j];
      array1[j] = y1;

      y2 = array2[i];
      array2[i] = array2[j];
      array2[j] = y2;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs3_double(array,array1,array2,left,j);
  if(i<right) qs3_double(array,array1,array2,i,right);
}


void qs2l(PRECISION *array, long *array1, long left, long right)
{
  register long i,j;
  PRECISION x,y;
  long l1;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      l1 = array1[i];
      array1[i] = array1[j];
      array1[j] = l1;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs2l(array,array1,left,j);
  if(i<right) qs2l(array,array1,i,right);
}

void qs2l_double(double *array, long *array1, long left, long right)
{
  register long i,j;
  double x,y;
  long l1;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      l1 = array1[i];
      array1[i] = array1[j];
      array1[j] = l1;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs2l_double(array,array1,left,j);
  if(i<right) qs2l_double(array,array1,i,right);
}

void qs3ll_double(double *array, long *array1, long *array2, long left, long right)
{
  register long i,j;
  double x,y;
  long l1,l2;
  
  i = left; j = right;
  x = array[(left+right)/2];
  
  do {
    while(array[i]<x && i<right) i++;
    while(x<array[j] && j>left) j--;
    
    if(i<=j) {
      y = array[i];
      array[i] = array[j];
      array[j] = y;

      l1 = array1[i];
      array1[i] = array1[j];
      array1[j] = l1;

      l2 = array2[i];
      array2[i] = array2[j];
      array2[j] = l2;

      i++; j--;
    }
  } while(i<=j);
  
  if(left<j) qs3ll_double(array,array1,array2,left,j);
  if(i<right) qs3ll_double(array,array1,array2,i,right);
}

void quick_sort(PRECISION *array, long count)
{
  qs(array,0,count-1);
}

void quick_sort_long(long *array, long count)
{
  qs_long(array,0,count-1);
}

void quick_sort_double(double *array, long count)
{
  qs_double(array,0,count-1);
}

void quick_sort3(PRECISION *array, PRECISION *array1, PRECISION *array2, long count)
{
  qs3(array,array1,array2,0,count-1);
}

void quick_sort3_float(float *array, float *array1, float *array2, long count)
{
  qs3_float(array,array1,array2,0,count-1);
}

void quick_sort3_double(double *array, double *array1, double *array2, long count)
{
  qs3_double(array,array1,array2,0,count-1);
}

void quick_sort2l(PRECISION *array, long *array1, long count)
{
  qs2l(array,array1,0,count-1);
}

void quick_sort2l_double(double *array, long *array1, long count)
{
  qs2l_double(array,array1,0,count-1);
}

void quick_sort3ll_double(double *array, long *array1, long *array2, long count)
{
  qs3ll_double(array,array1,array2,0,count-1);
}

int lin_regress(PRECISION *a, PRECISION *b, PRECISION *Xi2, PRECISION *x, PRECISION *y, PRECISION *sig, int nb_points)
{
  PRECISION S,Sx,Sy,Sxx,Sxy,Syy;
  int i;
  PRECISION delta;

  S = 0;
  Sx = 0;
  Sy = 0;
  Sxx = 0;
  Syy = 0;
  Sxy = 0;
  for(i=0;i<nb_points;i++)
    {
      S += 1.0/sig[i]/sig[i];
      Sx += x[i]/sig[i]/sig[i];
      Sy += y[i]/sig[i]/sig[i];
      Sxx += x[i]*x[i]/sig[i]/sig[i];
      Syy += y[i]*y[i]/sig[i]/sig[i];
      Sxy += x[i]*y[i]/sig[i]/sig[i];
    }

  delta = S*Sxx-Sx*Sx;
  *a = (Sxx*Sy-Sx*Sxy)/delta;
  *b = (S*Sxy-Sx*Sy)/delta;
  *Xi2 = Syy-2*(*a)*Sy-2*(*a)*(*b)*Sx+(*a)*(*a)*S+2*(*a)*(*b)*Sx-(*b)*(*b)*Sxx;
  
  return(0);
}

int replace_char(char *content, char cin, char cout)
{
  long i;

  for(i=0;i<strlen(content);i++)
    if(content[i]==cin)
      content[i] = cout;

  return(0);
}


int read_config_parameter_exists(char *config_file, char *keyword)
{
  FILE *fp;
  char line[1000];
  char keyw[200];
  char cont[200];
  int read;

  read = 0;
  if((fp=fopen(config_file,"r"))==NULL)
    {
      sprintf(errormessage,"cannot open file \"%s\"",config_file);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  while((fgets(line,1000,fp)!=NULL)&&(read==0))
    {
      sscanf(line,"%s",keyw);
      if(strcmp(keyw,keyword)==0)
	read = 1;
    }
  if(read==0)
    {
      sprintf(errormessage,"parameter \"%s\" does not exist in file \"%s\"",keyword,config_file);
      printWARNING(__FILE__,__func__,__LINE__,errormessage);
    }

  fclose(fp);
    
  return(read);
}


int read_config_parameter(char *config_file, char *keyword, char *content)
{
  FILE *fp;
  char line[1000];
  char keyw[200];
  char cont[200];
  int read;

  read = 0;
  if((fp=fopen(config_file,"r"))==NULL)
    {
      sprintf(errormessage,"cannot open file \"%s\"",config_file);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  strcpy(content,"---");
  while(fgets(line,1000,fp)!=NULL)
    {
      sscanf(line,"%s %s",keyw,cont);
      if(strcmp(keyw,keyword)==0)
	{
	  strcpy(content,cont);
	  read = 1;
	}
      /*      printf("KEYWORD : \"%s\"   CONTENT : \"%s\"\n",keyw,cont);*/
    }
  if(read==0)
    {
      sprintf(errormessage,"parameter \"%s\" does not exist in file \"%s\"",keyword,config_file);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      sprintf(content,"-");
      //  exit(0);
    }

  fclose(fp);
  
  return(read);
}

float read_config_parameter_float(char *config_file, char *keyword)
{
  float value;
  char content[SBUFFERSIZE];
  
  read_config_parameter(config_file,keyword,content);
  //printf("content = \"%s\"\n",content);
  value = atof(content);
  //printf("Value = %g\n",value);

  return(value);
}

long read_config_parameter_long(char *config_file, char *keyword)
{
  long value;
  char content[SBUFFERSIZE];
  
  read_config_parameter(config_file,keyword,content);
  value = atol(content);
  
  return(value);
}

long read_config_parameter_int(char *config_file, char *keyword)
{
  int value;
  char content[SBUFFERSIZE];
  
  read_config_parameter(config_file,keyword,content);
  value = atoi(content);
  
  return(value);
}




long file_number_lines(char *file_name)
{
  long cnt;
  int c;
  FILE *fp;

  if((fp=fopen(file_name,"r"))==NULL)
    {
      sprintf(errormessage,"cannot open file \"%s\"",file_name);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  cnt = 0;
  while((c=fgetc(fp))!=EOF)
    if(c=='\n')
      cnt++;
  fclose(fp);

  return(cnt);
}

FILE* open_file_w(char *filename)
{
  FILE *fp;

  if((fp=fopen(filename,"w"))==NULL)
    {
      sprintf(errormessage,"cannot create file \"%s\"",filename);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  return(fp);
}

FILE* open_file_r(char *filename)
{
  FILE *fp;

  if((fp=fopen(filename,"r"))==NULL)
    {
      sprintf(errormessage,"cannot read file \"%s\"",filename);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  return(fp);
}

int write_1D_array(PRECISION *array, long nbpoints, char *filename)
{
  FILE *fp;
  long ii;

  fp = open_file_w(filename);
  for(ii=0;ii<nbpoints;ii++)
    fprintf(fp,"%ld\t%f\n",ii,array[ii]);
  fclose(fp);
  
  return(0);
}

int read_1D_array(PRECISION *array, long nbpoints, char *filename)
{
  FILE *fp;
  long ii;
  long tmpl;
  
  fp = open_file_r(filename);
  for(ii=0;ii<nbpoints;ii++)
    {
      if(fscanf(fp,"%ld\t%f\n",&tmpl,&array[ii])!=2)
	{
	  printERROR(__FILE__,__func__,__LINE__,"fscanf error");
	  exit(0);
	}
    }
  fclose(fp);
  
  return(0);
}



/* test point */ 
int tp(char *word)
{
  printf("---- Test point %s ----\n",word);
  fflush(stdout);

  return(0);
}

int read_int_file(char *fname)
{
  int value;
  FILE *fp;
  
  if((fp = fopen(fname,"r"))==NULL)
    {
      value = 0;
    }
  else
    {
      if(fscanf(fp,"%d",&value)!=1)
	{
	  printERROR(__FILE__,__func__,__LINE__,"fscanf error");
	  exit(0);
	}
      fclose(fp);
    }

  return(value);
}

int write_int_file(char *fname, int value)
{
  FILE *fp;
  
  if((fp = fopen(fname,"w"))==NULL)
    {
      sprintf(errormessage,"cannot create file \"%s\"\n",fname);
      printERROR(__FILE__,__func__,__LINE__,errormessage);
      exit(0);
    }
  
  fprintf(fp,"%d\n",value);
  fclose(fp);

  return(value);
}

int write_float_file(char *fname, float value)
{
  FILE *fp;
  int mode = 0; // default, create single file
  
  if(variable_ID("WRITE2FILE_APPEND")!=-1)
    mode = 1;

  if(mode == 0)
    {
      if((fp = fopen(fname,"w"))==NULL)
	{
	  sprintf(errormessage,"cannot create file \"%s\"\n",fname);
	  printERROR(__FILE__,__func__,__LINE__,errormessage);
	  exit(0);
	}
      fprintf(fp,"%g\n",value);
      fclose(fp);
    }
  
  if(mode == 1)
   {
      if((fp = fopen(fname,"a"))==NULL)
	{
	  sprintf(errormessage,"cannot create file \"%s\"\n",fname);
	  printERROR(__FILE__,__func__,__LINE__,errormessage);
	  exit(0);
	}
      fprintf(fp," %g",value);
      fclose(fp);
   }

  return(0);
}

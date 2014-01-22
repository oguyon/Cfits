#include <fitsio.h>  /* required by every program that uses CFITSIO  */
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "../Cfits.h"
#include "../00CORE/00CORE.h"
#include "../COREMOD_tools/COREMOD_tools.h"
#include "../COREMOD_memory/COREMOD_memory.h"
#include "../COREMOD_arith/COREMOD_arith.h"
#include "../COREMOD_fft/COREMOD_fft.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

extern DATA data;

long brighter(char *ID_name, PRECISION value) /* number of pixels brighter than value */
{
  int ID;
  long ii,jj;
  long naxes[2];
  long brighter, fainter;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
    
  brighter = 0;
  fainter = 0;
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++){
      if(data.image[ID].arrayD[jj*naxes[0]+ii]>value)
	brighter++;
      else
	fainter++;
    }
  printf("brighter %ld   fainter %ld\n", brighter, fainter );

  return(brighter);
}

int img_nbpix_flux(char *ID_name)
{
  int ID;
  long ii,jj;
  long naxes[2];
  PRECISION value = 0;
  PRECISION *array;
  long nelements,i;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
  nelements = naxes[0]*naxes[1];
  
  array = (PRECISION*) malloc(naxes[1]*naxes[0]*sizeof(PRECISION));
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++)
      array[jj*naxes[0]+ii] = data.image[ID].arrayD[jj*naxes[0]+ii];

  quick_sort(array,nelements);
  
  for(i=0;i<nelements;i++)
    {
      value += array[i];
      printf("%ld  %20.18e\n",i,value);
    }

  free(array);
  return(0);
}

PRECISION img_percentile(char *ID_name, PRECISION p)
{
  int ID;
  long ii;
  long naxes[2];
  PRECISION value = 0;
  PRECISION *array;
  long nelements;
  long n;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
  nelements = naxes[0]*naxes[1];
  
  array = (PRECISION*) malloc(nelements*sizeof(PRECISION));
  for (ii = 0; ii < nelements; ii++) 
    array[ii] = data.image[ID].arrayD[ii];

  quick_sort(array,nelements);
  
  n = (long) (p*naxes[1]*naxes[0]);
  if(n>(nelements-1))
    n = (nelements-1);
  if(n<0)
    n = 0;
  value = array[n];
  /*  printf("1  percent  %20.18e\n",array[(long) (0.01*naxes[1]*naxes[0])]);
  printf("5  percent  %20.18e\n",array[(long) (0.05*naxes[1]*naxes[0])]);
  printf("10 percent  %20.18e\n",array[(long) (0.1*naxes[1]*naxes[0])]);
  printf("20 percent  %20.18e\n",array[(long) (0.2*naxes[1]*naxes[0])]);
  printf("50 percent  %20.18e\n",array[(long) (0.5*naxes[1]*naxes[0])]);
  printf("80 percent  %20.18e\n",array[(long) (0.8*naxes[1]*naxes[0])]);
  printf("90 percent  %20.18e\n",array[(long) (0.9*naxes[1]*naxes[0])]);
  printf("95 percent  %20.18e\n",array[(long) (0.95*naxes[1]*naxes[0])]);
  printf("99 percent  %20.18e\n",array[(long) (0.99*naxes[1]*naxes[0])]);*/
  free(array);
  return(value);
}

int img_histoc(char *ID_name, char *fname)
{
  FILE *fp;
  int ID;
  long ii,jj;
  long naxes[2];
  PRECISION value = 0;
  PRECISION *array;
  long nelements;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
  nelements = naxes[0]*naxes[1];
  
  array = (PRECISION*) malloc(naxes[1]*naxes[0]*sizeof(PRECISION));
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++)
      array[jj*naxes[0]+ii] = data.image[ID].arrayD[jj*naxes[0]+ii];

  quick_sort(array,nelements);
  
  if((fp=fopen(fname,"w"))==NULL)
    {
      printf("ERROR: cannot open file \"%s\"\n",fname);
      exit(0);
    }
  value = 0.0;
  for(ii=0;ii<nelements;ii++)
    {
      value += array[ii];
      if(ii>0.99*nelements)
	fprintf(fp,"%ld %g %g\n",nelements-ii,value,array[ii]);
    }

  fclose(fp);
  free(array);

  return(0);
}

int make_histogram(char *ID_name, char *ID_out_name, PRECISION min, PRECISION max, long nbsteps)
{
  int ID,ID_out;
  long ii,jj;
  long naxes[2];
  long n;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    

  create_2Dimage_ID(ID_out_name,nbsteps,1);
  ID_out = image_ID(ID_out_name);
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++)
      {
	n = (long) ((data.image[ID].arrayD[jj*naxes[0]+ii]-min)/(max-min)*nbsteps);
	if((n>0)&&(n<nbsteps))
	  data.image[ID_out].arrayD[n] += 1;
      }
  return(0);
}

PRECISION ssquare(char *ID_name)
{
  int ID;
  long ii,jj;
  long naxes[2];
  PRECISION ssquare;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
    
  ssquare = 0;
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++){
      ssquare = ssquare + data.image[ID].arrayD[jj*naxes[0]+ii]*data.image[ID].arrayD[jj*naxes[0]+ii];
    }
  return(ssquare);
}

PRECISION rms_dev(char *ID_name)
{
  int ID;
  long ii,jj;
  long naxes[2];
  PRECISION ssquare,rms;
  PRECISION constant;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
    
  ssquare = 0;
  constant = arith_image_total(ID_name)/naxes[0]/naxes[1];
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++){
      ssquare = ssquare + (data.image[ID].arrayD[jj*naxes[0]+ii]-constant)*(data.image[ID].arrayD[jj*naxes[0]+ii]-constant);
    }
  rms = sqrt(ssquare/naxes[1]/naxes[0]);
  return(rms);
}


// option "fileout" : output to file imstat.info.txt
int stats(char *ID_name, char *options)
{
  int ID;
  long ii,jj,j;
  PRECISION min,max;
  double rms;
  long nelements;
  double tot;
  PRECISION *array;
  long iimin,iimax;
  int atype;
  long tmp_long;
  char type[20];
  char vname[200];
  double xtot,ytot;
  double vbx,vby;
  FILE *fp;
  int mode = 0;
  
  // printf("OPTIONS = %s\n",options);
   if (strstr(options,"fileout")!=NULL)
    mode = 1;
   
  if(mode == 1)
    fp = fopen("imstat.info.txt","w");
  

  ID = image_ID_noaccessupdate(ID_name);
  if(ID!=-1)
    {
      nelements =  data.image[ID].nelement;
      
      atype = data.image[ID].atype;
      tmp_long = data.image[ID].nelement*TYPESIZE[atype];
      printf("\n");
      printf("Image size (->imsize0...):     [");
      printf("% ld",data.image[ID].size[0]);
      j = 0;
      sprintf(vname,"imsize%ld",j);
      create_variable_ID(vname,1.0*data.image[ID].size[j]);
      for(j=1;j<data.image[ID].naxis;j++)
	{
	  printf(" %ld",data.image[ID].size[j]);
	  sprintf(vname,"imsize%ld",j);
	  create_variable_ID(vname,1.0*data.image[ID].size[j]);
	}
      printf(" ]\n");

      if(atype==CHAR)
	sprintf(type,"CHAR");
      if(atype==INT)
	sprintf(type,"INT");
      if(atype==FLOAT)
	sprintf(type,"FLOAT");
      if(atype==DOUBLE)
	sprintf(type,"DOUBLE");
      if(atype==COMPLEX_FLOAT)
	sprintf(type,"CFLOAT");
      if(atype==COMPLEX_DOUBLE)
	sprintf(type,"CDOUBLE");
      printf("type:            %s\n",type);
      printf("Memory size:     %ld Kb\n",(long) tmp_long/1024);
      printf("Created:         %s",asctime(localtime(&data.image[ID].creation_time)));
      printf("Last access:     %s\n",asctime(localtime(&data.image[ID].last_access)));



      min = data.image[ID].arrayD[0];
      max = data.image[ID].arrayD[0];
      
      iimin = 0;
      iimax = 0;
      for (ii = 0; ii < nelements; ii++) 
	{
	  if (min > data.image[ID].arrayD[ii])
	    {
	      min = data.image[ID].arrayD[ii];
	      iimin = ii;
	    }
	  if (max < data.image[ID].arrayD[ii])
	    {	
	      max = data.image[ID].arrayD[ii];
	      iimax = ii;
	    }
	}

      array = (PRECISION*) malloc(nelements*sizeof(PRECISION));
      tot = 0.0;

      rms = 0.0;
      for (ii = 0; ii < nelements; ii++) 
	{
	  tot += data.image[ID].arrayD[ii];
	  rms += data.image[ID].arrayD[ii]*data.image[ID].arrayD[ii];
	  array[ii] = data.image[ID].arrayD[ii];
	}
      rms = sqrt(rms);

      printf("minimum         (->vmin)     %20.18e [ pix %ld ]\n", min, iimin);
      if(mode == 1)
	fprintf(fp,"minimum                  %20.18e [ pix %ld ]\n", min, iimin);
      create_variable_ID("vmin",min);
      printf("maximum         (->vmax)     %20.18e [ pix %ld ]\n", max, iimax);
      if(mode == 1)
	fprintf(fp,"maximum                  %20.18e [ pix %ld ]\n", max, iimax);
      create_variable_ID("vmax",max);
      printf("total           (->vtot)     %20.18e\n",tot);
      if(mode == 1)
	fprintf(fp,"total                    %20.18e\n",tot);
      create_variable_ID("vtot",tot);
      printf("rms             (->vrms)     %20.18e\n",rms);
      if(mode == 1)
	fprintf(fp,"rms                      %20.18e\n",rms);
      create_variable_ID("vrms",rms);
      printf("rms per pixel   (->vrmsp)    %20.18e\n",rms/sqrt(nelements));
      if(mode == 1)
	fprintf(fp,"rms per pixel            %20.18e\n",rms/sqrt(nelements));
      create_variable_ID("vrmsp",rms/sqrt(nelements));
      printf("rms dev per pix (->vrmsdp)   %20.18e\n",sqrt(rms*rms/nelements - tot*tot/nelements/nelements));
      create_variable_ID("vrmsdp",sqrt(rms*rms/nelements - tot*tot/nelements/nelements));
      printf("mean            (->vmean)    %20.18e\n",tot/nelements);
      if(mode == 1)
	fprintf(fp,"mean                     %20.18e\n",tot/nelements);
      create_variable_ID("vmean",tot/nelements);
      
      if(data.image[ID].naxis==2)
	{
	  xtot = 0.0;
	  ytot = 0.0;
	  for(ii=0;ii<data.image[ID].size[0];ii++)
	    for(jj=0;jj<data.image[ID].size[1];jj++)
	      {
		xtot += data.image[ID].arrayD[jj*data.image[ID].size[0]+ii]*ii;
		ytot += data.image[ID].arrayD[jj*data.image[ID].size[0]+ii]*jj;
	      }
	  vbx = xtot/tot;
	  vby = ytot/tot;
	  printf("Barycenter x    (->vbx)      %20.18f\n",vbx);
	  if(mode == 1)
	    fprintf(fp,"photocenterX             %20.18e\n",vbx);
	  create_variable_ID("vbx",vbx);
	  printf("Barycenter y    (->vby)      %20.18f\n",vby);	  
	  if(mode == 1)
	    fprintf(fp,"photocenterY             %20.18e\n",vby);
	  create_variable_ID("vby",vby);
	}

      quick_sort(array,nelements);
      printf("\n");
      printf("percentile values:\n");
      
      printf("1  percent      (->vp01)     %20.18e\n",array[(long) (0.01*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile01             %20.18e\n",array[(long) (0.01*nelements)]);
      create_variable_ID("vp01",array[(long) (0.01*nelements)]);

      printf("5  percent      (->vp05)     %20.18e\n",array[(long) (0.05*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile05             %20.18e\n",array[(long) (0.05*nelements)]);
      create_variable_ID("vp05",array[(long) (0.05*nelements)]);

      printf("10 percent      (->vp10)     %20.18e\n",array[(long) (0.1*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile10             %20.18e\n",array[(long) (0.10*nelements)]);
      create_variable_ID("vp10",array[(long) (0.1*nelements)]);

      printf("20 percent      (->vp20)     %20.18e\n",array[(long) (0.2*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile20             %20.18e\n",array[(long) (0.20*nelements)]);
      create_variable_ID("vp20",array[(long) (0.2*nelements)]);

      printf("50 percent      (->vp50)     %20.18e\n",array[(long) (0.5*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile50             %20.18e\n",array[(long) (0.50*nelements)]);
      create_variable_ID("vp50",array[(long) (0.5*nelements)]);

      printf("80 percent      (->vp80)     %20.18e\n",array[(long) (0.8*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile80             %20.18e\n",array[(long) (0.80*nelements)]);
      create_variable_ID("vp80",array[(long) (0.8*nelements)]);

      printf("90 percent      (->vp90)     %20.18e\n",array[(long) (0.9*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile90             %20.18e\n",array[(long) (0.90*nelements)]);
      create_variable_ID("vp90",array[(long) (0.9*nelements)]);

      printf("95 percent      (->vp95)     %20.18e\n",array[(long) (0.95*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile95             %20.18e\n",array[(long) (0.95*nelements)]);
      create_variable_ID("vp95",array[(long) (0.95*nelements)]);

      printf("99 percent      (->vp99)     %20.18e\n",array[(long) (0.99*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile99             %20.18e\n",array[(long) (0.99*nelements)]);
      create_variable_ID("vp99",array[(long) (0.99*nelements)]);

      printf("99.5 percent    (->vp995)    %20.18e\n",array[(long) (0.995*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile995            %20.18e\n",array[(long) (0.995*nelements)]);
      create_variable_ID("vp995",array[(long) (0.995*nelements)]);

      printf("99.8 percent    (->vp998)    %20.18e\n",array[(long) (0.998*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile998            %20.18e\n",array[(long) (0.998*nelements)]);
      create_variable_ID("vp998",array[(long) (0.998*nelements)]);

      printf("99.9 percent    (->vp999)    %20.18e\n",array[(long) (0.999*nelements)]);
      if(mode == 1)
	fprintf(fp,"percentile999            %20.18e\n",array[(long) (0.999*nelements)]);
      create_variable_ID("vp999",array[(long) (0.999*nelements)]);

      printf("\n");
      free(array);      
    }
  
  if(mode == 1)
    fclose(fp);
  
  return(0);
}

PRECISION img_min(char *ID_name)
{
  int ID;
  long ii;
  PRECISION min;

  ID = image_ID(ID_name);

  min = data.image[ID].arrayD[0];
  for (ii = 0; ii < data.image[ID].nelement; ii++) 
    if (min > data.image[ID].arrayD[ii])
      min = data.image[ID].arrayD[ii];

   return(min);
}

PRECISION img_max(char *ID_name)
{
  int ID;
  long ii;
  PRECISION max;

  ID = image_ID(ID_name);

  max = data.image[ID].arrayD[0];
    for (ii = 0; ii < data.image[ID].nelement; ii++) 
      if (max < data.image[ID].arrayD[ii])
	max = data.image[ID].arrayD[ii];

   return(max);
}

int profile(char *ID_name, char *outfile, PRECISION xcenter, PRECISION ycenter, PRECISION step, long nb_step)
{
  int ID;
  long ii,jj;
  long naxes[2];
  long nelements;
  PRECISION distance;
  PRECISION *dist;
  PRECISION *mean;
  PRECISION *rms;
  long *counts;
  FILE *fp;
  long i;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
  nelements = naxes[0] * naxes[1]; 
  dist = (PRECISION*) malloc(nb_step*sizeof(PRECISION));
  mean = (PRECISION*) malloc(nb_step*sizeof(PRECISION));
  rms = (PRECISION*) malloc(nb_step*sizeof(PRECISION));
  counts = (long*) malloc(nb_step*sizeof(long));
  
  //  if( Debug )
  //printf("Function profile. center = %f %f, step = %f, NBstep = %ld\n",xcenter,ycenter,step,nb_step);
  
  for (i=0;i<nb_step;i++)
    {
      dist[i] = 0.0;
      mean[i] = 0.0;
      rms[i] = 0.0;
      counts[i] = 0;
    }

  if ((fp=fopen(outfile,"w"))==NULL)
    printf("error : can't open file %s\n",outfile);
 
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++){
      distance = sqrt((1.0*ii-xcenter)*(1.0*ii-xcenter)+(1.0*jj-ycenter)*(1.0*jj-ycenter));
      i = (long) distance/step;
      if(i<nb_step)
	{
	  dist[i] += distance;
	  mean[i] += data.image[ID].arrayD[jj*naxes[0]+ii];
	  rms[i] += data.image[ID].arrayD[jj*naxes[0]+ii]*data.image[ID].arrayD[jj*naxes[0]+ii];
	  counts[i] += 1;
	}
    }

  for (i=0;i<nb_step;i++)
    {
      dist[i] /= counts[i];
      mean[i] /= counts[i];
      rms[i] = 0.0;
    }

  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++){
      distance = sqrt((1.0*ii-xcenter)*(1.0*ii-xcenter)+(1.0*jj-ycenter)*(1.0*jj-ycenter));
      i = (long) distance/step;
      if(i<nb_step)
	{
	  rms[i] += (data.image[ID].arrayD[jj*naxes[0]+ii]-mean[i])*(data.image[ID].arrayD[jj*naxes[0]+ii]-mean[i]);
	  //	  counts[i] += 1;
	}
    }

  printf("--------- TEST ----------\n");
  
  for (i=0;i<nb_step;i++)
    {
      //     dist[i] /= counts[i];
      // mean[i] /= counts[i];
      // rms[i] = sqrt(rms[i]-1.0*counts[i]*mean[i]*mean[i])/sqrt(counts[i]);
      rms[i] = sqrt(rms[i]/counts[i]);
      fprintf(fp,"%.18f %.18f %.18f %ld %ld\n",dist[i],mean[i],rms[i],counts[i],i);
    }
  

  fclose(fp);

  free(counts);
  free(dist);
  free(mean);
  free(rms);
  return(0);
}

int profile2im(char *profile_name, long nbpoints, long size, PRECISION xcenter, PRECISION ycenter, PRECISION radius, char *out)
{
  FILE *fp;
  long ID;
  PRECISION *profile_array;
  long i;
  long index;
  PRECISION tmp;
  long ii,jj;
  PRECISION r,x;

  ID = create_2Dimage_ID(out,size,size);
  profile_array = (PRECISION*) malloc(sizeof(PRECISION)*nbpoints);

  if((fp=fopen(profile_name,"r"))==NULL)
    {
      printf("ERROR: cannot open profile file \"%s\"\n",profile_name);
      exit(0);
    }
  for(i=0;i<nbpoints;i++)
    {
      if(fscanf(fp,"%ld %f\n",&index,&tmp)!=2)
	{
	  printf("ERROR: fscanf, %s line %d\n",__FILE__,__LINE__);
	  exit(0);
	}
      profile_array[i] = tmp;
    }
  fclose(fp);

  for(ii=0;ii<size;ii++)
    for(jj=0;jj<size;jj++)
      {
	r = sqrt((1.0*ii-xcenter)*(1.0*ii-xcenter)+(1.0*jj-ycenter)*(1.0*jj-ycenter))/radius;
	i = (long) (r*nbpoints);
	x = r*nbpoints-i; // 0<x<1

	if(i+1<nbpoints)
	  {
	    data.image[ID].arrayD[jj*size+ii] = (1.0-x)*profile_array[i]+x*profile_array[i+1];
	  }
	else
	  if(i<nbpoints)
	    data.image[ID].arrayD[jj*size+ii] = profile_array[i];
      }

  free(profile_array);

  return(0);
}

int printpix(char *ID_name, char *filename)
{
  int ID;
  long ii,jj,kk;
  long nelements;
  long nbaxis;
  long naxes[3];
  FILE *fp;

  long iistep = 1;
  long jjstep = 1;

  ID = variable_ID("_iistep");
  if(ID!=-1)
      {
      iistep = (long) (0.1+data.variable[ID].value);
      printf("iistep = %ld\n", iistep);
    }
  ID = variable_ID("_jjstep");
  if(ID!=-1)
    {
      jjstep = (long) (0.1+data.variable[ID].value);
       printf("jjstep = %ld\n", jjstep);
    }

  if((fp=fopen(filename,"w"))==NULL)
    {
      printf("ERROR: cannot open file \"%s\"\n",filename);
      exit(0);
    }

  ID = image_ID(ID_name);
  nbaxis = data.image[ID].naxis;
  if(nbaxis==2)
    {
      naxes[0] = data.image[ID].size[0];
      naxes[1] = data.image[ID].size[1];    
      nelements = naxes[0] * naxes[1]; 
      for (ii = 0; ii < naxes[0]; ii+=iistep)
	{
	  for (jj = 0; jj < naxes[1]; jj+=jjstep)
	    {
	      //  fprintf(fp,"%f ",data.image[ID].arrayD[jj*naxes[0]+ii]);
	      fprintf(fp,"%ld %ld %g\n",ii,jj,data.image[ID].arrayD[jj*naxes[0]+ii]);
	    }
	  fprintf(fp,"\n");
	}      
    }
  if(nbaxis==3)
    {
      naxes[0] = data.image[ID].size[0];
      naxes[1] = data.image[ID].size[1];    
      naxes[2] = data.image[ID].size[2];    
      nelements = naxes[0] * naxes[1]; 
      for (ii = 0; ii < naxes[0]; ii+=iistep) 
	for (jj = 0; jj < naxes[1]; jj+=jjstep)
	  for (kk = 0; kk < naxes[2]; kk++)
	    {
	      fprintf(fp,"%ld %ld %ld %f\n",ii,jj,kk,data.image[ID].arrayD[kk*naxes[1]*naxes[0]+jj*naxes[0]+ii]);
	    }
    
    }
  fclose(fp);

  return(0);
}

PRECISION background_photon_noise(char *ID_name)
{
  int ID;
  long ii,jj;
  long naxes[2];
  PRECISION value1,value2,value3,value;
  PRECISION *array;
  long nelements;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    
  nelements = naxes[0]*naxes[1];
  
  array = (PRECISION*) malloc(naxes[1]*naxes[0]*sizeof(PRECISION));
  for (jj = 0; jj < naxes[1]; jj++) 
    for (ii = 0; ii < naxes[0]; ii++)
      array[jj*naxes[0]+ii] = data.image[ID].arrayD[jj*naxes[0]+ii];

  quick_sort(array,nelements);
  
  /* uses the repartition function F of the normal distribution law */
  /* F(0) = 0.5 */
  /* F(-0.1 * sig) = 0.460172162723 */
  /* F(-0.2 * sig) = 0.420740290562 */
  /* F(-0.3 * sig) = 0.382088577811 */
  /* F(-0.4 * sig) = 0.34457825839 */
  /* F(-0.5 * sig) = 0.308537538726 */
  /* F(-0.6 * sig) = 0.27425311775 */
  /* F(-0.7 * sig) = 0.241963652223 */
  /* F(-0.8 * sig) = 0.211855398584 */
  /* F(-0.9 * sig) = 0.184060125347 */
  /* F(-1.0 * sig) = 0.158655253931 */
  /* F(-1.1 * sig) = 0.135666060946 */
  /* F(-1.2 * sig) = 0.115069670222 */
  /* F(-1.3 * sig) = 0.0968004845855 */

  /* calculation using F(-0.9*sig) and F(-1.3*sig) */
  value1 = array[(long) (0.184060125347*naxes[1]*naxes[0])]-array[(long) (0.0968004845855*naxes[1]*naxes[0])];
  value1 /= (1.3-0.9);
  printf("(-1.3 -0.9) %f\n",value1);

  /* calculation using F(-0.6*sig) and F(-1.3*sig) */
  value2 = array[(long) (0.27425311775*naxes[1]*naxes[0])]-array[(long) (0.0968004845855*naxes[1]*naxes[0])];
  value2 /= (1.3-0.6);
  printf("(-1.3 -0.6) %f\n",value2);

  /* calculation using F(-0.3*sig) and F(-1.3*sig) */
  value3 = array[(long) (0.382088577811*naxes[1]*naxes[0])]-array[(long) (0.0968004845855*naxes[1]*naxes[0])];
  value3 /= (1.3-0.3);
  printf("(-1.3 -0.3) %f\n",value3);

  value = value3;
  
  free(array);
  return(value);
}

int test_structure_function(char *ID_name, long NBpoints, char *ID_out)
{
  int ID,ID1,ID2;
  long ii1,ii2,jj1,jj2,i,ii,jj;
  long naxes[2];
  long nelements;
  PRECISION v1,v2;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    

  nelements = naxes[0]*naxes[1];

  ID1=create_2Dimage_ID("tmp1",naxes[0],naxes[1]);
  ID2=create_2Dimage_ID("tmp2",naxes[0],naxes[1]);
  
  for(i=0;i<NBpoints;i++)
    {
      ii1=(long) (data.INVRANDMAX*rand()*naxes[0]);
      jj1=(long) (data.INVRANDMAX*rand()*naxes[1]);
      ii2=(long) (data.INVRANDMAX*rand()*naxes[0]);
      jj2=(long) (data.INVRANDMAX*rand()*naxes[1]);
      v1=data.image[ID].arrayD[jj1*naxes[0]+ii1];
      v2=data.image[ID].arrayD[jj2*naxes[0]+ii2];
      ii=(ii1-ii2);
      if(ii<0)
	ii=-ii;
      jj=(jj1-jj2);
      if(jj<0)
	jj=-jj;
      data.image[ID1].arrayD[jj*naxes[0]+ii] += (v1-v2)*(v1-v2);
      data.image[ID2].arrayD[jj*naxes[0]+ii] += 1.0;
    }
  arith_image_div("tmp1","tmp2",ID_out);


  return(0);
}

int full_structure_function(char *ID_name, long NBpoints, char *ID_out)
{
  int ID,ID1,ID2;
  long ii1,ii2,jj1,jj2;
  long naxes[2];
  PRECISION v1,v2;
  long i=0;
  long STEP1=2;
  long STEP2=3;

  ID = image_ID(ID_name);
  naxes[0] = data.image[ID].size[0];
  naxes[1] = data.image[ID].size[1];    

  ID1=create_2Dimage_ID("tmp1",naxes[0],naxes[1]);
  ID2=create_2Dimage_ID("tmp2",naxes[0],naxes[1]);
  

  for(ii1=0;ii1<naxes[0];ii1+=STEP1)
    {
      printf(".");
      for(jj1=0;jj1<naxes[1];jj1+=STEP1)
	{
	  if(i<NBpoints)
	    {
	      i++;
	      fflush(stdout);
	      for(ii2=0;ii2<naxes[0];ii2+=STEP2)
		for(jj2=0;jj2<naxes[1];jj2+=STEP2)
		  if((ii2>ii1)&&(jj2>jj1))
		    {
		      v1=data.image[ID].arrayD[jj1*naxes[0]+ii1];
		      v2=data.image[ID].arrayD[jj2*naxes[0]+ii2];
		      data.image[ID1].arrayD[(jj2-jj1)*naxes[0]+ii2-ii1] += (v1-v2)*(v1-v2);
		      data.image[ID2].arrayD[(jj2-jj1)*naxes[0]+ii2-ii1] += 1.0;
		    }
	    }
	}
    }
  printf("\n");

  arith_image_div("tmp1","tmp2",ID_out);

  return(0);
}

int fft_structure_function(char *ID_in, char *ID_out)
{
  long ID;
  PRECISION value;
  long nelement;

  autocorrelation(ID_in,"stftmp");
  ID=image_ID("stftmp");
  nelement = data.image[ID].nelement;
  value=-data.image[ID].arrayD[0];
  
  arith_image_cstadd("stftmp",value,"stftmp1");
  delete_image_ID("stftmp");
  arith_image_cstmult("stftmp1",-2.0/sqrt(nelement),ID_out);
  delete_image_ID("stftmp1");

  return(0);
}

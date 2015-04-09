#include "cow-cfg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if (COW_RNPL)
#include <bbhutil.h>
#endif // COW_RNPL

#include "read-ff-sdf.h"

#define NUM_VARIABLES 6

int getEBFromSDF(int fileRank,
		 int timestep,
		 char *fileStub,
		 double *EB[NUM_VARIABLES],
		 int *startIndices,
		 int *sizes,
		 double *time)
{
#if (COW_RNPL)
   char* variableNames[NUM_VARIABLES] = {"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
   const static double global_bbox[6] = {-1, 1, -1, 1, -1, 1};
   const static int numFiles = NUM_VARIABLES;
   const static int numTrim = 7;

   int shape[3] = {1, 1, 1};
   int rank;
   char cnames[128]="";
   char fileName[128]="";
   int fileLevel = timestep+1; 
   int nf;
   int size, csize;
   double *data=NULL;
   double *coords=NULL;
   for (nf=0; nf<numFiles; nf++) {
      strcpy(fileName, fileStub);
      char fileNameEnd[50];
      sprintf(fileNameEnd,"%s_em_tl2_%d.sdf",variableNames[nf],fileRank);
      strcat(fileName, fileNameEnd);
      int i,j,k;
      if (nf==0) {
         if (!gft_read_shape(fileName, fileLevel,shape)) {
            gft_close_sdf_stream(fileName);
            return 1;
         }
         gft_read_rank(fileName, fileLevel, &rank);
         if (rank!=2) {
            gft_close_sdf_stream(fileName);
            printf("Error rank!=2\n");
            return 1;
         }

         size = 1;
         for (i=0; i<rank; i++) size *= shape[i];
         data = malloc(size*sizeof(double));
         csize = 0;
         for (i=0; i<rank; i++) csize += shape[i];
         coords = malloc(csize*sizeof(double));
      }
      int shapeNew[3] = {1, 1, 1};
      int ret = gft_read_full(fileName, /* IN: name of this file */
                 fileLevel, /* IN: time level */
                 shapeNew, /* OUT: shape of data set */
                 cnames, /* OUT: names of coordinates */
                 rank, /* IN: rank of data set */
                 time, /* OUT: time value */
                 coords, /* OUT: values of coordinates */
                 data ); /* OUT: actual data */
      if (!ret) return 1;
      for (i=0; i<rank; i++) if (shape[i]!=shapeNew[i]) {
         printf("Error %s and %s shapes do not match\n", variableNames[0],variableNames[nf]);
         return 1;
      }
      int sizeOut = 1;
      int csizeOut = 0; 
      int shapeOut[3];
      int num_trim[6];       
      double dx = coords[1]-coords[0];
      int offset = 0;
      for (i=0; i<rank; i++) { 
         num_trim[2*i] = 0;
         num_trim[2*i+1] = numTrim;
         if (coords[shape[i]-1+offset]>global_bbox[2*i+1]) num_trim[2*i+1] = num_trim[2*i+1]-1; 
         offset += shape[i];
      }
      for (i=0; i<rank; i++) {  
         shapeOut[i] = (shape[i]-num_trim[2*i]-num_trim[2*i+1]);
         sizeOut *= shapeOut[i];
         csizeOut += shapeOut[i];
         if (sizeOut<=0) break;
         sizes[i] = shapeOut[i];
      }
      if (sizeOut>0) { 
         double *dataOut = malloc(sizeOut*sizeof(double));
         double *coordsOut = malloc(csizeOut*sizeof(double));
         int n = 0;
         if (rank==3) for (k=num_trim[4]; k<(shape[2]-num_trim[5]); k++) for (j=num_trim[2]; j<(shape[1]-num_trim[3]); j++) for (i=num_trim[0]; i<(shape[0]-num_trim[1]); i++) {
            dataOut[n++] = data[i+j*shape[0]+k*shape[0]*shape[1]];
         }
         else if (rank==2) for (j=num_trim[2]; j<(shape[1]-num_trim[3]); j++) for (i=num_trim[0]; i<(shape[0]-num_trim[1]); i++) {
            dataOut[n++] = data[i+j*shape[0]];
         }
         n = 0; 
         offset = 0;
         for (i=0; i<rank; i++) { 
            startIndices[i] = (int) ((coords[num_trim[2*i]+offset]-global_bbox[2*i])/dx+0.1);
            for (j=num_trim[2*i]; j<(shape[i]-num_trim[2*i+1]); j++) coordsOut[n++] = coords[j+offset];
            offset += shape[i];
         }
         //gft_write_sdf_stream(outFileName, time, shapeOut, cnames, rank, coordsOut, dataOut);
         EB[nf] = dataOut;
         //free(dataOut);
         free(coordsOut);
      }
      gft_close_sdf_stream(fileName);
   }
   if (data) free(data);
   if (coords) free(coords);

#else
   printf("[read-ff-sdf] error: cow was not compiled with the rnpl library\n");
#endif // COW_RNPL

   return 0;
}

#ifndef READFFSDF_H
#define READFFSDF_H

int getEBFromSDF(int fileRank,
		 int timestep,
		 char *fileStub,
		 double *EB[6],
		 int *startIndices,
		 int *sizes,
		 double *time);

#endif // READFFSDF_H

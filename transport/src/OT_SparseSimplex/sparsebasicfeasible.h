#ifndef SPARSEBASICFEASIBLE_H_
#define SPARSEBASICFEASIBLE_H_

void dedewithchannels(int m, int n, int no_of_ones, int *basis, int *channels_byrow_over, int **channels_byrow);
void findblocks(int m, int n, int *basis, int *nblock, int *rowblock, int *colblock);

#endif


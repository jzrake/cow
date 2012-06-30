
#include <string.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"


cow_histogram *cow_histogram_new()
{
  cow_histogram *h = (cow_histogram*) malloc(sizeof(cow_histogram));
  cow_histogram hist = {
    .nbinsx = 1,
    .nbinsy = 1,
    .x0 = 0.0,
    .x1 = 1.0,
    .y0 = 0.0,
    .y1 = 1.0,
    .bedgesx = NULL,
    .bedgesy = NULL,
    .weight = NULL,
    .counts = NULL,
    .nickname = NULL,
    .fullname = NULL,
    .binmode = COW_BINNING_LINSPACE,
    .n_dims = 0
  } ;
  *h = hist;
  return h;
}

void cow_histogram_commit(cow_histogram *h)
{
  h->n_dims = h->nbinsy > 1 ? 2 : 1;
  if (h->n_dims == 1) {
    const double dx = (h->x1 - h->x0) / h->nbinsx;

    h->bedgesx = (double*) malloc((h->nbinsx+1)*sizeof(double));
    h->weight = (double*) malloc((h->nbinsx)*sizeof(double));
    h->counts = (long*) malloc((h->nbinsx)*sizeof(long));
    for (int n=0; n<h->nbinsx+1; ++n) {
      if (h->binmode == COW_BINNING_LOGSPACE) {
        h->bedgesx[n] = h->x0 * pow(h->x1 / h->x0, (double)n / h->nbinsx);
      }
      else {
        h->bedgesx[n] = h->x0 + n * dx;
      }
    }
    for (int n=0; n<h->nbinsx; ++n) {
      h->counts[n] = 0;
      h->weight[n] = 0.0;
    }
  }
  else if (h->n_dims == 2) {
    const int nbins = h->nbinsx * h->nbinsy;
    const double dx = (y1-y0) / h->nbinsx;
    const double dy = (y1-y0) / h->nbinsy;

    h->bedgesx = (double*) malloc((h->nbinsx+1)*sizeof(double));
    h->bedgesy = (double*) malloc((h->nbinsy+1)*sizeof(double));
    h->weight = (double*) malloc(nbins*sizeof(double));
    h->counts = (long*) malloc(nbins*sizeof(long));
    for (int n=0; n<h->nbinsx+1; ++n) {
      if (h->binmode == COW_BINNING_LOGSPACE) {
        h->bedgesx[n] = h->x0 * pow(h->x1/h->x0, (double)n / h->nbinsx);
      }
      else {
        h->bedgesx[n] = h->x0 + n * dx;
      }
    }
    for (int n=0; n<h->nbinsy+1; ++n) {
      if (h->binmode == COW_BINNING_LOGSPACE) {
        h->bedgesy[n] = h->y0 * pow(h->y1/h->y0, (double)n / h->nbinsy);
      }
      else {
        h->bedgesy[n] = h->y0 + n * dy;
      }
    }
    for (int n=0; n<nbins; ++n) {
      h->counts[n] = 0;
      h->weight[n] = 0.0;
    }
  }
}
void cow_histogram_del(cow_histogram *h)
{
  free(h->bedgesx);
  free(h->bedgesy);
  free(h->weight);
  free(h->counts);
  free(h->nickname);
  free(h->fullname);
  free(h);
}

void cow_histogram_setbinmode(cow_histogram *h, int binmode)
{
  h->binmode = binmode;
}
void cow_histogram_setnbins(cow_histogram *h, int dim, int nbins)
{
  switch (dim) {
  case 0: h->nbinsx = nbins; break;
  case 1: h->nbinsy = nbins; break;
  case COW_ALL_DIMS: h->nbinsx = h->nbinsy = nbins; break;
  default: break;
  }
}
void cow_histogram_setlower(cow_histogram *h, int dim, double v0)
{
  switch (dim) {
  case 0: h->x0 = v0; break;
  case 1: h->y0 = v0; break;
  case COW_ALL_DIMS: h->x0 = h->y0 = v0; break;
  default: break;
  }
}
void cow_histogram_setupper(cow_histogram *h, int dim, double v1)
{
  switch (dim) {
  case 0: h->x1 = v1; break;
  case 1: h->y1 = v1; break;
  case COW_ALL_DIMS: h->x1 = h->y1 = v1; break;
  default: break;
  }
}
void cow_histogram_setfullname(cow_histogram *h, const char *fullname)
{
  h->fullname = (char*) realloc(h->fullname, strlen(fullname)+1);
  strcpy(h->fullname, fullname);
}
void cow_histogram_setnickname(cow_histogram *h, const char *nickname)
{
  h->nickname = (char*) realloc(h->nickname, strlen(nickname)+1);
  strcpy(h->nickname, nickname);
}

void cow_histogram_addsample1(cow_histogram *h, double x, double w)
{
  for (int n=0; n<h->nbinsx; ++n) {
    if (h->bedgesx[n] < x && x < h->bedgesx[n+1]) {
      h->weight[n] += w;
      h->counts[n] += 1;
      return;
    }
  }
}
void cow_histogram_addsample2(cow_histogram *h, double x, double y, double w)
{
  int nx=-1, ny=-1;
  for (int n=0; n<h->nbinsx; ++n) {
    if (h->bedgesx[n] < x && x < h->bedgesx[n+1]) {
      nx = n;
      break;
    }
  }
  for (int n=0; n<h->nbinsy; ++n) {
    if (h->bedgesy[n] < y && y < h->bedgesy[n+1]) {
      ny = n;
      break;
    }
  }
  if (nx == -1 || ny == -1) {
    return;
  }
  else {
    h->counts[nx * h->nbinsy + ny] += 1;
    h->weight[nx * h->nbinsy + ny] += w;
    return;
  }
}
void cow_histogram_dumpascii(cow_histogram *h, const char *fn)
{
  FILE *file = fopen(fn, "w");
  if (file == NULL) return;
  if (h->n_dims == 1) {
    for (int n=0; n<h->nbinsx; ++n) {
      fprintf(file, "%f %f\n", 0.5*(h->bedgesx[n] + h->bedgesx[n+1]),
	      h->weight[n]/h->counts[n]);
    }
  }
  else if (h->n_dims == 2) {
    for (int nx=0; nx<h->nbinsx; ++nx) {
      for (int ny=0; ny<h->nbinsy; ++ny) {
	fprintf(file, "%f %f %f\n",
		0.5*(h->bedgesx[nx] + h->bedgesx[nx+1]),
		0.5*(h->bedgesy[ny] + h->bedgesy[ny+1]),
		h->counts[nx * h->nbinsy + ny] / h->weight[nx * h->nbinsy + ny]);
      }
    }
  }
  fclose(file);
}

/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#ifndef MRST_GEOMETRY_H_INCLUDED
#define MRST_GEOMETRY_H_INCLUDED

void compute_face_geometry(long long ndims, double *coords, long long nfaces,
                           unsigned *nodepos, long long *facenodes,
                           double *fnormals, double *fcentroids,
                           double *fareas);
void compute_cell_geometry(long long ndims, double *coords,
                           unsigned* nodepos, long long *facenodes, long long *neighbours,
                           double *fnormals,
                           double *fcentroids, long long ncells,
                           unsigned* facepos, long long *cellfaces,
                           double *ccentroids, double *cvolumes);

#endif /* MRST_GEOMETRY_H_INCLUDED */

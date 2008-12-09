// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id config-1.0.7-74-gc7a455f Mon Oct 29 23:28:05 2007 -0300$
#ifndef PETSCFEM_MSHPART_H
#define PETSCFEM_MSHPART_H

void metis_part(int nelemfat,Mesh *mesh,
		const int nelemsets,int *vpart,
		int *nelemsetptr,int *n2eptr,
		int *node2elem,int size,const int myrank,
		const int partflag,float *tpwgts,
		int max_partgraph_vertices,
		int iisd_subpart,
		int print_partitioning_statistics);

void match_graph(const dvector<double> &bw,
                 const dvector<double> &bflux,
                 dvector<int> &proc);

void perfo(const dvector<double> &bw,
           const dvector<double> &bflux,
           const dvector<int> &proc, double &perfo_max,
           double &perfo_sum);

#endif



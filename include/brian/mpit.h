#ifndef MPIT_H
#define MPIT_H

/*
 *
 * mpit.h: templated MPI interface
 * BRIAN Software Package Version 3.0
 *
 * $Id: mpit.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (27/11/11): for BRIAN2.2 by FK
 * 0.40 (16/12/13): documented
 * 0.50 (17/11/14): released for BRIAN2.7
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions for the message passing interface (MPI).
*/

#include <mpi.h>

extern bool mpiSend(char* x, unsigned int len, unsigned int proc, MPI_Request* req);
extern bool mpiSend(int *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiSend(float *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiSend(double *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiSend(unsigned int *x, unsigned int len, unsigned int proc, MPI_Request *req);

extern bool mpiRecv(int *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiRecv(float *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiRecv(double *x, unsigned int len, unsigned int proc, MPI_Request *req);
extern bool mpiRecv(unsigned int *x, unsigned int len, unsigned int proc, MPI_Request *req);

extern bool mpiSendBuf(void* x, unsigned int len, unsigned int proc);
extern void *mpiRecvBuf(unsigned int& len, unsigned int proc);

extern int mpiSum(int x, unsigned int len);
extern float mpiSum(float x, unsigned int len);
extern double mpiSum(double x, unsigned int len);
extern unsigned int mpiSum(unsigned int x, unsigned int len);

extern bool mpiGather(int *x, unsigned int n, int *y, int *len, int *off);
extern bool mpiGather(float *x, unsigned int n, float *y, int *len, int *off);
extern bool mpiGather(double *x, unsigned int n, double *y, int *len, int *off);
extern bool mpiGather(unsigned int *x, unsigned int n, unsigned int *y, int *len, int *off);

extern bool mpiScatter(int *x, int *len, int *off, int *y, unsigned int n);
extern bool mpiScatter(float *x, int *len, int *off, float *y, unsigned int n);
extern bool mpiScatter(double *x, int *len, int *off, double *y, unsigned int n);
extern bool mpiScatter(unsigned int *x, int *len, int *off, unsigned int *y, unsigned int n);

extern bool mpiAllgather(int* x, unsigned int n, int* y, int* len, int* off);
extern bool mpiAllgather(float* x, unsigned int n, float* y, int* len, int* off);
extern bool mpiAllgather(double* x, unsigned int n, double* y, int* len, int* off);
extern bool mpiAllgather(unsigned int* x, unsigned int n, unsigned int* y, int* len, int* off);

extern bool mpiWait(unsigned int n, MPI_Request *req, MPI_Status *sta = MPI_STATUSES_IGNORE);
extern bool mpiSync(int *x, unsigned int len, int *y);
extern bool mpiSync(unsigned int *x, unsigned int len, unsigned int *y);
extern unsigned int mpiSize();
extern unsigned int mpiRank();
extern bool mpiInit(int argc, char *argv[]);
extern void mpiFinalize();

#endif

#ifndef INC_CPPTRAJSTDIO_H
#define INC_CPPTRAJSTDIO_H
/*! \file CpptrajStdio.h
    \brief Interface between Cpptraj and Stdio.

    This is a useful abstraction that allows cpptraj to use several 
    often-used functions from the CSTDLIB stdio library without having
    to worry about details. For example, the mprintf function ensures
    that during parallel runs messages are only printed to the master
    thread, etc.
 */
#ifdef MPI
void SetStdioRank(int);
#endif
void mprintf(const char *, ...);
void mprinterr(const char *, ...);
void mflush();
void rprintf(const char *, ...);
void rprinterr(const char *, ...);
void SetWorldSilent(bool);
void SuppressAllOutput();
void SuppressErrorMsg(bool);
void FinalizeIO();
int OutputToFile(const char*);
int ErrToFile(const char*);
void* CpptrajStdout();
#endif

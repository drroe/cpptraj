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
/// Set parallel world rank; will initially supress output for non-master ranks.
void SetStdioRank(int);
#endif
/// Print to STDOUT
void mprintf(const char *, ...);
/// Print to STDERR
void mprinterr(const char *, ...);
/// Flush STDOUT
void mflush();
/// Print to STDOUT; prepend rank in MPI
void rprintf(const char *, ...);
/// Print to STDERR; prepend rank in MPI
void rprinterr(const char *, ...);
/// Can suppress/enable output.
void SetWorldSilent(bool);
/// Suppress all output/error forever. Use by pytraj.
void SuppressAllOutput();
/// Can suppress/enable error output.
void SuppressErrorMsg(bool);
/// Finalize any IO
void FinalizeIO();
/// Redirect output to given file.
int OutputToFile(const char*);
/// Redirect error output to given file.
int ErrToFile(const char*);
/// \retrun Current output FILE* stream as a void pointer.
void* CpptrajStdout();
/// \return Current error FILE* stream as a void pointer.
void* CpptrajStderr();
#endif

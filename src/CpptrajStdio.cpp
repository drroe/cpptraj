#include <cstdio>
#include <cstdarg>

/// If this is false, do not allow any output to be written. 
static bool enable_output_ = true;
/// If true, enable the Xprintf functions.
static bool enable_stdout_ = true;
/// If true, enable the Xprinterr functions.
static bool enable_stderr_ = true;
/// Where normal output (Xprintf) should be written.
static FILE* STDOUT_ = stdout;
/// Where error output (Xprinterr) should be written.
static FILE* STDERR_ = stderr;
#ifdef MPI
/// World rank in parallel, for rprintX functions. -1 means do not print rank prefix.
static int world_rank_ = -1;

/** Set the parallel world rank. Initially, only enable output for the master rank. */
void SetStdioRank(int rankIn) {
  world_rank_ = rankIn;
  if (world_rank_ > 0) {
    enable_output_ = false;
    enable_stdout_ = false;
    enable_stderr_ = false;
  }
  //printf("[%i] DEBUG: stdout= %i  stderr= %i\n", world_rank_, (int)enable_stdout_, (int)enable_stderr_);
}
#endif

//#ifdef PARALLEL_DEBUG_VERBOSE
// -----------------------------------------------------------------------------
/** Master prints message to STDOUT, others to mpidebugfile. */
/*void mprintf(const char*format, ...) {
  va_list args;
  va_start(args,format);
  if (Parallel::World().Master()) {
    vfprintf(STDOUT_, format, args);
    vfprintf(Parallel::mpidebugfile_, format, args);
  } else
    vfprintf(Parallel::mpidebugfile_, format, args);
  va_end(args);
}*/

/** Master prints message to STDERR, others to mpidebugfile. */
/*void mprinterr(const char *format, ...) {
  va_list args;
  va_start(args,format);
  if (Parallel::World().Master()) {
    vfprintf(stderr,format,args);
    vfprintf(Parallel::mpidebugfile_, format, args);
  } else
    vfprintf(Parallel::mpidebugfile_, format, args);
  va_end(args);
}*/
// -----------------------------------------------------------------------------
//#else /* PARALLEL_DEBUG_VERBOSE */

/** Print message to STDOUT. */
void mprintf(const char *format, ...) {
  if (enable_stdout_) {
    va_list args;
    va_start(args, format);
    vfprintf(STDOUT_, format,args);
    va_end(args);
  }
}

/** Print message to STDERR. */
void mprinterr(const char *format, ...) {
  if (enable_stderr_) {
    va_list args;
    va_start(args,format);
    vfprintf(STDERR_, format,args);
    va_end(args);
  }
}
//#endif

// mflush()
/** Call flush on STDOUT. */
void mflush() {
  if (enable_stdout_)
    fflush(STDOUT_);
}

// rprintf()
/** Print message to STDOUT for this worldrank */
void rprintf(const char *format, ...) {
  if (enable_stdout_) {
    va_list args;
    va_start(args, format);
#   ifdef MPI
    char buffer[1024];
    int nc = sprintf(buffer, "[%i]\t", world_rank_);
    nc += vsprintf(buffer + nc, format, args);
    fwrite(buffer, 1, nc, STDOUT_);
#   else
    vfprintf(STDOUT_, format,args);
#   endif
    va_end(args);
  }
}

// rprinterr()
/** Print message to STDERR for this worldrank */
void rprinterr(const char *format, ...) {
  va_list args;
  if (enable_stderr_) {
    va_start(args,format);
#   ifdef MPI
    char buffer[1024];
    int nc = sprintf(buffer, "[%i]\t", world_rank_);
    nc += vsprintf(buffer + nc, format, args);
    fwrite(buffer, 1, nc, STDERR_);
#   else
    vfprintf(STDERR_, format,args);
#   endif
    va_end(args);
  }
}

/** Change status of STDOUT output as long as SuppressAllOutput has not
  * been called.
  * \param silentIn if true, silence STDOUT output, otherwise enable.
  */
void SetWorldSilent(bool silentIn) {
  //printf("DEBUG: Calling SetWorldSilent %i\n", (int)silentIn);
  if (enable_output_) {
    enable_stdout_ = !silentIn;
  }
}

/** Suppress all STDOUT/STDERR output for the entire run. */
void SuppressAllOutput() { 
  //printf("DEBUG: Calling SuppressAllOutput\n");
  enable_output_ = false;
  enable_stdout_ = false;
  enable_stderr_ = false;
}

/** \param supressIn if true, silence STDERR output, otherwise enable. */
void SuppressErrorMsg(bool suppressIn) {
  //printf("DEBUG: Calling SuppressErrorMsg %i\n", (int)suppressIn);
  if (enable_output_) {
    enable_stderr_ = !suppressIn;
  }
}

/** Close output if not stdout/stderr */
void FinalizeIO() {
  if (STDOUT_ != stdout) {
    fclose(STDOUT_);
    STDOUT_ = stdout;
  }
  if (STDERR_ != stderr) {
    fclose(STDERR_);
    STDERR_ = stderr;
  }
}

/** Redirect output to file. If no name given assume STDOUT. */
int OutputToFile(const char* fname) {
  if (!enable_stdout_) return 0;
  if (STDOUT_ != stdout) {
    fclose(STDOUT_);
    STDOUT_ = stdout;
  }
  if (fname != 0) {
    mprintf("Info: Redirecting output to file '%s'\n", fname);
    STDOUT_ = fopen(fname, "wb");
    if (STDOUT_ == 0) {
      mprinterr("Error: Could not open output file '%s'\n", fname);
      return 1;
    }
  }
  return 0;
}

/** \return STDOUT_ cast to a void pointer. */
void* CpptrajStdout() {
  return (void*)STDOUT_;
}

/** Redirect errors to file. If no name given assume STDERR. */
int ErrToFile(const char* fname) {
  if (!enable_stderr_) return 0;
  if (STDERR_ != stderr) {
    fclose(STDERR_);
    STDERR_ = stderr;
  }
  if (fname != 0) {
    mprintf("Info: Redirecting error to file '%s'\n", fname);
    STDERR_ = fopen(fname, "wb");
    if (STDERR_ != 0) {
      mprinterr("Error: Could not open error file '%s'\n", fname);
      return 1;
    }
  }
  return 0;
}

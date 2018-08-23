#include "ProgressTimerBar.h"
#include "CpptrajStdio.h"

const int ProgressTimerBar::UNKNOWN_FRAMES_ = -1;

const char ProgressTimerBar::PROGCHAR_ = '+';

/** Set up progress bar for given number of iterations. */
void ProgressTimerBar::SetupProgress(int nIterationsIn, double intervalIn, double pctIntIn)
{
  total_it_ = nIterationsIn;
  tInterval_ = intervalIn;
  tgtTime_ = intervalIn;
  if (total_it_ > 0) {
    // Total number of iterations is known.
    pctInterval_ = pctIntIn;
    if (pctInterval_ <= 0.0)
      pctInterval_ = 10.0; // default
    else if (pctInterval_ > 100.0)
      pctInterval_ = 100.0; // TODO should be set to default?
    // Adjust percent interval if necessary.
    if ((double)total_it_ < pctInterval_)
      pctInterval_ = 100.0 / (double)total_it_;
    tgtPct_ = pctInterval_;
    // Determine next iteration that should be printed.
    tgt_it_ = (int)((double)total_it_ * (pctInterval_/100.0)) - 1;
    mprintf("%2.0f%%", 0.0);
  } else {
    // Total number of iterations unknown.
    total_it_ = UNKNOWN_FRAMES_;
    // If pctIntIn is negative, assume this is the number of iterations
    // between printing a character. Otherwise use a default.
    if (pctIntIn < 0.0)
      tgt_it_ = (int)-pctIntIn;
    else
      tgt_it_ = 200; // default;
    iInterval_ = tgt_it_;
    mprintf("\tProgress: '%c' = %i iterations.\n", PROGCHAR_, iInterval_);
    tgt_it_ = tgt_it_ - 1;
  }
  mflush();
  if (tInterval_ > 0.0) time_.Start();
}

/** Check if progress bar or elapsed time should be printed. This is intended
  * to be called at the beginning of a loop, i.e. before work has been done.
  */
void ProgressTimerBar::Update(int iteration) {
  if (total_it_ == UNKNOWN_FRAMES_) {
    if (iteration >= tgt_it_) {
      mprintf("%c", PROGCHAR_);
      //mprintf(" (%i) %c\n", iteration+1, PROGCHAR_); // DEBUG
      mflush();
      tgt_it_ += iInterval_;
      // TODO number of columns?
    }
    if (tInterval_ > 0.0) {
      double elapsed = time_.Elapsed();
      if (elapsed > tgtTime_) {
        mprintf("\n\t%i iterations in %.2f s (%.2f /s).\n",
                iteration+1, elapsed, (double)(iteration+1)/elapsed);
        mflush();
        tgtTime_ += tInterval_;
      }
    }
  } else {
    if (iteration > tgt_it_) {
      double pctComplete = ((double)iteration / (double)total_it_) * 100.0;
      mprintf(" %2.0f%%", pctComplete);
      //mprintf(" (%i) %2.0f%%\n", iteration+1, pctComplete); // DEBUG
      mflush();
      // Determine next percent.
      tgtPct_ += pctInterval_;
      tgt_it_ = (int)((tgtPct_/100.0) * (double)total_it_) - 1;
    }
    if (tInterval_ > 0.0) {
      double elapsed = time_.Elapsed();
      if (elapsed > tgtTime_) {
        tgtTime_ += tInterval_;
        int it_remaining = total_it_ - iteration;
        double it_per_s = ((double)(iteration+1)) / elapsed;
        double time_remaining = ((double)it_remaining) / it_per_s;
        mprintf("\n\t%i iterations in %.2f s, %.2f s remaining.\n",
                iteration+1, elapsed, time_remaining);
        mflush();
      }
    }
  }
}

/** Finish up progress bar. */
void ProgressTimerBar::Finish() const {
  if (total_it_ == UNKNOWN_FRAMES_) {
    mprintf("\n");
  } else {
    mprintf(" %2.0f%% Complete.\n", 100.0);
    //mprintf("%.2f s\n", time_.Elapsed()); // DEBUG
  }
}

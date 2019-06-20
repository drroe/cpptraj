#ifndef INC_ACCUMULATOR_H
#define INC_ACCUMULATOR_H
/** Abstract base class for accumulator using stable version of Welford's
  * online algorithm.
  */
class Accumulator
{
  public:
    Accumulator() : n_(0) {}
    virtual ~Accumulator() {}
    virtual void accumulate(double) = 0;
    virtual double mean() const = 0;
    virtual double variance() const = 0;
    virtual double popVariance() const = 0;
    virtual double StdDev() const = 0;
    unsigned int nData() const { return (unsigned int)n_; }
  protected:
    static inline void onlineAccumulate(double xIn, double nVals, double& mean, double& M2) {
      double delta = xIn - mean;
      mean += delta / nVals;
      M2 += delta * (xIn - mean);
    }

    double n_;
};

/** Accumulate scalar non-periodic values (degrees). */
class Accumulator_Scalar : public Accumulator
{
  public:
    Accumulator_Scalar() : mean_(0), M2_(0) {}

    void accumulate(double xIn) {
      n_++;
      onlineAccumulate(xIn, n_, mean_, M2_);
    }

    double mean() const { return mean_; }

    double variance() const {
      if (n_ < 2) return 0;
      return M2_ / (n_ - 1);
    }

    double popVariance() const {
      if (n_ < 1) return 0;
      return M2_ / n_;
    }

    double StdDev() const;
  private:
    double mean_;
    double M2_;
};

/** Accumulate scalar periodic torsion values. */
class Accumulator_Torsion : public Accumulator
{
  public:
    Accumulator_Torsion() : meanSp_(0), meanCp_(0), MS2_(0), MC2_(0) {}

    void accumulate(double);

    double mean() const;

    double variance() const { return popVariance(); }

    double popVariance() const;

    double StdDev() const;

  private:
    double meanSp_;
    double meanCp_;
    double MS2_;
    double MC2_;
};
#endif

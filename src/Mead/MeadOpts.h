#ifndef INC_MEAD_MEADOPTS_H
#define INC_MEAD_MEADOPTS_H
#ifdef HAS_MEAD
namespace Cpptraj {
namespace Mead {
/// Hold MEAD options
class MeadOpts {
  public:
    MeadOpts();
    /// Set MEAD verbosity level
    static void MeadVerbosity(int);

    void SetEpsIn(double e)       { epsin_ = e; }
    void SetEpsExt(double e)      { epsext_ = e; }
    void SetEpsVac(double e)      { epsvac_ = e; }
    void SetSolRad(double s)      { solrad_ = s; }
    void SetSterLn(double s)      { sterln_ = s; }
    void SetIonicStr(double i)    { ionicstr_ = i; }
    void SetTemperature(double t) { temp_ = t; }

    double EpsIn()       const { return epsin_; }
    double EpsExt()      const { return epsext_; }
    double EpsVac()      const { return epsvac_; }
    double SolRad()      const { return solrad_; }
    double SterLn()      const { return sterln_; }
    double IonicStr()    const { return ionicstr_; }
    double Temperature() const { return temp_; }

  private:
    double epsin_;    ///< Dielectric constant of molecular interior
    double epsext_;   ///< Dielectric constant  of solvent
    double epsvac_;   ///< Dielectric constant of vacuum.
    double solrad_;   ///< Solvent probe radius to determine boundary between epsin and epsext.
    double sterln_;   ///< Size of ion exclusion layer thickness added to atomic radii
    double ionicstr_; ///< Ionic strength (mol/L)
    double temp_;     ///< Temperature (K)
};
}
}
#endif
#endif

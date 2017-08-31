#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include "File.h"
#include "CoordinateInfo.h"
/// The base interface to NetCDF trajectory files.
class NetcdfFile : private File::Base {
  public:
    /// For determining NetCDF trajectory file type
    enum NCTYPE { NC_UNKNOWN = 0, NC_AMBERTRAJ, NC_AMBERRESTART, NC_AMBERENSEMBLE };
    /// \return NetCDF trajectory type of given file.
    static NCTYPE GetNetcdfConventions(File::Name const&);
#   ifndef BINTRAJ
    NetcdfFile() { }
#   else 
    NetcdfFile();
    // Set up NetCDF file for reading
    int NC_setupRead(File::Name const&);
    // Open NetCDF file
    int NC_open();
    // Members of Base that should be public
    using Base::Filename;
    using Base::Close;
    using Base::SetDebug;
  private:
    // ----- Inherited classes -------------------
    int InternalSetup();
    int InternalOpen();
    void InternalClose();
    // -------------------------------------------
    static NCTYPE GetNetcdfConventions(int);
    /// Check conventions version
    void CheckConventionsVersion();
    /// Create new file
    int CreateNewFile();
    /// Setup existing file
    int SetupExistingFile();
    /// DEBUG - Write start and count arrays to STDOUT
    void WriteIndices() const;
    /// DEBUG - Write all variable IDs to STDOUT
    void WriteVIDs() const;
    /// Convert given float array to double.
    inline void FloatToDouble(double*,const float*) const;
    /// Convert given double array to float.
    inline void DoubleToFloat(float*,const double*) const; 

    size_t start_[4];     ///< Array starting indices
    size_t count_[4];     ///< Array counts

    int ncid_;            ///< NetCDF file ID
    int ncatom_;          ///< Number of atoms
    int ncatom3_;         ///< Number of coordinates (# atoms * 3)
    int ncframe_;         ///< Total number of frames in file
    int remd_dimension_;  ///< Number of replica dimensions
    NCTYPE type_;         ///< NetCDF trajectory type
    // Variable IDs
    int TempVID_;         ///< Temperature variable ID.
    int coordVID_;        ///< Coordinates variable ID.
    int velocityVID_;     ///< Velocity variable ID.
    int frcVID_;          ///< Force variable ID.
    int spatialVID_;      ///< Spatial (x, y, z) variable ID
    int cell_spatialVID_; ///< Box length labels variable ID
    int cell_angularVID_; ///< Box angle labels variable ID
    int cellAngleVID_;    ///< Box angles variable ID.
    int cellLengthVID_;   ///< Box lengths variable ID.
    int timeVID_;         ///< Time variable ID.
    int indicesVID_;      ///< Variable ID for replica indices.
    int repidxVID_;       ///< Variable ID for overall replica index.
    int crdidxVID_;       ///< Variable ID for overall coordinate index.
    // Dimension IDs
    int ensembleDID_;     ///< Ensemble dimenison ID
    int frameDID_;        ///< Frames dimension ID
    int atomDID_;         ///< Atom dimension ID
    int spatialDID_;      ///< Spatial dimension ID (3)
    int labelDID_;        ///< Box angle labels dimension ID (alpha, beta, gamma)
    int cell_spatialDID_; ///< Box lengths dimension ID
    int cell_angularDID_; ///< Box angles dimension ID
#   endif /* BINTRAJ */
};
#ifdef BINTRAJ
// ----- Inline Functions ------------------------------------------------------
/** Convert float coords to double coords
  * NOTE: natom3 needs to match up with size of Coord!
  */
void NetcdfFile::FloatToDouble(double* X, const float* Coord) const {
  for (int i=0; i < ncatom3_; ++i)
    X[i] = (double)Coord[i];
}

/** Convert double coords to float coords */
void NetcdfFile::DoubleToFloat(float* Coord, const double* X) const {
  for (int i=0; i < ncatom3_; ++i)
    Coord[i] = (float)X[i];
}
#endif /* BINTRAJ */
#endif

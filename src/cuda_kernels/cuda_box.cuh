#ifndef INC_CUDABOX_CUH
#define INC_CUDABOX_CUH
#include "../ImageOption.h"
/// Used to hold box parameters in device memory
template <class T> class cuda_box {
  public:
    /// CONSTRUCTOR
    cuda_box() : boxDev_(0), ucellDev_(0), fracDev_(0), last_err_(cudaSuccess), err_desc_(0) {}
    /// DESTRUCTOR
    ~cuda_box();
    /// Allocate and copy box parameters (lengths, ucell, frac) to device
    int Setup(ImageOption::Type, T const*, T const*, T const*);
    /// \return Box lengths device address
    T const* BoxDev() const { return boxDev_; }
    /// \return Box unit cell vectors device address
    T const* UcellDev() const { return ucellDev_; }
    /// \return Box fractional cell vectors device address
    T const* FracDev() const { return fracDev_; }
    /// \return Last error code encountered
    cudaError_t LastErr() const { return last_err_; }
    /// \return Description of last error code encountered
    const char* LastErrDesc() const { return err_desc_; }
  private:
    /// \return 1 on error, 0 on success
    int check(cudaError_t err, const char* desc) {
      if (err != cudaSuccess) {
        last_err_ = err;
        err_desc_ = desc;
        return 1;
      }
      return 0;
    }

    T* boxDev_;   ///< Address on device for orthogonal box lengths (X Y Z)
    T* ucellDev_; ///< Address on device for box unit cell vectors (3x3 row major matrix)
    T* fracDev_;  ///< Address on device for box fractioncal cell vectors (3x3 row major matrix)
    cudaError_t last_err_; ///< Will be set to last error code if CUDA error occurs
    const char* err_desc_; ///< Will be set to last error description if CUDA error occurs.
};

/** Allocate and copy box parameters. */
template <class T> int cuda_box<T>::Setup(ImageOption::Type imageType,
                                          T const* box, T const* ucell, T const* frac)
{
  if (imageType == ImageOption::ORTHO) {
    if (check(cudaMalloc(((void**)(&boxDev_)), 3 * sizeof(T)), "Allocating box")) return 1;
    if (check(cudaMemcpy(boxDev_,box, 3 * sizeof(T), cudaMemcpyHostToDevice), "Copying box")) return 1;
  } else if (imageType == ImageOption::NONORTHO) {
    if (check(cudaMalloc(((void**)(&ucellDev_)), 9 * sizeof(T)), "Allocating ucell")) return 1;
    if (check(cudaMalloc(((void**)(&fracDev_)), 9 * sizeof(T)), "Allocating frac")) return 1;
    if (check(cudaMemcpy(ucellDev_,ucell, 9 * sizeof(T), cudaMemcpyHostToDevice), "Copying ucell")) return 1;
    if (check(cudaMemcpy(fracDev_,frac, 9 * sizeof(T), cudaMemcpyHostToDevice), "Copying frac")) return 1;
  }
  return 0;
}

/** Destructor */
template <class T> cuda_box<T>::~cuda_box() {
  if (boxDev_ != 0) cudaFree(boxDev_);
  if (ucellDev_ != 0) cudaFree(ucellDev_);
  if (fracDev_ != 0) cudaFree(fracDev_);
}
#endif

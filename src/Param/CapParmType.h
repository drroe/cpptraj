#ifndef INC_PARAM_CAPPARMTYPE_H
#define INC_PARAM_CAPPARMTYPE_H
namespace Cpptraj {
namespace Param {
/// Hold parameters for a water 'cap'
class CapParmType {
  public:
    CapParmType() : natcap_(0), cutcap_(0), xcap_(0), ycap_(0), zcap_(0) {}
    bool HasWaterCap() const { return cutcap_ > 0.0; }
    int NatCap()       const { return natcap_; }
    double CutCap()    const { return cutcap_; }
    double xCap()      const { return xcap_;   }
    double yCap()      const { return ycap_;   }
    double zCap()      const { return zcap_;   }
    void Clear() { natcap_ = 0; cutcap_ = 0.0; xcap_ = 0.0; ycap_ = 0.0; zcap_ = 0.0; }
    void SetNatcap(int n)    { natcap_ = n; }
    void SetCutCap(double c) { cutcap_ = c; }
    void SetXcap(double x)   { xcap_ = x;   }
    void SetYcap(double y)   { ycap_ = y;   }
    void SetZcap(double z)   { zcap_ = z;   }
  private:
    int natcap_;    ///< last atom before the start of the cap of waters
    double cutcap_; ///< the distance from the center of the cap to the outside
    double xcap_;   ///< X coordinate for the center of the cap
    double ycap_;   ///< Y coordinate for the center of the cap
    double zcap_;   ///< Z coordinate for the center of the cap
};
}
}
#endif

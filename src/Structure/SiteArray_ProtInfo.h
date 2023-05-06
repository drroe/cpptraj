#ifndef INC_STRUCTURE_PROTONATOR_PROTINFO_H
#define INC_STRUCTURE_PROTONATOR_PROTINFO_H
/** Class used to save residue name info for protonation/deprotonation. */
class Protonator::ProtInfo {
  public:
    enum PstateType { PROTONATED = 0, DEPROTONATED, UNKNOWN };

    ProtInfo() : state_(UNKNOWN) {}

  private:
    std::string names_[2]; ///< 0 is protonated name, 1 is deprotonated name
    PstateType state_;     ///< Which state is the default
};
#endif

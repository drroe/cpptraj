#ifndef INC_EXEC_MUTATE_H
#define INC_EXEC_MUTATE_H
#include "Exec.h"
/// <Enter description of Exec_Mutate here>
class Exec_Mutate : public Exec {
  public:
    Exec_Mutate() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Mutate(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int debug_;
};
#endif

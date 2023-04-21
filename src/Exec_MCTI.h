#ifndef INC_EXEC_MCTI_H
#define INC_EXEC_MCTI_H
#include "Exec.h"
/// <Enter description of Exec_MCTI here>
class Exec_MCTI : public Exec {
  public:
    Exec_MCTI() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_MCTI(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif

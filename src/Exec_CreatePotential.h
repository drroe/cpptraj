#ifndef INC_EXEC_CREATEPOTENTIAL_H
#define INC_EXEC_CREATEPOTENTIAL_H
#include "Exec.h"
/// Create a potential function 
class Exec_CreatePotential : public Exec {
  public:
    Exec_CreatePotential() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CreatePotential(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif

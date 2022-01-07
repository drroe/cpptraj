#ifndef INC_EXEC_REPLICATE_H
#define INC_EXEC_REPLICATE_H
#include "Exec.h"
/// <Enter description of Exec_Replicate here>
class Exec_Replicate : public Exec {
  public:
    Exec_Replicate() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Replicate(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif

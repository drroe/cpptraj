#ifndef INC_EXEC_OUTPUT_H
#define INC_EXEC_OUTPUT_H
#include "Exec.h"
/// Can be used to change where cpptraj writes output/error messages. 
class Exec_Output : public Exec {
  public:
    Exec_Output() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Output(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif

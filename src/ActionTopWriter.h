#ifndef INC_ACTIONTOPWRITER_H
#define INC_ACTIONTOPWRITER_H
#include <string>
#include <map>
// Forward declares
class ArgList;
class Topology;
class DataSetList;
class DataSet_Topology;
/// Class to hold common functionality for actions that will write modified topologies.
class ActionTopWriter {
  public:
    ActionTopWriter();

    static const char* Keywords();
    static const char* Options();

    /// Parse arguments.
    int InitTopWriter(ArgList&, const char*, int, DataSetList*);

    // FIXME Remove this after debug
    int InitTopWriter(ArgList& a, const char* c, int d) {
      return InitTopWriter(a, c, d, 0);
    }

    /// Write options to stdout.
    void PrintOptions() const;
    /// \return DataSet_Topology for containing modified version of given Topology
    DataSet_Topology* CreateTopSet(Topology const&);
    /// Write the Topology to file(s)
    int WriteTops(Topology const&) const;
  private:
    std::string prefix_;      ///< Prefix for writing topology as <prefix>.<originalname>
    std::string parmoutName_; ///< Output topology file name
    std::string parmOpts_;    ///< Topology file write args
    std::string typeStr_;     ///< Label for kind of topology being written.
    DataSetList* masterDSL_;  ///< Pointer to the master DataSetList
    int debug_;               ///< Debug level to pass to Topology file writer.
    typedef std::map<const Topology*,DataSet_Topology*> Imap;
    typedef std::pair<const Topology*,DataSet_Topology*> Ipair;
    Imap modifiedTopMap_;    ///< Map original parm addresses to modified Topologies.
};
#endif

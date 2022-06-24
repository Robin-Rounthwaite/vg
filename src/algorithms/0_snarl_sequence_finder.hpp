#include <gbwtgraph/gbwtgraph.h>
#include "../handle.hpp"
#include "../subgraph.hpp"


namespace vg {
namespace algorithms {

class SnarlSequenceFinder {
  public:
    virtual ~SnarlSequenceFinder()
    {
      cerr << "SnarlSequenceFinder destructor called." << endl;
      debug_check_destructor = true;
    };

    SnarlSequenceFinder(const PathHandleGraph & graph, const SubHandleGraph &snarl,
                   const gbwtgraph::GBWTGraph &gbwt_graph, const id_t &source_id, 
                   const id_t &sink_id, const bool &backwards);
    // {
    //   cerr << "SnarlSequenceFinder constructor called." << endl;
    // };

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
    find_gbwt_haps();

    pair<unordered_set<string>, unordered_set<handle_t>> find_exhaustive_paths();
    
    vector<pair<step_handle_t, step_handle_t>> find_embedded_paths();

  protected:
    // member variables:
    // the handle graph containing the snarl
    const PathHandleGraph &_graph;
    // a subhandlegraph with only the nodes in the snarl
    const SubHandleGraph &_snarl;
    // GBWT graph with snarls to normalize, includes the embedded threads needed for the
    // GBWTPathFinder approach.
    const gbwtgraph::GBWTGraph &_gbwt_graph;
    const id_t &_source_id; 
    const id_t &_sink_id; 
    const bool &_backwards;
    //todo: either make this an optional argument for the class (and pull it from the 
    //todo:    SnarlNormalizer args), or remove as a debug var.
    bool _full_log_print = false; // for printing info that isn't necessarily something gone wrong.
    bool debug_check_destructor = false; //only true after destructed.

    vector<vector<handle_t>> 
    find_haplotypes_not_at_source();
};
}
}
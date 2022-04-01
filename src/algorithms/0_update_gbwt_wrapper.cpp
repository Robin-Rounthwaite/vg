// #include "0_oo_normalize_snarls.hpp"
// #include "0_snarl_sequence_finder.hpp"
// #include <string>

// // #include <deps/seqan/include/seqan/align.h>
// // #include <deps/seqan/include/seqan/graph_align.h>
// // #include <deps/seqan/include/seqan/graph_msa.h>
// #include <seqan/align.h>
// #include <seqan/graph_align.h>
// #include <seqan/graph_msa.h>

// #include <gbwtgraph/gbwtgraph.h>
// #include "../gbwt_helper.hpp"

// #include "../gbwt_helper.hpp"
// #include "../handle.hpp"
// #include "../msa_converter.hpp"
// #include "../snarls.hpp"
// #include "../vg.hpp"

// #include "../types.hpp"
// #include "extract_containing_graph.hpp"

// #include "multipath_mapper.hpp"



#include "0_update_gbwt_wrapper.hpp"
#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"
// #include "../handle.hpp"
// #include "../types.hpp"


////////////////////////////////////////////////////////////
/// Use the gbwt_changelog to perform paralellized rebuild_gbwt():
////////////////////////////////////////////////////////////
namespace vg {
namespace algorithms {


gbwt::GBWT apply_gbwt_changelog(const gbwtgraph::GBWTGraph &gbwt_graph, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog, const gbwt::GBWT& gbwt, const int& threads, const bool& debug_print)
{
    // is the changelog in _graph ids? Or gbwt ids? If gbwt ids, use the gbwt weakly_connected_components. If it's the graph, I need to compute the connected components before changing it.
    // Uh, but actually the gbwt algorithm is gonna be expecting gbwt ids. So basically, use the gbwt info.
    
    // vector<unordered_set<nid_t>> components = gbwtgraph::weakly_connected_components(&_gbwt_graph);
    vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&gbwt_graph);
    // vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&_graph); // this was using the handlegraph algorithm, but I think we want the gbwt's view of the connected components.

    std::unordered_map<nid_t, size_t> node_to_job = get_node_to_job(components);

    std::vector<RebuildJob> jobs = divide_changelog_into_jobs(node_to_job, components, gbwt_changelog); 

    RebuildParameters rebuild_parameters = set_update_gbwt_parameters(threads, debug_print);
    
    gbwt::GBWT output_gbwt = rebuild_gbwt(gbwt, jobs, node_to_job, rebuild_parameters);
    //todo: remove temporary non-parallelized rebuild_gbwt.
    // gbwt::GBWT output_gbwt = rebuild_gbwt(_gbwt, _gbwt_changelog);
    return output_gbwt;
}

std::unordered_map<nid_t, size_t> get_node_to_job(const vector<unordered_set<nid_t>>& weakly_connected_components)
{
    // cerr << "weakly_connected_components.size()" << weakly_connected_components.size() << endl;
    std::unordered_map<nid_t, size_t> node_to_job;
    for (int i=0; i!= weakly_connected_components.size(); i++)
    {
        // cerr << "i" << i << endl;
        for (nid_t node : weakly_connected_components[i])
        {
            // cerr << "node in i" << node << endl;
            node_to_job[node] = i;
        }
    }
    return node_to_job;
}

std::vector<RebuildJob> divide_changelog_into_jobs(const std::unordered_map<nid_t, size_t>& node_to_job, const vector<unordered_set<nid_t>>& weakly_connected_components, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog)
{
    std::vector<RebuildJob> jobs(weakly_connected_components.size());
    // cerr << "jobs.size()" << jobs.size() << endl;   
    // cerr << "jobs[0].size() " << jobs[0].mappings.size() << endl;

    for (pair<gbwt::vector_type, gbwt::vector_type> change : gbwt_changelog)
    {
        // RebuildJob job;
        // cerr << "change.first.front() " << change.first.front() << endl;
        // cerr << "change.first.front() in node_id form. " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(change.first.front())) << endl;
        // cerr << "The node id from gbwt::Node::id: " << gbwt::Node::id(change.first.front()) << endl;
        // cerr << "does ndoe_to_job have the id when in change.first.front() in node_id form? " << node_to_job.at(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(change.first.front()))) << endl;
        // cerr << "node_to_job.contains(change.first.front())" << node_to_job.contains(change.first.front()) << endl;
        RebuildJob& cur_job = jobs[node_to_job.at(gbwt::Node::id(change.first.front()))];
        cur_job.mappings.push_back(change); //todo: make sure that I'm still correct in placing the gbwt::node in here, not the node id. I'm pretty sure I'm right.
        cur_job.total_size++;
        // cerr << "uh oh" << endl;
    }
    // cerr << "jobs[0].size() " << jobs[0].mappings.size() << endl;
    return jobs;
}

RebuildParameters set_update_gbwt_parameters(const int& threads, const bool& debug_print)
{
    RebuildParameters parameters;
    parameters.num_jobs = threads; //todo: add option for this to constructor.
    parameters.show_progress = debug_print;
    return parameters;
    //todo: implement Jouni's advice for the batch_size and sample_interval. Advice:
    /*
    * From Jouni:
    * The defaults should be fine in most cases.
    * If the threads are particularly long and memory usage is not a problem, you may want to increase the batch size to ~20x the length of the longest threads.
    */
    // parameters.batch_size = ???;
    // parameters.sample_interval = ???;
    // /// Maximum number of parallel construction jobs.
    // size_t num_jobs = 1;

    // /// Print progress information to stderr.
    // bool show_progress = false;

    // /// Size of the GBWT construction buffer in nodes.
    // gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;

    // /// Sample interval for locate queries.
    // gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
}

}
}
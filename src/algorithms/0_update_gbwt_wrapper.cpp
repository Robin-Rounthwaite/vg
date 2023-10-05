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


/// @brief Given a changelog of pair<node_ids before changes, node_ids after changes>,
/// will record the changes in the gbwt. The graph must have identical haplotype content
/// before/after. (so the graph changes may only cause disconnects between sequences that
/// were never explicitly connected by haplotypes).
/// @param gbwt_graph 
/// @param gbwt_changelog //todo: figure out if encoded in gbwt nodes or if in normal graph nodes.
/// @param gbwt 
/// @param threads 
/// @param debug_print 
/// @return 
gbwt::GBWT apply_gbwt_changelog(const gbwtgraph::GBWTGraph &gbwt_graph, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog, const gbwt::GBWT& gbwt, const int& threads, const bool& debug_print)
{
    // //todo: begin debug code
    cerr<< "************ does gbwt graph have node 358816? " << gbwt_graph.has_node(358816) << endl;;
    cerr << "gbwt.contains(358816): " << gbwt.contains(gbwt::Node::encode(358816, false)) << endl;
    // cerr << "gbwt.contains(12): " << gbwt.contains(gbwt::Node::encode(12, false)) << endl;
    // cerr << "gbwt.contains(16): " << gbwt.contains(gbwt::Node::encode(16, false)) << endl;
    // cerr << "gbwt.contains(9): " << gbwt.contains(gbwt::Node::encode(9, false)) << endl;
    // cerr << "gbwt.contains(17): " << gbwt.contains(gbwt::Node::encode(17, false)) << endl;
    // cerr << "gbwt.contains(6): " << gbwt.contains(gbwt::Node::encode(6, false)) << endl;
    // cerr << "gbwt.contains(18): " << gbwt.contains(gbwt::Node::encode(18, false)) << endl;
    // cerr << "gbwt_changelog overview: " << endl;
    // for (auto mapping : gbwt_changelog)
    // {
    //     cerr << "from: ";
    //     for (auto unedited_region : mapping.first)
    //     {
    //         cerr << gbwt::Node::id(unedited_region) << " ";
    //     }
    //     cerr << endl;
    //     cerr << "  to: ";
    //     for (auto edited_region : mapping.second)
    //     {
    //         cerr << gbwt::Node::id(edited_region) << " ";
    //     }
    //     cerr << endl;
    // }
    // //todo: end debug code

    vg::RebuildJob::mapping_type first_change = gbwt_changelog.front();
    // std::vector<vg::RebuildJob::mapping_type> debug_changelog;
    // debug_changelog.push_back(first_change);
    // cerr << "hanging change in apply_gbwt_changelog: " << first_change.first.back() << " " << first_change.second.back() << endl;

    // vector<unordered_set<nid_t>> components = gbwtgraph::weakly_connected_components(&_gbwt_graph);
    cerr << "weakly_connected_components" << endl;
    vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&gbwt_graph);
    // cerr << "size of components: " << components.size() << endl;
    // cerr << "gbwt_graph.has_node(1): " << gbwt_graph.has_node(1) << endl;
    gbwt_graph.has_node(1);
    // vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&_graph); // this was using the handlegraph algorithm, but I think we want the gbwt's view of the connected components.

    cerr << "get_node_to_job" << endl;
    std::unordered_map<nid_t, size_t> node_to_job = get_node_to_job(components);

    cerr << "divide_changelog_into_jobs" << endl;
    // std::vector<RebuildJob> jobs = divide_changelog_into_jobs(node_to_job, components, debug_changelog); 
    std::vector<RebuildJob> jobs = divide_changelog_into_jobs(node_to_job, components, gbwt_changelog); 

    cerr << "set_update_gbwt_parameters" << endl;
    RebuildParameters rebuild_parameters = set_update_gbwt_parameters(threads, debug_print);
    
    cerr << "rebuild_gbwt" << endl;
    gbwt::GBWT output_gbwt = rebuild_gbwt(gbwt, jobs, node_to_job, rebuild_parameters);
    // //todo: remove temporary non-parallelized rebuild_gbwt.
    // cerr << "output_gbwt.contains(12): " << output_gbwt.contains(gbwt::Node::encode(12, false)) << endl;
    // cerr << "output_gbwt.contains(16): " << output_gbwt.contains(gbwt::Node::encode(16, false)) << endl;
    // cerr << "output_gbwt.contains(9): " << output_gbwt.contains(gbwt::Node::encode(9, false)) << endl;
    // cerr << "output_gbwt.contains(17): " << output_gbwt.contains(gbwt::Node::encode(17, false)) << endl;
    // cerr << "output_gbwt.contains(6): " << output_gbwt.contains(gbwt::Node::encode(6, false)) << endl;
    // cerr << "output_gbwt.contains(18): " << output_gbwt.contains(gbwt::Node::encode(18, false)) << endl;

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
            if (node == 358816)
            {
                cerr << "OI! 358816 is about to be put in node_to_job." << endl;
            }
            // cerr << "node in i" << node << endl;
            node_to_job[node] = i;
        }
    }
    return node_to_job;
}

std::vector<RebuildJob> divide_changelog_into_jobs(const std::unordered_map<nid_t, size_t>& node_to_job, const vector<unordered_set<nid_t>>& weakly_connected_components, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog)
{
    std::vector<RebuildJob> jobs(weakly_connected_components.size());
    cerr << "jobs.size() " << jobs.size() << endl;   
    cerr << "weakly_connected_components.size() " << weakly_connected_components.size() << endl;   
    cerr << "jobs[0].size() " << jobs[0].mappings.size() << endl;
    cerr << "node_to_job.at(358816) " << node_to_job.at(358816) << endl;
    cerr << "jobs[node_to_job.at(358816)].mappings.size() " << jobs[node_to_job.at(358816)].mappings.size() << endl;
    jobs[node_to_job.at(358816)].mappings.size();
    cerr << "jobs[1].size() " << jobs[1].mappings.size() << endl;
    cerr << "jobs[234234].size() " << jobs[234234].mappings.size() << endl;
    for (pair<gbwt::vector_type, gbwt::vector_type> change : gbwt_changelog)
    {
        // cerr << "running jobs[node_to_job.at(358816)].mappings.size(): " <<endl;
        // jobs[node_to_job.at(358816)].mappings.size();

        // RebuildJob job;
        // cerr << "change.first.front() " << change.first.front() << endl;
        // // cerr << "change.first.front() in node_id form. " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(change.first.front())) << endl;
        // cerr << "The node id from gbwt::Node::id: " << gbwt::Node::id(change.first.front()) << endl;
        // // cerr << "does ndoe_to_job have the id when in change.first.front() in node_id form? " << node_to_job.at(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(change.first.front()))) << endl;
        // cerr << "node_to_job.contains(change.first.front()) " << !(node_to_job.find(change.first.front()) == node_to_job.end()) << endl;
        // cerr << "contents of node_to_job, of size " << node_to_job.size() << endl;
        // for (auto it = node_to_job.begin(); it != node_to_job.end(); it++)
        // {
        //     cerr << it->first << " " << it->second << endl;
        // }
        try
        {
           jobs[node_to_job.at(gbwt::Node::id(change.first.front()))];
        }
        catch (...)
        {
            cerr << " this is the buggy node: " << gbwt::Node::id(change.first.front()) << endl;
            exit(1);
        }
        RebuildJob& cur_job = jobs[node_to_job.at(gbwt::Node::id(change.first.front()))];
        // RebuildJob& cur_job = jobs[node_to_job.at(change.first.front())];
        cur_job.mappings.push_back(change); //todo: make sure that I'm still correct in placing the gbwt::node in here, not the node id. I'm pretty sure I'm right.
        cur_job.total_size++;
    }
    // cerr << "jobs[0].size() " << jobs[0].mappings.size() << endl;
    return jobs;
}

RebuildParameters set_update_gbwt_parameters(const int& threads, const bool& debug_print)
{
    RebuildParameters parameters;
    parameters.num_jobs = threads;
    parameters.show_progress = true; //todo: remove thisi line, uncomment below.
    // parameters.show_progress = debug_print;
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
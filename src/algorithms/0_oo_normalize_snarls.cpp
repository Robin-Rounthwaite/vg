#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?

#include "0_oo_normalize_snarls.hpp"
#include "0_snarl_sequence_finder.hpp"
#include <string>
#include <tuple>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

#include <gbwtgraph/gbwtgraph.h>

#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../msa_converter.hpp"
#include "../snarls.hpp"
#include "../vg.hpp"

#include "../types.hpp"
#include "extract_containing_graph.hpp"

#include "multipath_mapper.hpp"

#include <bdsg/packed_graph.hpp>


/*
TODO: allow for snarls that have haplotypes that begin or end in the middle of the snarl

TODO: allow normalization of multiple adjacent snarls in one combined realignment.

TODO: test that extract_gbwt haplotypes successfully extracts any haplotypes that start/end in the middle of
TODO:    the snarl.
*/

// todo: add cyclic snarls to the ones to skip, if cyclic snarls turns out to be frequent. Right now, I want to know when I have a cyclic snarl.

int _big_snarl_alignment_job = 900;

namespace vg {
namespace algorithms{
/**
 * To "normalize" a snarl, SnarlNormalizer extracts all the sequences in the snarl as
 * represented in the gbwt, and then realigns them to create a replacement snarl. 
 * This process hopefully results in a snarl with less redundant sequence, and with 
 * duplicate variation combined into a single variant.
*/

SnarlNormalizer::SnarlNormalizer(MutablePathDeletableHandleGraph &graph,
                                 const gbwt::GBWT &gbwt,
                                 const gbwtgraph::GBWTGraph &gbwt_graph,
                                 const int max_handle_size,
                                 const int max_region_size,
                                 const int threads,
                                 const int max_alignment_size, /*= MAX_INT*/
                                 const string path_finder, /*= "GBWT"*/
                                 const string alignment_algorithm, /*= "sPOA"*/
                                 const bool disable_gbwt_update, /*= false*/
                                 const bool debug_print /*= false*/)
: _graph(graph), _gbwt(gbwt), _gbwt_graph(gbwt_graph), _max_alignment_size(max_alignment_size),
      _max_handle_size(max_handle_size), _max_region_size(max_region_size), _threads(threads), _path_finder(path_finder),
       _alignment_algorithm(alignment_algorithm), _disable_gbwt_update(disable_gbwt_update), _debug_print(debug_print){}

/// @brief note: assumes that each region is pair<leftmost_id, rightmost_id>, i.e. that 
///     backwards is false (which it is fo return values of get_ and split_normalize_regions).
/// @param split_normalize_regions 
std::vector<vg::RebuildJob::mapping_type> SnarlNormalizer::parallel_normalization(vector<pair<id_t, id_t>> split_normalize_regions)
{
    // cerr << "list of all to-delete handles:" << endl;
    // for (auto deletable : _nodes_to_delete)
    // {
    //     cerr << deletable << endl;
    //     cerr << "sequence: " << _graph.get_sequence(_graph.get_handle(deletable)) << endl;
    // }
    
    //todo: add back the exhaustive path finder option?
    assert(_path_finder=="GBWT");

    //normalized snarls will contain all the normalized snarls before inserting them into the graph.
    // This is because inserting the snarls in a parallel fashion is difficult. 
    // (and would probably involve locking individual subgraphs using an adaptation of 
    // Adam's lock code)
    // each entry in normalized_snarls will be the arguments necessary to call integrate_snarl.
    // integrate_snarl(SubHandleGraph &old_snarl,  const HandleGraph &to_insert_snarl, vector<pair<step_handle_t, step_handle_t>>& embedded_paths,  const id_t source_id, const id_t sink_id, const bool backwards)
    //todo: change to make more memory efficient? 
    //todo: I could also at least re-derive the SubHandleGraph if I needed to.
    // vector< pair< vg::VG, pair<id_t,id_t> > > normalized_snarls;
    vector< tuple< SubHandleGraph, shared_ptr<MutablePathDeletableHandleGraph>, std::vector<std::pair<vg::step_handle_t, vg::step_handle_t>>, id_t, id_t, bool, vector<pair<gbwt::vector_type, string>> >> normalized_snarls;
    

    // //todo: debug_code
    _debug_print=true;
    // split_normalize_regions.clear();
    // split_normalize_regions.push_back(make_pair(1898335, 1898346));
    // split_normalize_regions.push_back(make_pair(2624390, 2624421));
    int num_snarls_normalized = 0;
        // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    cerr << "starting clock for internal normalize stats." << endl;



    //todo: begin debug_code:
    // split_normalize_regions.clear();
    // split_normalize_regions.push_back(make_pair(2369282, 2369288));
    // split_normalize_regions.push_back(make_pair(2555912, 2555931));
    // // vector<pair<id_t, id_t>> debug_split_normalize_regions;
    // int count = 0;
    // for (auto i = split_normalize_regions.begin(); i!= split_normalize_regions.end(); i++)
    // {
    //     count++;
    //     if (i->first < 2556000 && i->first > 2555000) 
    //     {
    //         cerr << count << endl;
    //         cerr << i->first << " " << i->second << endl;
    //     }
    //     if (i->first == 2555912 && i->second==2555931)
    //     {
    //         cerr << "count is " << count << endl;
    //         split_normalize_regions.push_back(*(--i));
    //         split_normalize_regions.push_back(*(i));
    //     }
    //     break;
    // }
    //todo: end debug_code:


    omp_set_num_threads(_threads);
    #pragma omp parallel for
    for (auto region : split_normalize_regions)
    {
        #pragma omp critical(print_progress)
        {
        if (num_snarls_normalized%10000 == 0)
        // if (num_snarls_normalized%1 == 0)
        {
            auto cur_time = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> elapsed = cur_time - start;
            
            cerr << "normalizing " << num_snarls_normalized+1 << "/" << split_normalize_regions.size() << " regions by time " << elapsed.count() << endl;
            
        }
        num_snarls_normalized++;
        }
        pair<bool, bool> sequence_added_because_empty_node = make_pair(false, false);

        // _debug_print=true;
        if (_debug_print)
        {
            cerr << "region: " << region.first << " " << region.second << endl;
            cerr << "about to extract_subgraph" << endl;
        }
        // cerr << "does the graph contain node 18? " << _graph.has_node(18) << endl;
        SubHandleGraph snarl = extract_subgraph(_graph, region.first, region.second);

        //get original snarl size for comparison stats
        int original_snarl_size = 0;
        snarl.for_each_handle([&](handle_t handle){
            original_snarl_size += snarl.get_sequence(handle).size();
        });

        // check that snarl is normalizable:
        bool passes_normalize_tests = test_snarl(snarl, region, original_snarl_size);
        // if (passes_normalize_tests == false) { cerr <<  "failed normalize tests" << endl; continue; }
        if (passes_normalize_tests == false) { if(_debug_print){cerr <<  "failed normalize tests" << endl;} continue; }

        if (_debug_print)
        {
            cerr << "about to extract_haplotypes" << endl;
        }
        // stop_inclusive ensures that the step_handle_t doesn't follow the standard
        // convention of going one past the furthest handle on the path that we want to
        // look at. Instead, it points to the furthest handle itself That way we don't get
        // steps pointing to nodes outside the region of this current normalization, so we
        // don't risk edits from other snarls messing up that step. 
        bool stop_inclusive = true;
        //declare arguments to be filled by extract_haplotypes:
            // haplotypes is of format:
            // 0: a set of all the haplotypes which stretch from source to sink, in string format.
            //   - it's a set, so doesn't contain duplicates
            // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
            // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
        tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<id_t>> haplotypes;
        vector<pair<step_handle_t, step_handle_t>> embedded_paths;
        vector<pair<gbwt::vector_type, string>> source_to_sink_gbwt_paths;
        //extract the haplotypes:
        extract_haplotypes(snarl, region, sequence_added_because_empty_node, haplotypes, embedded_paths, source_to_sink_gbwt_paths, stop_inclusive);

        if (_debug_print)
        {
            cerr << "about to test haplotypes" << endl;
        }
        // further checks that snarl is normalizable:
        bool passes_haplotype_tests = test_haplotypes(haplotypes, region, original_snarl_size);
        if (passes_haplotype_tests == false) { cerr <<  "failed haplotype tests" << endl; continue; }

        if (_debug_print)
        {
            cerr << "about to align haplotypes" << endl;
        }
        // align haplotypes to form new snarl:
        shared_ptr<MutablePathDeletableHandleGraph> new_snarl;
        // bdsg::PackedGraph new_snarl_graph;
        // MutablePathDeletableHandleGraph* new_snarl = &new_snarl_graph;
        if (_alignment_algorithm == "TCoffee")
        {
            cerr << "TCoffee is worse than sPOA in nearly every way, and is no longer supported." << endl; //todo, delete this option."
            exit(1);
            // new_snarl = align_source_to_sink_haplotypes(get<0>(haplotypes));
        }
        else if (_alignment_algorithm == "sPOA")
        {
            cerr << "get<0>(haplotypes), new_snarl, false) " << get<0>(haplotypes).size() << endl;
            cerr << "new_snarl interaction: " << endl;
            new_snarl->for_each_handle([&](handle_t handle){
                cerr << new_snarl->get_id(handle) << endl;
            });

            cerr << "about to run poa_source_to_sink_haplotypes. " << endl;
            new_snarl = poa_source_to_sink_haplotypes(get<0>(haplotypes), false);
            cerr << "finished running poa_source_to_sink_haplotypes. " << endl;
            new_snarl->for_each_handle([&](handle_t handle){
                cerr << new_snarl->get_id(handle) << endl;
            });
            cerr << "done checking handles after poa fxn" << endl;
            if (!new_snarl)
            {
                //note: this snippet should never run, because it's also handled by the "if (hap.size() > max_spoa_length)" condition a in test_haplotypes.
                cerr << "ERROR: poa source to sink haplotypes run unsuccessful because haps too long. But since this situation is handled in a previous check, you should never see this message." << endl;
                exit(1);
            }
        }
        else
        {
            cerr << "error:[vg normalize] _alignment_algorithm variable must be set as either T-Coffee or sPOA." << endl;
            exit(1);
        }

        //remove any excess characters that were added for alignment (i.e. if the
        //get_parallel_normalize process introduced an empty node for source or sink in
        //the original snarl.)
        pair<vector<handle_t>, vector<handle_t>> to_insert_snarl_defining_handles =
            debug_get_sources_and_sinks(*new_snarl);
        if (sequence_added_because_empty_node.first)
        {
            string new_node_seq = new_snarl->get_sequence(to_insert_snarl_defining_handles.first.back());
            id_t new_node_id = new_snarl->get_id(to_insert_snarl_defining_handles.first.back());
            // cerr << "1 contents of replace node using sequence: " << new_node_id << " " <<  new_node_seq << " " << new_node_seq.substr(1) << endl;
            replace_node_using_sequence(new_node_id, new_node_seq.substr(1), *new_snarl);
            sequence_added_because_empty_node.first = false;
        }
        if (sequence_added_because_empty_node.second)
        {
            string new_node_seq = new_snarl->get_sequence(to_insert_snarl_defining_handles.second.back());
            id_t new_node_id = new_snarl->get_id(to_insert_snarl_defining_handles.second.back());
            // cerr << "2 contents of replace node using sequence: " << new_node_id << " new_node_seq: " << new_node_seq << " new_node_seq.substr(1): " << new_node_seq.substr(1) << endl;
            replace_node_using_sequence(new_node_id, new_node_seq.substr(0, new_node_seq.size() - 1), *new_snarl);
            sequence_added_because_empty_node.second = false;
        }
        
        // get new snarl size for comparison stats
        int new_snarl_size = 0;
        new_snarl->for_each_handle([&](const handle_t handle) {
            new_snarl_size += new_snarl->get_sequence(handle).size();
        });
        

        //make sure all the handles are the proper size of <=_maximum_handle_size.
        // force_maximum_handle_size(*new_snarl);

        if (_debug_print)
        {
            cerr << "about to store the normalized snarl for later." << endl;
        }

        //record snarl size change statistic

        #pragma omp critical(save_snarl)
        {
        _snarl_size_changes[make_pair(region.first, region.second)] = make_pair(original_snarl_size, new_snarl_size);
        // std::vector<std::tuple<vg::SubHandleGraph, vg::VG, std::vector<std::pair<vg::step_handle_t, vg::step_handle_t>>, vg::id_t, vg::id_t, bool, std::vector<std::pair<gbwt::vector_type, std::string>>>> normalized_snarls
        normalized_snarls.emplace_back(make_tuple(snarl, new_snarl, embedded_paths, region.first, region.second, false, source_to_sink_gbwt_paths));
        }
        // pair<handle_t, handle_t> new_left_right = integrate_snarl(snarl, new_snarl, embedded_paths, region.first, region.second, false);
    }

    auto align_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> align_elapsed = align_time - start;
    cerr << "all snarls normalized at time " << align_elapsed.count() << ". About to start integrating snarls." << endl;

    //integrate all the normalized snarls formed in the parallel loop above.
    for (auto snarl : normalized_snarls)
    {
        
        if (_debug_print)
        {
            cerr << "======================about to integrate snarl " << get<3>(snarl) << " " << get<4>(snarl) << "======================" << endl;
        }

        // cerr << "integration" << endl;
        pair<handle_t, handle_t> new_left_right = integrate_snarl(get<0>(snarl), *get<1>(snarl), get<2>(snarl), get<3>(snarl), get<4>(snarl), get<5>(snarl));
        // make a subhandlegraph of the normalized snarl to find the new gbwt paths in the graph.
        // cerr << "extraction" << endl;
        SubHandleGraph integrated_snarl = extract_subgraph(_graph, _graph.get_id(new_left_right.first), _graph.get_id(new_left_right.second));
        // cerr << "log changes" << endl;
        log_gbwt_changes(get<6>(snarl), integrated_snarl);

    }
    auto integration_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> integration_elapsed = integration_time - align_time;
    std::chrono::duration<double> total_elapsed = integration_time - start;
    cerr << "Integration took " << integration_elapsed.count() << " seconds. Total time: " << total_elapsed.count() << endl;

    print_parallel_statistics();
    return _gbwt_changelog;
}

void SnarlNormalizer::print_parallel_statistics()
{
    cerr << "---------summary of skipped regions:---------" << endl;
    // vector<int> _sizes_of_snarls_skipped_because_gbwt_misses_handles;
    if (_snarls_skipped_because_gbwt_misses_handles.size() >0)
    {
        cerr << "number of regions skipped because gbwt misses handles: " << _snarls_skipped_because_gbwt_misses_handles.size() << endl; 
    }

    if (_snarls_skipped_because_gbwt_misses_edges.size() >0)
    {
        cerr << "number of regions skipped because gbwt misses edges: " << _snarls_skipped_because_gbwt_misses_edges.size() << endl; 
    }

    // vector<int> _sizes_of_snarls_skipped_because_cyclic;
    if (_snarls_skipped_because_cyclic.size() >0)
    {
        cerr << "number of regions skipped because cyclic: " << _snarls_skipped_because_cyclic.size() << endl; 
    }

    // vector<int> _sizes_of_snarls_skipped_because_haplotype_ends_in_middle;
    if (_snarls_skipped_because_haplotype_ends_in_middle.size() >0)
    {
        cerr << "number of regions skipped because haplotype ends in middle: " << _snarls_skipped_because_haplotype_ends_in_middle.size() << endl; 
    }

    // vector<int> _sizes_of_snarls_skipped_because_alignment_too_many_threads;
    if (_snarls_skipped_because_alignment_too_many_threads.size() >0)
    {
        cerr << "number of regions skipped because alignment too many threads: " << _snarls_skipped_because_alignment_too_many_threads.size() << endl; 
    }

    // vector<int> _sizes_of_snarls_skipped_because_haps_too_long_for_spoa;
    if (_snarls_skipped_because_haps_too_long_for_spoa.size() >0)
    {
        cerr << "number of regions skipped because haps too long for spoa: " << _snarls_skipped_because_haps_too_long_for_spoa.size() << endl; 
    }

    cerr << "---------summary of normalized regions:---------" << endl;
    int num_top_snarls_tracked = 10; // hardcoded for debugging convenience. //todo: change?

    // Detect the top best (shrinking-size, not actually necessarily best) snarl changes:
    auto best_compare = [](const pair<pair<id_t, id_t>, pair<int, int>> left, const pair<pair<id_t, id_t>, pair<int, int>> right) 
    {
        return (right.second.second - right.second.first) >= (left.second.second - left.second.first);
    };
    vector<pair<pair<id_t, id_t>, pair<int, int>>> best_snarl_changes;


    for (pair<pair<id_t, id_t>, pair<int, int>> region : _snarl_size_changes)
    {

        if (best_snarl_changes.size() < num_top_snarls_tracked)
        {
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
        }
        else if ((region.second.second - region.second.first) < (best_snarl_changes.back().second.second - best_snarl_changes.back().second.first)) //todo: check that this is proper comparison.
        {
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
            best_snarl_changes.pop_back();
        }
    }


    // general statistics
    int size_of_all_shrinking_snarls_pre_norm = 0;
    int size_of_all_shrinking_snarls_post_norm = 0;
    int size_of_all_growing_snarls_pre_norm = 0;
    int size_of_all_growing_snarls_post_norm = 0;
    
    // Detect the top worst (growing-size, not actually necessarily worst) snarl changes:
    int snarls_that_grow_after_norm = 0;
    int snarls_that_shrink_after_norm = 0;
    auto worst_compare = [](const pair<pair<id_t, id_t>, pair<int, int>> left, const pair<pair<id_t, id_t>, pair<int, int>> right) 
    {
        return (right.second.second - right.second.first) <= (left.second.second - left.second.first);
    };
    vector<pair<pair<id_t, id_t>, pair<int, int>>> worst_snarl_changes;
    for (pair<pair<id_t, id_t>, pair<int, int>> region : _snarl_size_changes)
    {
        
        if (region.second.second - region.second.first > 0)
        {
            snarls_that_grow_after_norm += 1;
            size_of_all_growing_snarls_pre_norm += region.second.first;
            size_of_all_growing_snarls_post_norm += region.second.second;
        }
        else if (region.second.second - region.second.first < 0)
        {
            snarls_that_shrink_after_norm += 1;
            size_of_all_shrinking_snarls_pre_norm += region.second.first;
            size_of_all_shrinking_snarls_post_norm += region.second.second;
        }
        // cerr << "deciding whether or not to insert following value: " << region.second.first << " , " << region.second.second << endl;
        if (worst_snarl_changes.size() < num_top_snarls_tracked)
        {
            // cerr << "inserted because the worst_snarl_changes isn't big enough." << endl;
            worst_snarl_changes.insert(upper_bound(worst_snarl_changes.begin(), worst_snarl_changes.end(), region, worst_compare), region);
        }
        else if ((region.second.second - region.second.first) > (worst_snarl_changes.back().second.second - worst_snarl_changes.back().second.first))
        {
            // cerr << "inserted and replaced an item. shrinkage of Item replaced: " << worst_snarl_changes.back().second.second - worst_snarl_changes.back().second.first << endl;
            // cerr << "Shrinkage of item added: " << region.second.second - region.second.first << endl;
            worst_snarl_changes.insert(upper_bound(worst_snarl_changes.begin(), worst_snarl_changes.end(), region, worst_compare), region);
            worst_snarl_changes.pop_back();
        }
    }
    cerr << "total snarls that shrink in size: " << snarls_that_shrink_after_norm << endl;
    cerr << "percent change of snarls that shrink in size: " << (((double)size_of_all_shrinking_snarls_post_norm - (double)size_of_all_shrinking_snarls_pre_norm)/(double)snarls_that_shrink_after_norm)*100 << "%" << endl;  
    cerr << "top " << num_top_snarls_tracked << " snarls *reduction* in size (probably desirable): " << endl;
    cerr << "leftmost_id  rightmost_id  original_size  normalized_size  (change_in_size;  percent change)" << endl;
    
    // for (auto i = _snarl_size_changes.end(); i != next(_snarl_size_changes.end(), -num_top_snarls_tracked); i--)
    for (auto region : best_snarl_changes)
    {
        cerr << region.first.first << "   " << region.first.second << "   " << region.second.first << "   " << region.second.second << "   (" << region.second.second - region.second.first << ";" << "   " << ((((double)region.second.second - (double)region.second.first)/(double)region.second.first)*100.0) << "%)" << endl;
    }

    cerr << "total snarls that grow in size: " << snarls_that_grow_after_norm << endl;
    cerr << "percent change of snarls that grow in size: " << (((double)size_of_all_growing_snarls_post_norm - (double)size_of_all_growing_snarls_pre_norm)/(double)snarls_that_grow_after_norm)*100 << "%" << endl;  
    cerr << "top " << num_top_snarls_tracked << " snarls *increase* in size (probably undesirable): " << endl;
    cerr << "leftmost_id  rightmost_id  original_size  normalized_size  (change_in_size;  percent change)" << endl;
    for (auto region : worst_snarl_changes)
    {
        cerr << region.first.first << "   " << region.first.second << "   " << region.second.first << "   " << region.second.second << "   (" << region.second.second - region.second.first << ";" << "   " << ((((double)region.second.second - (double)region.second.first)/(double)region.second.first)*100.0) << "%)" << endl;
    }

}

/// @brief 
/// @param snarl 
/// @param region 
/// @param snarl_size 
/// @return bool snarl_passes: if the snarl is normalizable, returns true. Else, false. 
bool SnarlNormalizer::test_snarl(const SubHandleGraph& snarl, const pair<id_t, id_t>& region, const int snarl_size)
{
    //make sure that all handles in the snarl are represented by the gbwt.
    bool all_handles_in_gbwt = true;
    // bool all_edges_in_gbwt = true;
    // try 
    // {
    //     snarl.for_each_handle([&](handle_t handle){
    //         cerr << "handle exists: " << snarl.get_id(handle) << endl;
    //     }, true);

    // }
    // catch (...)
    // {
    //     cerr << "error in for_each_handle. Here is the region: " << region.first << " " << region.second << endl;
    // }
    // cerr << endl;
    snarl.for_each_handle([&](handle_t handle){
        // cerr << "handle exists in second for each handle: " << snarl.get_id(handle) << endl;
        // cerr << "snarl.get_id(handle)" << " " << snarl.get_id(handle) << endl;
        const gbwt::SearchState handle_state = _gbwt_graph.get_state(_gbwt_graph.get_handle(snarl.get_id(handle)));
        // cerr << "_gbwt_graph.get_handle(snarl.get_id(handle))" << " " << _gbwt_graph.get_handle(snarl.get_id(handle)) << endl;
        // cerr << "_gbwt_graph.get_bd_state(_gbwt_graph.get_handle(snarl.get_id(handle)))" << " " << _gbwt_graph.get_bd_state(_gbwt_graph.get_handle(snarl.get_id(handle))) << endl;
        // cerr << "does the gbwt contain this node? " << _gbwt.contains(gbwt::Node::encode(snarl.get_id(handle), false)) << endl;
        // cerr << "does the gbwt_graph contian this node? " << _gbwt_graph.has_node(snarl.get_id(handle)) << endl;

        if (handle_state.empty())
        {
            if (_debug_print)
            {
                cerr << "handle " << snarl.get_id(handle) << " not in gbwt." << endl;
            }

            // cerr << "debug state is empty." << endl;
            all_handles_in_gbwt = false;
            return;
        }

        // set<id_t> right_neighbors;
        // _graph.follow_edges(handle, false, [&](handle_t next_handle)
        // {
        //     right_neighbors.emplace(_graph.get_id(next_handle));
        // });

        // if (snarl.get_id(handle) == 2624390)
        // {
        //     cerr << "at handle 2624390, looking at right neighbors: " << endl;
        //     for (auto n : right_neighbors)
        //     {
        //         cerr << n << " ";
        //     }
        //     cerr << endl;
        // }
        //todo: uncomment if I decide to check that all edges in the graph are in the gbwt:
        // set<id_t> gbwt_right_neighbors;


        // _gbwt_graph.follow_paths(handle_state,
        // [&](const gbwt::SearchState next_search) -> bool {
        //     gbwt_right_neighbors.emplace(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)));
        //     // if (right_neighbors.find(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node))) == right_neighbors.end())
        //     // {
        //     //     if (snarl.get_id(handle) == 2624390)
        //     //     {
                    
        //     //         cerr << "looking for node id " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << " equivalence: " << (right_neighbors.find(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(handle_state.node))) == right_neighbors.end()) << endl;
        //     //         // cerr << "nodes in right_neighbors: " << endl;
        //     //         // for (auto node : right_neighbors)
        //     //         // {
        //     //         //     cerr << node << " ";
        //     //         // }
        //     //         // cerr << endl;
        //     //     }
                
        //     //     if (_debug_print)
        //     //     {
        //     //         cerr << "handle " << snarl.get_id(handle) << " not connect to " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << " in gbwt." << endl;
        //     //     }
        //     //     all_edges_in_gbwt = false;
        //     // }
        //     return true;
        // });
        // _graph.follow_edges(handle, false, [&](handle_t next_handle)
        // {
        //     if (gbwt_right_neighbors.find(_graph.get_id(next_handle)) ==  gbwt_right_neighbors.end())
        //     {
        //         all_edges_in_gbwt = false;
        //     }
        // });
        //todo: end uncoment.

    // cerr << "all_edges_in_gbwt " << all_edges_in_gbwt << endl << endl << endl; 
    }, true);
    if (!all_handles_in_gbwt)
    {
        if (_debug_print)
        {
            cerr << "test failed because !all_handles_in_gbwt. " << endl;
        }
        #pragma omp critical(error_update_1)
        {
        _sizes_of_snarls_skipped_because_gbwt_misses_handles.push_back(snarl_size);
        _snarls_skipped_because_gbwt_misses_handles.push_back(region);
        }
        // cerr << "!all handles in gbwt!" << "number of snarls skipped for this reason is: " << _snarls_skipped_because_gbwt_misses_handles.size() << endl;
        return false;
    }

    // if (!all_edges_in_gbwt)
    // {
    //     if (_debug_print)
    //     {
    //         cerr << "test failed because !all_edges_in_gbwt. " << endl;
    //     }
    //     #pragma omp critical(error_update_1_5)
    //     {
    //     _sizes_of_snarls_skipped_because_gbwt_misses_edges.push_back(snarl_size);
    //     _snarls_skipped_because_gbwt_misses_edges.push_back(region);
    //     }
    //     // cerr << "!all handles in gbwt!" << "number of snarls skipped for this reason is: " << _snarls_skipped_because_gbwt_misses_handles.size() << endl;
    //     return false;
    // }

    // check for cyclic snarls
    if (!handlealgs::is_acyclic(&snarl)) {
        if (_debug_print)
        {
            cerr << "test failed because !handlealgs::is_acyclic(&snarl). " << endl;
        }
        #pragma omp critical(error_update_2)
        {
        _sizes_of_snarls_skipped_because_cyclic.push_back(snarl_size);
        _snarls_skipped_because_cyclic.push_back(region);
        }
        return false;
    }

    return true;
}

bool SnarlNormalizer::test_haplotypes(const tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<id_t>>& haplotypes, const pair<id_t, id_t>& region, const int original_snarl_size)
{
    // this if statement only permits snarls that have no haplotypes beginning or 
    //      ending in the middle of the snarl. 
    //todo: Get rid of this once I can align these types of haps.
    if (!(get<1>(haplotypes).empty()))
    {
        #pragma omp critical (test_haplotypes_1)
        {
        _sizes_of_snarls_skipped_because_haplotype_ends_in_middle.push_back(original_snarl_size);
        _snarls_skipped_because_haplotype_ends_in_middle.push_back(region);
        }
        return false;
    }
    // also, limit the number of haplotypes to be aligned, since normalize slows down 
    // greatly at higher thread counts. //todo: identify best _max_alignment_size for spoa/abpoa
    if (get<0>(haplotypes).size() > _max_alignment_size)
    {
        #pragma omp critical (test_haplotypes_2)
        {
        _sizes_of_snarls_skipped_because_alignment_too_many_threads.push_back(original_snarl_size);
        _snarls_skipped_because_alignment_too_many_threads.push_back(region);
        }
        return false;
    }

    int max_spoa_length = 750; // somewhere between 500-1000 bases, sPOA starts to struggle. //todo: That's why I'll eventually want abPOA to take over.
    for (string hap : get<0>(haplotypes))
    {
        if (hap.size() > max_spoa_length)
        {
            #pragma omp critical (test_haplotypes_3)
            {
            _sizes_of_snarls_skipped_because_haps_too_long_for_spoa.push_back(original_snarl_size);
            _snarls_skipped_because_haps_too_long_for_spoa.push_back(region);
            }
            return false;
        }
    }
    return true;
}

//todo: make the return value void, and instead be returned by reference in arguments passed to extract_haplotypes.
void SnarlNormalizer::extract_haplotypes(const SubHandleGraph& snarl, const pair<id_t, id_t>& region, 
        pair<bool, bool>& sequence_added_because_empty_node, 
        tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<id_t>>& haplotypes, 
        vector<pair<step_handle_t, step_handle_t>>& embedded_paths, 
        vector<pair<gbwt::vector_type, string>>& source_to_sink_gbwt_paths, 
        const bool stop_inclusive/*=false*/)
{
    bool backwards = false;

    // if (_debug_print == true)
    // {
    //     cerr << "all nodes in snarl according to extract_haplotypes:" << endl;

    //     snarl.for_each_handle([&](handle_t handle)
    //     {
    //         cerr << snarl.get_id(handle) << " " << snarl.get_sequence(handle) << endl;
    //         cerr << "lefts: ";
    //         snarl.follow_edges(handle, true, [&](handle_t left){
    //             cerr << snarl.get_id(left) << " ";
    //         });
    //         cerr << endl;
    //         cerr << "rights: ";
    //         snarl.follow_edges(handle, false, [&](handle_t right){
    //             cerr << snarl.get_id(right) << " ";
    //         });
    //         cerr << endl;
    //     });
    //     cerr << "end all nodes in snarl according to extract_haplotypes." << endl;


        // //todo: debug_print:
        // gbwt::SearchState debug_state = _gbwt_graph.get_state(_gbwt_graph.get_handle(2555915));
        // _gbwt_graph.follow_paths(debug_state,
        //                          [&](const gbwt::SearchState next_search) -> bool {
        //                              cerr << "(directly from the query node of 2555915): adjacent handles to state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(debug_state.node)) << " is " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
        //                              return true;
        //                          });

        // gbwt::SearchState debug_state_2 = _gbwt_graph.get_state(_gbwt_graph.get_handle(2555917));
        // _gbwt_graph.follow_paths(debug_state_2,
        //                          [&](const gbwt::SearchState next_search) -> bool {
        //                              cerr << "(directly from the query node of 2555917): adjacent handles to state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(debug_state.node)) << " is " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
        //                              return true;
        //                          });
        // //todo: end debug_print

    // }
    
    // extract threads.
    // haplotypes is of format:
    // 0: a set of all the haplotypes which stretch from source to sink, in string format.
    //   - it's a set, so doesn't contain duplicates
    // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
    // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, region.first, region.second, false);

    //todo: here is where the exhaustive path finder would be used, if it was working.
    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<id_t>>
        gbwt_haplotypes = sequence_finder.find_gbwt_haps();

    unordered_set<id_t> nodes_in_snarl;
    snarl.for_each_handle([&](handle_t handle){
        nodes_in_snarl.emplace(snarl.get_id(handle));
    });

    //check that all handles touched by find_gbwt_haps is equivalent to all the 
    //  handles in the subgraph:
    int debug_sequence_not_in_gbwt = 0;
    for(id_t node_id : nodes_in_snarl)
    {
        if (get<2>(gbwt_haplotypes).find(node_id) == nodes_in_snarl.end()){
            const gbwt::BidirectionalState debug_state = _gbwt_graph.get_bd_state(_gbwt_graph.get_handle(node_id));
            if (debug_state.empty())
            {
                //i.e., there's a node in the graph that isn't in the gbwt.
                cerr << "since this situation is handled in a previous check, you should never see this message." << endl;
                exit(1);
            }
            else
            {
                cerr << "ERROR: the node " << node_id << " in the graph at the snarl with leftmost node id " << region.first << " and rightmost node id " << region.second << " is not touched by find_gbwt_haps(), but the node still exists in the gbwt. This is not supposed to happen. Exiting the program." << endl;
                exit(1);
            }
        }
    }

    // Convert the haplotypes from vector<handle_t> format to string format.
    get<0>(haplotypes) = format_handle_haplotypes_to_strings(_graph, _gbwt_graph, get<0>(gbwt_haplotypes));

    // //todo: debug_print
    // cerr << "all haploytypes as seen in extract_haplotypes." << endl;
    // for (auto hap : get<0>(haplotypes))
    // {
    //     cerr << hap << endl;
    // }
    // cerr << "all haploytypes as seen in extract_haplotypes' non-source_to_sink." << endl;
    // for (auto hap : get<1>(haplotypes))
    // {
    //     for (handle_t handle : hap)
    //     {
    //         cerr << _graph.get_sequence(handle) << " ";

    //     }
    //     cerr << endl;
    // }
    // //todo: end debug_print
    // cerr << "hap original text: " << *get<0>(haplotypes).begin() << endl;
    // check to see if the leftmost or rightmost node is empty. If so, treat the blank
    // node as containing a character "A". (this is important for dealing with how
    // get_parallel_normalize_regions sometimes adds a blank node when isolating adjacent
    // snarls for parallelization).
    if (snarl.get_sequence(snarl.get_handle(region.first)).size() == 0)
    {
        // //check that this is a node we expected to find as empty.
        // if (_nodes_to_delete.find(region.first) == _nodes_to_delete.end())
        // {
        //     cerr << "ERROR: a source or sink node for a snarl is empty in an unexpected manner (i.e., it wasn't introduced into the graph during get_parallel_normalize_regions)." << endl;
        //     exit(1);
        // }
        // cerr << "before adding teh extra characters" << endl;
        // for (auto hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }
        auto old_haplotypes = get<0>(haplotypes);
        get<0>(haplotypes).clear();
        for (auto old_hap : old_haplotypes)
        {
            get<0>(haplotypes).emplace("A" + old_hap);
            sequence_added_because_empty_node.first = true;
        }
        // cerr << "aftter adding teh extra characters" << endl;
        // for (auto hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }
    }
    // else if (_nodes_to_delete.find(region.first) != _nodes_to_delete.end())
    // {
    //     cerr << "earlier error of the two. " << endl;
    //     cerr << "here is the node to delete: " << region.second << endl;
    //     cerr << "here is the _nodes_to_delete.size(): " << _nodes_to_delete.size() << endl;
    //     cerr << "here is the _nodes_to_delete: " <<endl;
    //     for (auto node : _nodes_to_delete)
    //     {
    //         cerr << node << " ";
    //     }
    //     cerr << endl;
    //     cerr << "ERROR: found a node_to_delete that isn't of length zero. This shouldn't happen." << endl;
    //     exit(1);
    // }

    // check to see if the leftmost or rightmost node is empty. If so, treat the blank
    // node as containing a character "A". (this is important for dealing with how
    // get_parallel_normalize_regions sometimes adds a blank node when isolating adjacent
    // snarls for parallelization).
    if (snarl.get_sequence(snarl.get_handle(region.second)).size() == 0)
    {
        //check that this is a node we expected to find as empty.
        // if (_nodes_to_delete.find(region.second) == _nodes_to_delete.end())
        // {
        //     cerr << "ERROR: a source or sink node for a snarl is empty in an unexpected manner (i.e., it wasn't introduced into the graph during get_parallel_normalize_regions)." << endl;
        //     exit(1);
        // }
        // cerr << "before adding teh extra characters" << endl;
        // for (auto hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }

        auto old_haplotypes = get<0>(haplotypes);
        get<0>(haplotypes).clear();
        // for (auto it = get<0>(haplotypes).begin(); it != get<0>(haplotypes).end(); it++)
        // {
            
        // }
        for (auto old_hap : old_haplotypes)
        {
            get<0>(haplotypes).emplace(old_hap + "A");
            sequence_added_because_empty_node.second = true;
            // _sequence_added_because_empty_node.second = true;
        }
        // cerr << "aftter adding teh extra characters" << endl;
        // for (auto hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }

    }
    // else if (_nodes_to_delete.find(region.second) != _nodes_to_delete.end())
    // {
    //     cerr << "here is the node to delete: " << region.second << endl;
    //     cerr << "here is the _nodes_to_delete.size(): " << _nodes_to_delete.size() << endl;
    //     cerr << "here is the _nodes_to_delete: " <<endl;
    //     for (auto node : _nodes_to_delete)
    //     {
    //         cerr << node << " ";
    //     }
    //     cerr << endl;
    //     cerr << "ERROR: found a node_to_delete that isn't of length zero. This shouldn't happen." << endl;
    //     exit(1);
    // }
    // cerr << "hap new text: " << *get<0>(haplotypes).begin() << endl;


    //todo: possibly remove the duplicate storage of gbwt info in source_to_sink_gbwt_paths, by finding a way to only pass the gbwt info to the "log_gbwt_changes" function. (currently, get<0>haplotypes will also include any source-to-sink paths embedded in the graph.)
    //deep copy of gbwt_haplotypes.
    for (vector<handle_t> hap_handles : get<0>(gbwt_haplotypes))
    {
        string hap_str;
        gbwt::vector_type hap_ids;
        for (handle_t handle : hap_handles) 
        {
            hap_ids.emplace_back(_gbwt_graph.handle_to_node(handle));
            hap_str += _graph.get_sequence(_graph.get_handle(_gbwt_graph.get_id(handle), _gbwt_graph.get_is_reverse(handle)));
        }
        
        pair<gbwt::vector_type, string> hap = make_pair(hap_ids, hap_str);
        source_to_sink_gbwt_paths.emplace_back(hap);
    }
    get<1>(haplotypes) = get<1>(gbwt_haplotypes);
    get<2>(haplotypes) = get<2>(gbwt_haplotypes);

    // Get the embedded paths in the snarl from _graph, to move them to new_snarl. Any
    // embedded paths not in gbwt are aligned in the new snarl. if stop_inclusive (e.g.
    // for parallelization), then final path handle will point to the actual handle in the
    // graph that is the furthest we want for the alignment, not one beyond. (to ensure
    // that it doesn't extend to a region that might be edited by a different snarl).
    embedded_paths = sequence_finder.find_embedded_paths(stop_inclusive);


    // TODO: once haplotypes that begin/end in the middle of the snarl have been
    //       accounted for in the code, remove next chunk of code that finds 
    //       source-to-sink paths.
    // find the paths that stretch from source to sink:
    for (auto path : embedded_paths) 
    {
        if (_graph.get_id(_graph.get_handle_of_step(path.first)) == region.first &&
            _graph.get_id(_graph.get_handle_of_step(
                _graph.get_previous_step(path.second))) == region.second)  {
            // get the sequence of the source to sink path, and add it to the
            // paths to be aligned.
            string path_seq;
            step_handle_t cur_step = path.first;
            while (cur_step != path.second) {
                path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
                cur_step = _graph.get_next_step(cur_step);
            }
            if (backwards) //guaranteed false, here. (but not where code was copied from)
            {
                get<0>(haplotypes).emplace(reverse_complement(path_seq));
            }
            else {
                get<0>(haplotypes).emplace(path_seq);
            }
        }
    }
    return;
}

/**
 * Iterates over all top-level snarls in _graph, and normalizes them.
 * @param snarl_stream file stream from .snarl.pb output of vg snarls
*/
tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> SnarlNormalizer::normalize_snarls(const vector<const Snarl *>& snarl_roots) {

    // vector<pair<id_t, id_t>> normalize_regions = get_normalize_regions(snarl_roots);

    // for (auto region : normalize_regions)
    // {
    //     cerr << "original normalize regions: " << region.first << " " << region.second << endl;
    // }

    // vector<pair<id_t, id_t>> split_normalize_regions = split_sources_and_sinks(normalize_regions);

    // for (auto region : split_normalize_regions)
    // {
    //     cerr << "split normalize regions: " << region.first << " " << region.second << endl;
    // }


    // parallel_normalization(split_normalize_regions);

    // // int num_snarls_normalized = 0;
    // // int total_num_snarls_skipped = 0;
    
    // // /**
    // //  * We keep an error record to observe when snarls are skipped because they aren't 
    // //  * normalizable under current restraints. Bools:
    // //  *      0) snarl exceeds max number of threads that can be efficiently aligned,
    // //  *      1) snarl has haplotypes starting/ending in the middle,
    // //  *      2)  some handles in the snarl aren't connected by a thread,
    // //  *      3) snarl is cyclic.
    // //  * There are two additional ints for tracking the snarl size. Ints:
    // //  *      4) number of bases in the snarl before normalization
    // //  *      5) number of bases in the snarl after normalization.
    // //  * Further error records:
    // //  *      6) snarl is trivial (either one or two nodes only), so we skipped normalizing them.
    // //  *      7) snarl has handles not represented in the gbwt, and so would be dropped if normalized.
    // //  *      8) snarl alignment includes sequences too long to suitably fit into sPOA. Need to implement abPOA.
    // // */ 
    // // int error_record_size = 8;
    // // vector<int> one_snarl_error_record(error_record_size, 0);
    // // vector<int> full_error_record(error_record_size, 0);

    // // int snarl_num = 0;
    // // for (auto region : split_normalize_regions) 
    // // {
    // //     snarl_num++;
    // //     if (_debug_print)
    // //     {
    // //         // cerr << "normalizing region number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
    // //         cerr << "normalizing region number " << snarl_num << " with source at: " << region.first << " and sink at: " << region.second << endl;
    // //     }
    // //     else if (snarl_num==1 || snarl_num%10000 == 0)
    // //     {
    // //         // cerr << "normalizing region number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
    // //         cerr << "normalizing region number " << snarl_num << " with source at: " << region.first << " and sink at: " << region.second << endl;
    // //     }
    // //         one_snarl_error_record = normalize_snarl(region.first, region.second, false, snarl_num);
    // //         if (!(one_snarl_error_record[0] || one_snarl_error_record[1] ||
    // //                 one_snarl_error_record[2] || one_snarl_error_record[3] ||
    // //                 one_snarl_error_record[6] || one_snarl_error_record[7] ||
    // //                 one_snarl_error_record[8])) {
    // //             // if there are no errors, then we've successfully normalized a snarl.
    // //             num_snarls_normalized += 1;
    // //             // track the change in size of the snarl.
    // //         } else {
    // //             // else, there was an error. Track which errors caused the snarl to not
    // //             // normalize.
    // //             // note: the ints 4 and 5 are ignored here b/c they're for
    // //             // recording the changing size of snarls that are successfully normalized.
    // //             for (int i = 0; i < error_record_size; i++) {
    // //                 if ( i != 4 && i != 5)
    // //                 {
    // //                     full_error_record[i] += one_snarl_error_record[i];
    // //                 }
    // //             }
    // //             total_num_snarls_skipped += 1;
    // //         }
        

    // // }
    
    // // print_statistics(normalize_regions, num_snarls_normalized, total_num_snarls_skipped, full_error_record);
    tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> gbwt_update_items = make_tuple(_gbwt_graph, _gbwt_changelog, _gbwt);
    return gbwt_update_items;
}

void SnarlNormalizer::print_statistics(const vector<pair<id_t, id_t>>& normalize_regions, const int num_snarls_normalized, const int total_num_snarls_skipped, const vector<int>& full_error_record){
    cerr << "Finished normalization. Generating statistics..." << endl;
    int post_norm_net_snarl_size = 0;
    unordered_set<id_t> post_norm_touched_border_nodes;
    unordered_set<pair<id_t, id_t>> snarls_found_after_norm;
    for (auto region : normalize_regions)
    {
        if (_skipped_snarls.find(region) != _skipped_snarls.end()){
            continue;
        }
        snarls_found_after_norm.emplace(region);
        
        //iterate through all regions after normalization, to measure total amount of sequence within those regions.
        SubHandleGraph region_graph = extract_subgraph(_graph, region.first, region.second);

        region_graph.for_each_handle([&](const handle_t handle){

            if (region_graph.get_id(handle) == region.first || region_graph.get_id(handle) == region.second)
            {
                // if this node is a border node, only count it if it hasn't been counted already. 
                if(post_norm_touched_border_nodes.find(region_graph.get_id(handle)) == post_norm_touched_border_nodes.end())
                {
                    post_norm_touched_border_nodes.emplace(region_graph.get_id(handle));
                    post_norm_net_snarl_size += region_graph.get_sequence(handle).size();
                }
            }
            else
            {
                // if this node isn't a border node, it hasn't been counted. Add it to the count.
                post_norm_net_snarl_size += region_graph.get_sequence(handle).size();
                
            }
        });
    }

    int num_top_snarls_tracked = 10; // hardcoded for debugging convenience. //todo: change?

    // Detect the top best snarl changes:
    auto best_compare = [](const pair<pair<id_t, id_t>, pair<int, int>> left, const pair<pair<id_t, id_t>, pair<int, int>> right) 
    {
        return (right.second.second - right.second.first) >= (left.second.second - left.second.first);
    };
    vector<pair<pair<id_t, id_t>, pair<int, int>>> best_snarl_changes;


    for (pair<pair<id_t, id_t>, pair<int, int>> region : _snarl_size_changes)
    {

        if (best_snarl_changes.size() < num_top_snarls_tracked)
        {
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
        }
        else if ((region.second.second - region.second.first) < (best_snarl_changes.back().second.second - best_snarl_changes.back().second.first)) //todo: check that this is proper comparison.
        {
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
            best_snarl_changes.pop_back();
        }
    }
    int size_of_all_shrinking_snarls_pre_norm = 0;
    int size_of_all_shrinking_snarls_post_norm = 0;
    int size_of_all_growing_snarls_pre_norm = 0;
    int size_of_all_growing_snarls_post_norm = 0;
    // Detect the top worst snarl changes:
    int snarls_that_grow_after_norm = 0;
    int snarls_that_shrink_after_norm = 0;
    auto worst_compare = [](const pair<pair<id_t, id_t>, pair<int, int>> left, const pair<pair<id_t, id_t>, pair<int, int>> right) 
    {
        return (right.second.second - right.second.first) <= (left.second.second - left.second.first);
    };
    vector<pair<pair<id_t, id_t>, pair<int, int>>> worst_snarl_changes;
    for (pair<pair<id_t, id_t>, pair<int, int>> region : _snarl_size_changes)
    {
        
        if (region.second.second - region.second.first > 0)
        {
            snarls_that_grow_after_norm += 1;
            size_of_all_growing_snarls_pre_norm += region.second.first;
            size_of_all_growing_snarls_post_norm += region.second.second;
        }
        else if (region.second.second - region.second.first < 0)
        {
            snarls_that_shrink_after_norm += 1;
            size_of_all_shrinking_snarls_pre_norm += region.second.first;
            size_of_all_shrinking_snarls_post_norm += region.second.second;
        }
        // cerr << "deciding whether or not to insert following value: " << region.second.first << " , " << region.second.second << endl;
        if (worst_snarl_changes.size() < num_top_snarls_tracked)
        {
            // cerr << "inserted because the worst_snarl_changes isn't big enough." << endl;
            worst_snarl_changes.insert(upper_bound(worst_snarl_changes.begin(), worst_snarl_changes.end(), region, worst_compare), region);
        }
        else if ((region.second.second - region.second.first) > (worst_snarl_changes.back().second.second - worst_snarl_changes.back().second.first))
        {
            // cerr << "inserted and replaced an item. shrinkage of Item replaced: " << worst_snarl_changes.back().second.second - worst_snarl_changes.back().second.first << endl;
            // cerr << "Shrinkage of item added: " << region.second.second - region.second.first << endl;
            worst_snarl_changes.insert(upper_bound(worst_snarl_changes.begin(), worst_snarl_changes.end(), region, worst_compare), region);
            worst_snarl_changes.pop_back();
        }
    }



    //todo: replace error_record system with object-wide variables that are edited whenever needed in normalize_snarl(). If I feel like making the code easier to read and edit.
    // float percent_snarl_sequence_reduced = static_cast<float>((post_norm_net_snarl_size/_pre_norm_net_snarl_size)*100.0);
    float percent_snarl_sequence_change = ((static_cast<float>(post_norm_net_snarl_size)-static_cast<float>(_pre_norm_net_snarl_size))/static_cast<float>(_pre_norm_net_snarl_size))*100.0;
    // cerr.precision(2);
    cerr << endl
         << "normalization arguments:" << endl
         << "aligner (-A): " << _alignment_algorithm << endl
         << "max normalization region size (-k): " << _max_region_size << " snarls" << endl //todo: change to "bases" if/when I change _max_region_size's metric.
        //  << "max snarl spacing (-i): " << _max_snarl_spacing << endl //todo: add units. Handles? bases? I don't recall.
         << "max number of threads in an alignment (-m): " << _max_alignment_size << endl
         << "~" << endl
         << "normalized " << num_snarls_normalized << " snarls, skipped "
         << total_num_snarls_skipped << " snarls because. . .\nthey exceeded the size limit ("
         << full_error_record[0] << " snarls),\n"
         << "had haplotypes starting/ending in the middle of the snarl ("
         << full_error_record[1] << "),\n"
         << "the snarl was cyclic (" << full_error_record[3] << " snarls)" << endl
         << "or the snarl contained handles not represented by the GBWT haplotypes (" << full_error_record[7] << " snarls).\n" 
         << "fraction of snarls skipped as a result of incomplete gbwt: " << ((double)_skipped_snarls.size()/((double)_unskipped_snarl_num + (double)_skipped_snarls.size()))*100 << "%" << endl 
         << "average size of a skipped-because-gbwt snarl: " << (double)_skipped_snarl_sizes/(double)_skipped_snarls.size() << " bases" << endl
         << "average size of an unskipped snarl: " << (double)_unskipped_snarl_sizes/(double)_unskipped_snarl_num << " bases" << endl
         << "total quantity of skipped-because-gbwt sequence: " << _skipped_snarl_sizes << " bases" << endl
         << "total quantity of unskipped sequence: " << _unskipped_snarl_sizes << " bases" << endl
         << "fraction of sequence skipped as a result of incomplete gbwt: " << ((double)_skipped_snarl_sizes/((double)_unskipped_snarl_sizes + (double)_skipped_snarl_sizes))*100 << "%" << endl 

        //  << "or the snarl was trivial - composed of only one or two nodes (" //removed because trivial snarls are now removed in the clustering stage, which are not tracked..
        //  << full_error_record[6] << " snarls)."
         << endl;
    // cerr << "amount of sequence in handles that have no corresponding path in the gbwt: "
    //      << _sequence_not_touched_by_gbwt << " bases in "<< _handles_not_touched_by_gbwt << " handles." << endl;
    cerr << "amount of sequence in normalized snarls before normalization: "
         << _pre_norm_net_snarl_size << " bases" << endl;
    cerr << "amount of sequence in normalized snarls after normalization: "
         << post_norm_net_snarl_size << " bases" << endl;
    cerr << "total sequence change: "
         << post_norm_net_snarl_size - _pre_norm_net_snarl_size << " bases" << endl;
    cerr << "percent sequence change: "
         <<  percent_snarl_sequence_change << "%" << endl;
    cerr << "total snarls that shrink in size: " << snarls_that_shrink_after_norm << endl;
    cerr << "percent change of snarls that shrink in size: " << (((double)size_of_all_shrinking_snarls_post_norm - (double)size_of_all_shrinking_snarls_pre_norm)/(double)snarls_that_shrink_after_norm)*100 << "%" << endl;  
    cerr << "top " << num_top_snarls_tracked << " snarls *reduction* in size (probably desirable): " << endl;
    cerr << "leftmost_id  rightmost_id  original_size  normalized_size  (change_in_size;  percent change)" << endl;
    
    // for (auto i = _snarl_size_changes.end(); i != next(_snarl_size_changes.end(), -num_top_snarls_tracked); i--)
    for (auto region : best_snarl_changes)
    {
        cerr << region.first.first << "   " << region.first.second << "   " << region.second.first << "   " << region.second.second << "   (" << region.second.second - region.second.first << ";" << "   " << ((((double)region.second.second - (double)region.second.first)/(double)region.second.first)*100.0) << "%)" << endl;
    }

    cerr << "total snarls that grow in size: " << snarls_that_grow_after_norm << endl;
    cerr << "percent change of snarls that grow in size: " << (((double)size_of_all_growing_snarls_post_norm - (double)size_of_all_growing_snarls_pre_norm)/(double)snarls_that_grow_after_norm)*100 << "%" << endl;  
    cerr << "top " << num_top_snarls_tracked << " snarls *increase* in size (probably undesirable): " << endl;
    cerr << "leftmost_id  rightmost_id  original_size  normalized_size  (change_in_size;  percent change)" << endl;
    
    // for (auto i = _snarl_size_changes.end(); i != next(_snarl_size_changes.end(), -num_top_snarls_tracked); i--)
    for (auto region : worst_snarl_changes)
    {
        cerr << region.first.first << "   " << region.first.second << "   " << region.second.first << "   " << region.second.second << "   (" << region.second.second - region.second.first << ";" << "   " << ((((double)region.second.second - (double)region.second.first)/(double)region.second.first)*100.0) << "%)" << endl;
    }

    cerr << "number of snarls calling for abpoa: " << _alignments_calling_for_abpoa.size() << endl;

}

////////////////////////////////////////////////////////////
/// Normalization of a single snarl:
////////////////////////////////////////////////////////////


/**
 * Normalize a single snarl defined by a source and sink. Only extracts and realigns 
 * sequences found in the gbwt. 
 * @param source_id the source of the snarl of interest.
 * @param sink_id the sink of the snarl of interest.
 * @param error_record an empty vector of 6 integers.
*/
// Returns: none.
// TODO: allow for snarls that have haplotypes that begin or end in the middle of the
// snarl.
vector<int> SnarlNormalizer::normalize_snarl(const id_t source_id, const id_t sink_id, const bool backwards, const int snarl_num) {
    id_t leftmost_id;
    id_t rightmost_id;
    if (backwards) {
        leftmost_id = sink_id;
        rightmost_id = source_id;
    }
    else {
        leftmost_id = source_id;
        rightmost_id = sink_id;
    }

    _snarl_size_changes[make_pair(leftmost_id, rightmost_id)] = make_pair(0, 0);
    /**
     * We keep an error record to observe when snarls are skipped because they aren't 
     * normalizable under current restraints. Bools:
     *      0) snarl exceeds max number of threads that can be efficiently aligned,
     *      1) snarl has haplotypes starting/ending in the middle,
     *      2)  some handles in the snarl aren't connected by a thread,
     *      3) snarl is cyclic.
     * There are two additional ints for tracking the snarl size. Ints:
     *      4) number of bases in the snarl before normalization
     *      5) number of bases in the snarl after normalization.
     *      6) snarl is trivial (either one or two nodes only), so we skipped normalizing them.
     *      7) snarl has handles not represented in the gbwt, and so would be dropped if normalized.
     *      8) snarl alignment includes sequences too long to suitably fit into sPOA. Need to implement abPOA.
    */ 
    vector<int> error_record(9, 0);
    SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);

    int snarl_size = 0;
    snarl.for_each_handle([&](handle_t handle){
        snarl_size += snarl.get_sequence(handle).size();
    });

    // bool debug_still_looking = true;
    
    snarl.for_each_handle([&](handle_t handle){
        const gbwt::BidirectionalState debug_state = _gbwt_graph.get_bd_state(_gbwt_graph.get_handle(snarl.get_id(handle)));
        if (debug_state.empty())
        {
            // cerr << "there is an empty debug_state at: " << snarl.get_id(handle) << endl;
            // if (snarl_size <= 20)
            // {
            //     cerr << "skippable snarl of small size " << snarl_size << " is at leftmost: " << leftmost_id << " and rightmost: " << rightmost_id << endl; 
            //     debug_still_looking = false;
            // }
            // _snarls_skipped_because_gbwt_misses_handles++;
            error_record[7] = true;
            return;
        }
        // else
        // {
        //     return false;
        // }
    }, true);

    if (error_record[7])
    {
        _skipped_snarl_sizes += snarl_size;
        _skipped_snarls.emplace(make_pair(leftmost_id, rightmost_id));
        return error_record;
    }

    //todo: debug_statement: Evaluate connections of all nodes in subgraph.
    // snarl.for_each_handle([&](const handle_t handle){
    //     cerr << "examining left neighbors of handle " << snarl.get_id(handle) << ":" << endl;
    //     snarl.follow_edges(handle, false, [&](const handle_t next) {
    //         cerr << "     " << snarl.get_id(next) << " ";
    //     });
    //     cerr << endl;
    // });

    if (!handlealgs::is_acyclic(&snarl)) {
        cerr << "snarl at " << source_id << " is cyclic. Skipping." << endl;
        error_record[3] = true;
        _skipped_snarl_sizes += snarl_size;
        _skipped_snarls.emplace(make_pair(leftmost_id, rightmost_id));
        return error_record;
    }

    // only normalize non-trivial snarls (i.e. not composed of just a source and sink.):
    // int num_handles_in_snarl = 0;
    // snarl.for_each_handle([&](const handle_t handle){
    //     num_handles_in_snarl++;
        // if (num_handles_in_snarl >= 3)
        // {
        //     return;
        // }
    // });
    // if (num_handles_in_snarl <= _max_region_size+1) //note: this calculation of trivial snarls is incorrect. This is because the batch size is not a correct indicator of trivial-snarl-batch size. Sometimes, a snarl cluster is smaller than _max_region_size, because it was cut off early (because there was too large a gap between snarls, for example). All undesired trivial snarls should have been removed in the clustering stage.
    // {
    //     // cerr << "trivial, so skipping." << endl;
    //     // if (_debug_print)
    //     // {
    //     //     cerr << "snarl with source " << source_id << " and sink " << sink_id << " has"
    //     //         << " only " << num_handles_in_snarl << " nodes. Skipping normalization of"
    //     //         << " trivial snarl." << endl;
    //     // }
    //     error_record[6] += 1;
    //     return error_record;
    // }

    // cerr << "num_handles_in_snarl: " << num_handles_in_snarl << endl;

    // check to make sure that the gbwt _graph has threads connecting all handles:
    // ( needs the unordered_set from extract_gbwt haplotypes to be equal to the number of
    // handles in the snarl).
    unordered_set<id_t> nodes_in_snarl;
    snarl.for_each_handle([&](const handle_t handle) {
        nodes_in_snarl.emplace(snarl.get_id(handle));
        // count the number of bases in the snarl (buggy counter that includes doublecounting of border nodes between snarls in a chain).
        error_record[4] += snarl.get_sequence(handle).size();

        // count the number of bases in the snarl (fixed counter, no doublecounting).
        if(snarl.get_id(handle) == leftmost_id || snarl.get_id(handle) == rightmost_id)
        {
            // if this node is a border node, only count it if it hasn't been counted already. 
            if(_touched_border_nodes.find(snarl.get_id(handle)) == _touched_border_nodes.end())
            {
                _pre_norm_net_snarl_size+=snarl.get_sequence(handle).size();
                _touched_border_nodes.emplace(snarl.get_id(handle));
            }
            _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].first += snarl.get_sequence(handle).size(); 
        }
        else
        {
            // if this node isn't a border node, it hasn't been counted. Add it to the count.
            _pre_norm_net_snarl_size+=snarl.get_sequence(handle).size();
            _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].first += snarl.get_sequence(handle).size(); 

        }
    });
    // initialize the post-normalization snarl size to be the same as pre-norm. This is so
    // that, if the program drops this snarl without normalization, the snarl is properly 
    // registered as unchanged in size.
    _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].second = _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].first;

    // extract threads
    // haplotypes is of format:
    // 0: a set of all the haplotypes which stretch from source to sink, in string format.
    //   - it's a set, so doesn't contain duplicates
    // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
    // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
    tuple<unordered_set<string>, vector<vector<handle_t>>, unordered_set<id_t>> haplotypes;
    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, source_id, sink_id, backwards);
    vector<pair<gbwt::vector_type, string>> source_to_sink_gbwt_paths;
    if (_path_finder == "GBWT") {
        tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<id_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();


        // cerr << "sizes of gbwt_haplotypes output: " << get<0>(gbwt_haplotypes).size() << " " << get<1>(gbwt_haplotypes).size() << " " << get<2>(gbwt_haplotypes).size() << endl;

        // cerr << "each gbwt hap output: " << endl;
        // for (vector<handle_t> hap : get<0>(gbwt_haplotypes))
        // {
        //     for (handle_t handle : hap)
        //     {
        //         cerr << _gbwt_graph.get_id(handle) << " ";
        //     }
        //     cerr << "(";
        //     for (handle_t handle : hap)
        //     {
        //         cerr << _gbwt_graph.get_sequence(handle);
        //     }
        //     cerr << ")" << endl;
        // }
            

        // cerr << "check that all handles touched by find_gbwt_haps are all the handles in the subgraph:" << endl;
        int debug_sequence_not_in_gbwt = 0;
        for(id_t node_id : nodes_in_snarl)
        {
            if (get<2>(gbwt_haplotypes).find(node_id) == nodes_in_snarl.end()){
                const gbwt::BidirectionalState debug_state = _gbwt_graph.get_bd_state(_gbwt_graph.get_handle(node_id));
                if (debug_state.empty())
                {
                    cerr << "since this situation is handled in a previous check, you should never see this message." << endl;
                    exit(1);

                    // debug_sequence_not_in_gbwt += _graph.get_sequence(_graph.get_handle(node_id)).size();
                    // _handles_not_touched_by_gbwt++;
                    // _sequence_not_touched_by_gbwt += _graph.get_sequence(_graph.get_handle(node_id)).size();
                }
                else
                {
                    cerr << "ERROR: the node " << node_id << " in the graph at snarl " << snarl_num << " with source " << source_id << " and sink " << sink_id << " is not touched by find_gbwt_haps(), but the node still exists in the gbwt. This is not supposed to happen. Exiting the program." << endl;
                    exit(1);
                }
                // cerr << "WARNING: the node " << node_id << " in the graph at snarl " << snarl_num << " with source " << source_id << " and sink " << sink_id << " that are not touched by the GBWT. The sequence information in these handles will be dropped from the normalized graph." << endl;
            }
        }
        // cerr << "**** total sequence that was removed from snarl " << snarl_num << " with source " << source_id << " and sink " << sink_id << " is: " << debug_sequence_not_in_gbwt << endl;
        // cerr << "sizes of fields in gbwt_haplotypes: " << get<0>(gbwt_haplotypes).size() << " " << get<1>(gbwt_haplotypes).size() << " " << get<2>(gbwt_haplotypes).size() << endl;
        // //todo: comment out debug
        // for (id_t node_id : get<2>(gbwt_haplotypes))
        // {
        //     if (node_id == 2605470)
        //     {
        //         cerr << "while iterating through the touched handles, " << endl;
        //         cerr << "found the handle missing from the updated gbwt. It's touched in snarl number " << snarl_num << " with leftmost_id " << leftmost_id << " and rightmost_id " << rightmost_id << endl;
        //     }
        // }
        // snarl.for_each_handle([&](const handle_t handle) {
        //     if (snarl.get_id(handle) == 2605470)
        //     {
        //         cerr << "while iterating through all handles in the snarl, " << endl;
        //         cerr << "found the handle missing from the updated gbwt. It's touched in snarl number " << snarl_num << " with leftmost_id " << leftmost_id << " and rightmost_id " << rightmost_id << endl;
        //     }

        // });


        // cerr << "check that all handles touched by find_gbwt_haps are all the handles in the subgraph:" << endl;
        // for (handle_t handle : get<2>(gbwt_haplotypes))
        // {
        //     cerr << _graph.get_id(handle) << endl;
        // }
        // cerr << endl;

        // cerr << "various sizes: " << get<0>(gbwt_haplotypes).size() << " " << get<1>(gbwt_haplotypes).size() << " " << get<2>(gbwt_haplotypes).size() << endl;

        // cerr << "naive? gbwt haplotypes extract: " << endl;
        // for (auto hap : get<0>(gbwt_haplotypes)) 
        // {
        //     cerr << "new hap" << endl;
        //     for (auto handle : hap)
        //     {
        //         cerr << "handle id: " << _graph.get_id(handle) << " seq: " << _graph.get_sequence(handle) << endl;
        //     }
        // }
        // Convert the haplotypes from vector<handle_t> format to string format.
        get<0>(haplotypes) = format_handle_haplotypes_to_strings(_graph, _gbwt_graph, get<0>(gbwt_haplotypes));

        //todo: possibly remove the duplicate storage of gbwt info in source_to_sink_gbwt_paths, by finding a way to only pass the gbwt info to the "log_gbwt_changes" function. (currently, get<0>haplotypes will also include any source-to-sink paths embedded in the graph.)
        //deep copy of gbwt_haplotypes.
        for (vector<handle_t> hap_handles : get<0>(gbwt_haplotypes))
        {
            string hap_str;
            gbwt::vector_type hap_ids;
            for (handle_t handle : hap_handles) 
            {
                // if (_gbwt_graph.get_id(handle) == 7405162)
                // {
                //     cerr << "id of node: " << _gbwt_graph.get_id(handle) << endl;
                //     cerr << "sequence of node: " << _gbwt_graph.get_sequence(handle) << endl;
                // }
                // cerr << "id of node: " << _gbwt_graph.get_id(handle) << endl;
                // cerr << "sequence of node: " << _gbwt_graph.get_sequence(handle) << endl;
                hap_ids.emplace_back(_gbwt_graph.handle_to_node(handle));
                // hap_str += _gbwt_graph.get_sequence(handle);
                hap_str += _graph.get_sequence(_graph.get_handle(_gbwt_graph.get_id(handle), _gbwt_graph.get_is_reverse(handle)));
                // cerr << "hap_str: " << hap_str << endl;
            }
            
            pair<gbwt::vector_type, string> hap = make_pair(hap_ids, hap_str);
            source_to_sink_gbwt_paths.emplace_back(hap);
        }
        // for (auto item : source_to_sink_gbwt_paths)
        // {
        //     for (auto nid : item.first)
        //     {
        //         cerr << "is the handle is-reverse? of the handles in 'before': " << _graph.get_is_reverse(_graph.get_handle(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(nid)), _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(nid)))) << endl;

        //     }

        // }
        get<1>(haplotypes) = get<1>(gbwt_haplotypes);
        get<2>(haplotypes) = get<2>(gbwt_haplotypes);
        // cerr << "haplotypes after formatting to strings: " << endl;
        // for (auto hap : get<0>(haplotypes)) 
        // {
        //     cerr << "hap: " << hap << endl;
        // }
        
    } else if (_path_finder == "exhaustive") {
        //todo: to enable support for exhaustive, make tests, run them, and also set up support for when I log changes for the gbwt update.
        cerr << "'exhaustive' path finder currently unsupported. Use 'GBWT'. '" << "'." << endl;
        exit(1);

        // pair<unordered_set<string>, unordered_set<handle_t>> exhaustive_haplotypes =
        //     sequence_finder.find_exhaustive_paths();
        // get<0>(haplotypes) = exhaustive_haplotypes.first;
        // get<2>(haplotypes) = exhaustive_haplotypes.second;
    } else {
        cerr << "path_finder type must be 'GBWT' or 'exhaustive', not '" << _path_finder
             << "'." << endl;
        exit(1);
    }




    // Print a heads-up about snarls that require an alignment with a greater number of 
    // threads than _big_snarl_alignment_job, so the user knows if they are hung up on a 
    // large job.
    if (get<0>(haplotypes).size() > _big_snarl_alignment_job)
    {
        cerr << "WARNING: aligning a snarl requiring a large number (>" << _big_snarl_alignment_job <<") of threads. Number of threads: " << get<0>(haplotypes).size() << endl;
    }
    // Record start time, for measuring alignment time for a big snarls:
    //todo: also add compute time (vs wall clock) measure.
    auto _big_snarl_time_start = chrono::high_resolution_clock::now();    

    // TODO: this if statement only permits snarls that satsify requirements, i.e.
    // TODO:    there are no haplotype begins/ends in the middle
    // TODO:    of the snarl. Get rid of this once alignment issue is addressed!
    // TODO: also, limits the number of haplotypes to be aligned, since snarl starting at
    // TODO:    2049699 with 258 haplotypes is taking many minutes.
    if (get<1>(haplotypes).empty() && get<0>(haplotypes).size() <= _max_alignment_size)
        // the following bool check was to ensure that all the handles in the handlegraph 
        // are touched by the gbwt. Turns out, though, that this isn't necessary. If we 
        // assume that all seq info is in gbwt, the gbwt is all we need to worry about.:
        // && get<2>(haplotypes).size() == handles_in_snarl.size()) {
        {
        // Get the embedded paths in the snarl from _graph, to move them to new_snarl.
        // Any embedded paths not in gbwt are aligned in the new snarl.
        vector<pair<step_handle_t, step_handle_t>> embedded_paths =
            sequence_finder.find_embedded_paths();

        //todo: debug_statement
        // cerr << "strings in path_seq before adding haplotypes: " << endl;
        // for (auto path : get<0>(haplotypes))
        // {
        //     cerr << path << endl;
        // }

        
        // TODO: once haplotypes that begin/end in the middle of the snarl have been
        // TODO:    accounted for in the code, remove next chunk of code that finds 
        // TODO: source-to-sink paths.
        // find the paths that stretch from source to sink:
        // cerr << "~~~~~~~~~~source: " << source_id << "sink: " << sink_id << endl;
        for (auto path : embedded_paths) 
        {

            // cerr << "checking path of name " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << " with source " << _graph.get_id(_graph.get_handle_of_step(path.first)) << " and sink " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) <<  endl;
            // cerr << "SOURCE info: prev step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << "prev prev step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(_graph.get_previous_step(path.second)))) << " source: " << _graph.get_id(_graph.get_handle_of_step(path.second)) << " next step: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_next_step(path.second))) << endl;
            // cerr << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << " " << source_id << " source bool: " <<  (_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == source_id) << endl;
            if (_graph.get_id(_graph.get_handle_of_step(path.first)) == source_id &&
                _graph.get_id(_graph.get_handle_of_step(
                    _graph.get_previous_step(path.second))) == sink_id)  {
                // cerr << "path_seq added to haplotypes. " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << endl;

                // cerr << "******************************************\nadding path of name " <<
                // _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) <<
                // endl; 
                // get the sequence of the source to sink path, and add it to the
                // paths to be aligned.
                string path_seq;
                step_handle_t cur_step = path.first;
                while (cur_step != path.second) {
                    // cerr << "while adding path, looking at node " << _graph.get_id(_graph.get_handle_of_step(cur_step)) << " with seq " << _graph.get_sequence(_graph.get_handle_of_step(cur_step)) << endl;
                    path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
                    cur_step = _graph.get_next_step(cur_step);
                }
                // cerr << "path seq:" << path_seq << endl;
                if (backwards) {
                    // cerr << "path seq emplaced (in reverse):" << reverse_complement(path_seq)  << endl;
                    // int init_hap_size = get<0>(haplotypes).size(); // Note: just for debug purposes.
                    get<0>(haplotypes).emplace(reverse_complement(path_seq));
                    // cerr << "was path_seq a new string? " << get<0>(haplotypes).size() - init_hap_size << endl;
                }
                else {
                    // cerr << "path seq emplaced (in forward):" << path_seq  << endl;
                    // int init_hap_size = get<0>(haplotypes).size(); // Note: just for debug purposes.
                    get<0>(haplotypes).emplace(path_seq);
                    // cerr << "was path_seq a copy? " << get<0>(haplotypes).size() - init_hap_size << endl;

                }
            }
        }
        
        int max_spoa_length = 750; // somewhere between 500-1000 bases, sPOA starts to struggle. That's why I'll eventually want abPOA to take over.
        for (string hap : get<0>(haplotypes))
        {
            if (hap.size() > max_spoa_length)
            {
                _skipped_snarls.emplace(make_pair(leftmost_id, rightmost_id));
                _alignments_calling_for_abpoa.push_back(snarl_num);
                return error_record;
            }
        }
        
        // cerr << "haps in haplotypes: " << endl;
        // for (string hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }
        // Align the new snarl:
        shared_ptr<MutableHandleGraph> new_snarl;
        // if (_alignment_algorithm == "TCoffee")
        // {
        //     new_snarl = align_source_to_sink_haplotypes(get<0>(haplotypes));
        // }
        // else if (_alignment_algorithm == "sPOA")
        if (_alignment_algorithm == "sPOA")
        {
            new_snarl = poa_source_to_sink_haplotypes(get<0>(haplotypes), false);
            if (!new_snarl)
            {
                //note: this snippet probably never needs to run. It's also handled by the "if (hap.size() > max_spoa_length)" condition a few lines above.
                _skipped_snarls.emplace(make_pair(leftmost_id, rightmost_id));
                _alignments_calling_for_abpoa.push_back(snarl_num);
                return error_record;
            }
            // if (leftmost_id == 996838)
            // {
            //     new_snarl = poa_source_to_sink_haplotypes(get<0>(haplotypes), snarl_num, true);
            // }
            // else
            // {
            //     new_snarl = poa_source_to_sink_haplotypes(get<0>(haplotypes), snarl_num);
            // }
        }
        // else if (_alignment_algorithm == "kalign") //todo: implement use of kalign. Then, update the error message in the else statement.
        // {

        // }
        else
        {
            cerr << "error:[vg normalize] _alignment_algorithm variable must be set as either T-Coffee or sPOA." << endl;
            exit(1);
        }

        //preprocess new_snarl for log_gbwt_changes:
        bool single_stranded = handlealgs::is_single_stranded(&(*new_snarl));
        shared_ptr<MutableHandleGraph> single_stranded_snarl;
        if (!single_stranded) 
        {
            handlealgs::split_strands(&(*new_snarl), &(*single_stranded_snarl));
            // handlealgs::SplitStrandOverlay(new_snarl)
        }
        else
        {
            single_stranded_snarl=new_snarl;
        }

        //todo: skipping dagification because I require the input snarl to be a DAG, and I don't think alignments of sequences should produce non-DAGs.
        // bool dag = handlealgs::is_directed_acyclic(&new_snarl);
        // if (!dag)
        // {
        //     handlealgs::dagify(single_stranded_snarl, dagified_snarl, );
        // }
        // else
        // {

        // }

        // count the number of bases in the snarl.
        // (and reinitialize the post-normalization at 0, so that we can properly count 
        // postnormalized snarl size.)
        _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].second = 0;
        new_snarl->for_each_handle([&](const handle_t handle) {
            error_record[5] += new_snarl->get_sequence(handle).size();
            _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].second += new_snarl->get_sequence(handle).size(); //todo: robin-review-and-fix this so that I can have correct change by initializing this here as zero, otherwise set it as the original size of the snarl when setting pair.first value.. 
        });
        // force_maximum_handle_size(*new_snarl); //Note: the maximum handle size is now enforced by the MSAConverter
        
        // integrate the new_snarl into the _graph, removing the old snarl as you go.
        // //todo: debug_statement
        // integrate_snarl(new_snarl, embedded_paths, sink_id, source_id);
        pair<handle_t, handle_t> new_left_right = integrate_snarl(snarl, *new_snarl, embedded_paths, source_id, sink_id, backwards);
        _unskipped_snarls.emplace(make_pair(leftmost_id, rightmost_id));

        // make a subhandlegraph of the normalized snarl to find the new gbwt paths in the graph.
        SubHandleGraph integrated_snarl = extract_subgraph(_graph, _graph.get_id(new_left_right.first), _graph.get_id(new_left_right.second));

        log_gbwt_changes(source_to_sink_gbwt_paths, integrated_snarl);

        // integrated_snarl.for_each_handle([&](const handle_t handle) {
        //     if (integrated_snarl.get_id(handle) == 2605470)
        //     {
        //         cerr << "while iterating through all handles in the snarl, " << endl;
        //         cerr << "found the handle missing from the updated gbwt. It's touched in snarl number " << snarl_num << " with leftmost_id " << leftmost_id << " and rightmost_id " << rightmost_id << endl;
        //     }

        // });


        // Print a heads-up about snarls that require an alignment with a greater number of 
        // threads than _big_snarl_alignment_job, so the user knows if they are hung up on a 
        // large job.
        if (get<0>(haplotypes).size() > _big_snarl_alignment_job)
        {
            // Record end time
            auto _big_snarl_time_finish = std::chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = _big_snarl_time_finish - _big_snarl_time_start;
            cerr << "big snarl with " << get<0>(haplotypes).size() << " threads for alignment finished normalization." << endl;
            cerr << "Elapsed time normalizing snarl after sequence extraction: " << elapsed.count() << " s\n";
        }
    } else {
        if (!get<1>(haplotypes).empty()) {
            cerr << "found a snarl with source " << source_id << " and sink "
                 << sink_id
                 << " with haplotypes that start or end in the middle. Skipping." << endl;
            cerr << "There are " << sizeof(get<1>(haplotypes)) << " haplotypes of that description." << endl;
            // vector<string> string_haps = format_handle_haplotypes_to_strings(get<1>(haplotypes).front());
            // cerr << "First example: " << get<1>(haplotypes) << endl;
            error_record[1] = true;
        }
        if (get<0>(haplotypes).size() > _max_alignment_size) {
            cerr << "found a snarl with source " << source_id << " and sink "
                 << sink_id << " with too many haplotypes (" << get<0>(haplotypes).size()
                 << "). Is greater than max, " << _max_alignment_size <<" Skipping." << endl;
            error_record[0] = true;
        }
        // if (get<2>(haplotypes).size() != handles_in_snarl.size()) {
        //     cerr << "some handles in the snarl with source " << source_id
        //          << " and sink " << sink_id
        //          << " aren't accounted for by the gbwt_graph. "
        //             "Skipping."
        //          << endl;
        //     cerr << "handles in snarl:" << handles_in_snarl.size() << "number of handles touched by gbwt graph: " << get<2>(haplotypes).size() << endl;
        //     cerr << "these handles are:" << endl << "\t";
        //     for (auto handle : handles_in_snarl) {
        //         if (get<2>(haplotypes).find(handle) == get<2>(haplotypes).end()) {
        //             cerr << _graph.get_id(handle) << " ";
        //         }
        //     }
        //     cerr << endl;
        //     error_record[2] = true;
        // }
    }
    // todo: decide if we should only normalize snarls that decrease in size.
    if (error_record[5] > error_record[4]) {
        if (_debug_print)
        {
            cerr << "**************************in UNIT-TEST for normalize_snarl: **************************" << endl;
            cerr << "NOTE: normalized a snarl which *increased* in sequence quantity, "
                    "with source: " << source_id << " and sink: " << sink_id << endl
                << "\tsize before: " << error_record[4] << " size after: " << error_record[5]
                << endl;
        }
    }
    else if (!error_record[0] && error_record[5] <= 0 && !error_record[1]) {
        cerr << "normalized snarl size is <= zero: " << error_record[5] << endl;
        cerr << "snarl number: " << snarl_num << endl;
    }
    _unskipped_snarl_sizes+=snarl_size; //todo: remove for increased efficiency? Or at least turn into a rolling calculation of averages. (just a rolling sum + a tracker of total number skipped). 
    _unskipped_snarl_num++;

    return error_record;

}


// Given a vector of haplotypes of format vector< handle_t >, returns a vector of
// haplotypes of
//      format string (which is the concatenated sequences in the handles).
// Arguments:
//      > haplotypes. haplotypte_handle_vectors: a vector of haplotypes in vector<
//      handle_t > format. the handles are from the _gbwt_graph.
// Returns: a vector of haplotypes of format string (which is the concatenated sequences
// in the handles).
unordered_set<string> SnarlNormalizer::format_handle_haplotypes_to_strings(const HandleGraph& graph, const gbwtgraph::GBWTGraph & gbwt_graph,
    const vector<vector<handle_t>> &haplotype_handle_vectors) {
    unordered_set<string> haplotype_strings;
    for (vector<handle_t> haplotype_handles : haplotype_handle_vectors) {
        string hap;
        for (handle_t handle : haplotype_handles) {
            // hap += _gbwt_graph.get_sequence(handle);
            hap += graph.get_sequence(graph.get_handle(gbwt_graph.get_id(handle), gbwt_graph.get_is_reverse(handle)));
        }
        haplotype_strings.emplace(hap);
    }
    return haplotype_strings;
}

// TODO: eventually change to deal with haplotypes that start/end in middle of snarl.
// Aligns haplotypes to create a new _graph using MSAConverter's seqan converter.
//      Assumes that each haplotype stretches from source to sink.
// Arguments:
//      source_to_sink_haplotypes: a vector of haplotypes in string format (concat of
//      handle sequences).
// Returns:
//      VG object representing the newly realigned snarl.
unique_ptr<MutablePathDeletableHandleGraph> SnarlNormalizer::align_source_to_sink_haplotypes(
    const unordered_set<string>& source_to_sink_haplotypes) {
    // cerr << "align_source_to_sink_haplotypes" << endl;
    // cerr << " haplotypes in source_to_sink_haplotypes: " << endl;
    // for (string hap : source_to_sink_haplotypes) {
    //     cerr << hap << endl;
    // }
    // cerr << "number of strings to align: " << source_to_sink_haplotypes.size() << endl;
    // TODO: make the following comment true, so that I can normalize haplotypes that
    // TODO:    aren't source_to_sink by adding a similar special character to strings in
    // TODO:    the middle of the snarl.
    // modify source_to_sink_haplotypes to replace the leading and
    // trailing character with a special character. This ensures that the leading char of
    // the haplotype becomes the first character in the newly aligned snarl's source - it
    // maintains the context of the snarl.

    // store the source/sink chars for later reattachment to source and sink.
    string random_element;
    for (auto hap : source_to_sink_haplotypes){
        random_element = hap;
        break;
    }
    string source_char(1, random_element.front());
    string sink_char(1, random_element.back());

    // cerr << "strings in path_seq before replacing final character: " << endl;
    // for (auto path : source_to_sink_haplotypes)
    // {
    //     cerr << path << endl;
    // }

    // replace the source and sink chars with X, to force match at source and sink.
    unordered_set<string> edited_source_to_sink_haplotypes;
    // for (auto it = source_to_sink_haplotypes.begin(); it != source_to_sink_haplotypes.end(); it++)
    for (auto hap : source_to_sink_haplotypes)
    {
        // cerr << "hap before replace: " << hap << endl;
        hap.replace(0, 1, "X");
        hap.replace(hap.size() - 1, 1, "X");
        // cerr << "hap after replace: " << hap << endl;
        edited_source_to_sink_haplotypes.emplace(hap);
    }
    // cerr << "source_char: " << source_char << endl;
    // cerr << "sink_char: " << sink_char << endl;

    // //todo: debug_statement
    // source_to_sink_haplotypes.emplace_back("XX");

    // /// make a new scoring matrix with _match=5, _mismatch = -3, _gap_extend = -1, and
    // _gap_open = -3, EXCEPT that Q has to be matched with Q (so match score between Q
    // and Q =len(seq)+1)
    // // 1. Define type and constants.
    // //
    // // Define types for the score value and the scoring scheme.
    // typedef int TValue;
    // typedef seqan::Score<TValue, seqan::ScoreMatrix<seqan::Dna5, seqan::Default> >
    // TScoringScheme;
    // // Define our gap scores in some constants.
    // int const gapOpenScore = -1;
    // int const gapExtendScore = -1;

    // static int const _data[TAB_SIZE] =
    //     {
    //         1, 0, 0, 0, 0,
    //         0, 1, 0, 0, 0,
    //         0, 0, 1, 0, 0,
    //         0, 0, 0, 1, 0,
    //         0, 0, 0, 0, 0
    //     };

    // create seqan multiple_sequence_alignment object
    //// seqan::Align<seqan::DnaString>   align;
    seqan::Align<seqan::CharString> align;

    seqan::resize(rows(align), edited_source_to_sink_haplotypes.size());
    int i = 0;
    for (auto hap : edited_source_to_sink_haplotypes) {
        assignSource(row(align, i), hap.c_str());
        i++;
    }

    globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));

    vector<string> row_strings;
    for (auto &row : rows(align)) {
        string row_string;
        auto it = begin(row);
        auto itEnd = end(row);
        for (; it != itEnd; it++) {
            row_string += *it;
        }
        // todo: debug_statement
        // cerr << "ROW_STRING: " << row_string << endl;
        // edit the row so that the proper source and sink chars are added to the
        // haplotype instead of the special characters added to ensure correct alignment
        // of source and sink.
        // cerr << "row_string before: " << row_string << endl;
        row_string.replace(0, 1, source_char);
        row_string.replace(row_string.size() - 1, 1, sink_char);
        row_strings.push_back(row_string);
        // cerr << "row_string after: " << row_string << endl;
    }

    stringstream ss;
    for (string seq : row_strings) {
        // todo: debug_statement
        // cerr << "seq in alignment:" << seq << endl;
        ss << endl << seq;
    }
    // ss << align;
    MSAConverter myMSAConverter = MSAConverter();
    myMSAConverter.load_alignments(ss, "seqan");
    unique_ptr<MutablePathDeletableHandleGraph> snarl = myMSAConverter.make_graph();

    // snarl.clear_paths();

    pair<vector<handle_t>, vector<handle_t>> source_and_sink =
        debug_get_sources_and_sinks(*snarl);


    // TODO: throw exception(?) instead of cerr, or remove these messages if I'm confident
    // TODO:    code works.
    if (source_and_sink.first.size() != 1) {
        cerr << "WARNING! Snarl realignment has generated "
             << source_and_sink.first.size() << " source nodes." << endl;
    }

    if (source_and_sink.second.size() != 1) {
        cerr << "WARNING! Snarl realignment has generated "
             << source_and_sink.second.size() << " sink nodes." << endl;
    }

    return snarl;
}

/** For each handle in a given _graph, divides any handles greater than max_size into
 * parts that are equal to or less than the size of max_size.
 *
 * @param  {MutableHandleGraph} _graph : the _graph in which we want to force a maximum
 * handle size for all handles.
 * @param  {size_t} max_size          : the maximum size we want a handle to be.
 */
void SnarlNormalizer::force_maximum_handle_size(MutableHandleGraph &graph) {
    // forcing each handle in the _graph to have a maximum sequence length of max_size:
    graph.for_each_handle([&](handle_t handle) {
        // all the positions we want to make in the handle are in offsets.
        vector<size_t> offsets;

        size_t sequence_len = graph.get_sequence(handle).size();
        int number_of_divisions = floor(sequence_len / _max_handle_size);

        // if the handle divides evenly into subhandles of size _max_handle_size, we don't need to
        // make the last cut (which would be at the very end of the handle - cutting off
        // no sequence).
        if (sequence_len % _max_handle_size == 0) {
            number_of_divisions--;
        }

        // calculate the position of all the divisions we want to make.
        for (int i = 1; i <= number_of_divisions; i++) {
            offsets.push_back(i * _max_handle_size);
        }

        // divide the handle into parts.
        graph.divide_handle(handle, offsets);
    });
}

// Given a start and end node id, construct an extract subgraph between the two nodes
// (inclusive). Arguments:
//      graph: a pathhandlegraph containing the snarl with embedded paths.
//      source_id: the source of the snarl of interest.
//      sink_id: the sink of the snarl of interest.
// Returns:
//      a SubHandleGraph containing only the handles in _graph that are between start_id
//      and sink_id.
SubHandleGraph SnarlNormalizer::extract_subgraph(const HandleGraph &graph,
                                                 const id_t leftmost_id,
                                                 const id_t rightmost_id) {
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    // this bool tracks when we find the rightmost_id while extending the snarl. If we 
    // never find it, then we must have been passed leftmost_id and rightmost_id in 
    // reverse order. In that case, we'll throw an error.
    bool found_rightmost = false;

    unordered_set<id_t> visited;  // to avoid counting the same node twice.
    unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

    // initialize with leftmost_handle (because we move only to the right of leftmost_handle):
    handle_t leftmost_handle = graph.get_handle(leftmost_id);
    subgraph.add_handle(leftmost_handle);
    visited.insert(graph.get_id(leftmost_handle));

    // look only to the right of leftmost_handle
    graph.follow_edges(leftmost_handle, false, [&](const handle_t handle) {
        // mark the nodes to come as to_visit
        if (visited.find(graph.get_id(handle)) == visited.end()) {
            to_visit.insert(graph.get_id(handle));
        }
    });

    /// explore the rest of the snarl:
    while (to_visit.size() != 0) {
        // remove cur_handle from to_visit
        unordered_set<id_t>::iterator cur_index = to_visit.begin();
        handle_t cur_handle = graph.get_handle(*cur_index);

        to_visit.erase(cur_index);

        /// visit cur_handle
        visited.insert(graph.get_id(cur_handle));

        subgraph.add_handle(cur_handle);

        if (graph.get_id(cur_handle) != rightmost_id) { // don't iterate past rightmost node!
            // look for all nodes connected to cur_handle that need to be added
            // looking to the left,
            graph.follow_edges(cur_handle, true, [&](const handle_t handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
            // looking to the right,
            graph.follow_edges(cur_handle, false, [&](const handle_t handle) {
                // mark the nodes to come as to_visit
                if (visited.find(graph.get_id(handle)) == visited.end()) {
                    to_visit.insert(graph.get_id(handle));
                }
            });
        }
        else
        {
            // cerr << "found the righmost node id. Here's the cur_handle: " << graph.get_id(cur_handle) << endl;
            found_rightmost = true;
        }
    }
    if (!found_rightmost)
    {
        // we never found rightmost! We probably were fed a rightmost_id that was actually to the left of the leftmost_id. Throw error.
        cerr << "error:[vg normalize] in function extract_subgraph, was passed snarl with leftmost_id" << leftmost_id;
        cerr << " and rightmost_id " << rightmost_id;
        cerr << ". However, rightmost_id was not found to the right of leftmost_id.";
        cerr << " Were the ids swapped?" << endl;
        exit(1);
    }
    return subgraph;
}

// Integrates the snarl into the _graph, replacing the snarl occupying the space between
// source_id and sink_id.
//      In the process, transfers any embedded paths traversing the old snarl into the new
//      snarl.
// Arguments:
//      _graph: the _graph in which we want to insert the snarl.
//      to_insert_snarl: a *separate* handle_graph from _graph, often generated from
//      MSAconverter. embedded_paths: a vector of paths, where each is a pair.
//                        pair.first is the first step_handle of interest in the
//                        old_embedded_path, and pair.second is the step_handle *after*
//                        the last step_handle of interest in the old_embedded_path (can
//                        be the null step at the end of the path.)
//                        Note: these paths will be altered to represent the way they
//                        overlap in the new snarl. Otherwise, they would be invalidated.
//      source_id: the source of the old (to be replaced) snarl in _graph
//      sink_id: the sink of the old (to be replaced) snarl in _graph.
// Return: a pair of node ids, representing source and sink of the newly integrated snarl.
pair<handle_t, handle_t> SnarlNormalizer::integrate_snarl(SubHandleGraph &old_snarl, 
    const HandleGraph &to_insert_snarl,
    vector<pair<step_handle_t, step_handle_t>>& embedded_paths, 
    const id_t source_id, const id_t sink_id, const bool backwards) 
{
    
    // TODO: debug_statement: Check to make sure that newly made snarl has only one start
    // and end.
    // TODO:     (shouldn't be necessary once we've implemented alignment with
    // leading/trailing special chars.) Identify old and new snarl start and sink
    pair<vector<handle_t>, vector<handle_t>> to_insert_snarl_defining_handles =
        debug_get_sources_and_sinks(to_insert_snarl);

    if (to_insert_snarl_defining_handles.first.size() > 1 ||
        to_insert_snarl_defining_handles.second.size() > 1) {
        cerr << "ERROR: newly made snarl from a snarl with source " << source_id
             << " has more than one start or end. # of starts: "
             << to_insert_snarl_defining_handles.first.size()
             << " # of ends: " << to_insert_snarl_defining_handles.second.size() << endl;
        exit(1);
    }


    /// Replace start and end handles of old _graph snarl with to_insert_snarl start and
    /// end, and delete rest of old _graph snarl:

    // add to_insert_snarl into _graph without directly attaching the snarl to the _graph
    // (yet).
    vector<handle_t> to_insert_snarl_topo_order =
        handlealgs::lazier_topological_order(&to_insert_snarl);

    // Construct a parallel new_snarl_topo_order to identify
    // paralogous nodes between to_insert_snarl and the new snarl inserted in _graph.
    vector<handle_t> new_snarl_topo_order;

    // integrate the handles from to_insert_snarl into the _graph, and keep track of their
    // identities by adding them to new_snarl_topo_order.
    for (handle_t to_insert_snarl_handle : to_insert_snarl_topo_order) {
        // cerr << "About to insert snarl handle from normalized graph of id, seq: "
        //      << to_insert_snarl.get_id(to_insert_snarl_handle) << " "
        //      << to_insert_snarl.get_sequence(to_insert_snarl_handle) << endl;

        handle_t graph_handle =
            _graph.create_handle(to_insert_snarl.get_sequence(to_insert_snarl_handle));
        // cerr << "here is the new snarl handle: " 
        //      << _graph.get_id(graph_handle) << " "
        //      << _graph.get_sequence(graph_handle) << endl;
        new_snarl_topo_order.push_back(graph_handle);
    }
    // cerr << "finished inserting the snarls from to_insert_snarl into normalized graph." << endl;

    // Connect the newly made handles in the _graph together the way they were connected
    // in to_insert_snarl:
    for (int i = 0; i < to_insert_snarl_topo_order.size(); i++) {
        to_insert_snarl.follow_edges(
            to_insert_snarl_topo_order[i], false, [&](const handle_t snarl_handle) {
                // get topo_index of nodes to be connected to _graph start handle
                auto it = find(to_insert_snarl_topo_order.begin(),
                               to_insert_snarl_topo_order.end(), snarl_handle);
                int topo_index = it - to_insert_snarl_topo_order.begin();

                // connect _graph start handle
                _graph.create_edge(new_snarl_topo_order[i],
                                   new_snarl_topo_order[topo_index]);
            });
    }

    // save the source and sink values of new_snarl_topo_order, since topological order is
    // not necessarily preserved by move_path_to_snarl. Is temporary b/c we need to
    // replace the handles with ones with the right id_t label for source and sink later
    // on.
    id_t temp_snarl_leftmost_id = _graph.get_id(new_snarl_topo_order.front());
    id_t temp_snarl_rightmost_id = _graph.get_id(new_snarl_topo_order.back());
    if (new_snarl_topo_order.size() == 1)
    {
        // in case the normalized snarl is only one handle in size, split it into two.
        // This allows the front to be renamed after the source, and the end after the sink.
        std::pair<handle_t, handle_t> split_handle = _graph.divide_handle(new_snarl_topo_order.back(), 1);
        temp_snarl_leftmost_id = _graph.get_id(split_handle.first);
        temp_snarl_rightmost_id = _graph.get_id(split_handle.second);
    }
    // cerr << "the temp source id: " << temp_snarl_leftmost_id << endl;
    // cerr << "the temp sink id: " << temp_snarl_rightmost_id << endl;

    // Add the neighbors of the source and sink of the original snarl to the new_snarl's
    // source and sink.
    // source integration:
    if (!backwards)
    {
    _graph.follow_edges(
        _graph.get_handle(source_id), true, [&](const handle_t prev_handle) {
            _graph.create_edge(prev_handle, _graph.get_handle(temp_snarl_leftmost_id));
        });
    _graph.follow_edges(
        _graph.get_handle(sink_id), false, [&](const handle_t next_handle) {
            _graph.create_edge(_graph.get_handle(temp_snarl_rightmost_id), next_handle);
        });
    }
    else 
    {
        _graph.follow_edges(
        _graph.get_handle(source_id), false, [&](const handle_t next_handle) {
            _graph.create_edge(_graph.get_handle(temp_snarl_rightmost_id), next_handle);
        });
    _graph.follow_edges(
        _graph.get_handle(sink_id), true, [&](const handle_t prev_handle) {
            _graph.create_edge(prev_handle, _graph.get_handle(temp_snarl_leftmost_id));
        });
    }

    //todo: uncomment this for when I run normalize on region-segregated normalized snarls. It only has bugs when the regions are not separated.
    // For each path of interest, move it onto the new_snarl.
    for (int i = 0; i != embedded_paths.size(); i++)
    {
        pair<bool, bool> path_spans_left_right;
        path_spans_left_right.first = (_graph.get_id(_graph.get_handle_of_step(embedded_paths[i].first)) == source_id);
        path_spans_left_right.second = (_graph.get_id(_graph.get_handle_of_step(embedded_paths[i].second)) == sink_id); // not get_previous_step because stop_inclusive=true for extract haplotype paths.

        embedded_paths[i] = move_path_to_new_snarl(embedded_paths[i], temp_snarl_leftmost_id, temp_snarl_rightmost_id, path_spans_left_right, !backwards, make_pair(source_id, sink_id));
    }

    // Destroy the old snarl.
    old_snarl.for_each_handle([&](const handle_t handle) 
    {
        _graph.destroy_handle(handle);
    });

    // Replace the source and sink handles with ones that have the original source/sink id
    // (for compatibility with future iterations on neighboring top-level snarls using the
    // same snarl manager. Couldn't replace it before b/c we needed the old handles to
    // move the paths.
    handle_t new_leftmost_handle;
    handle_t new_rightmost_handle;
    if (!backwards) 
    {
        // cerr << "!backwards" << endl;
        // cerr << "overwriting node id " << temp_snarl_leftmost_id <<  " with " << source_id << " (which is source_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_leftmost_id)) << endl;
        new_leftmost_handle = overwrite_node_id(temp_snarl_leftmost_id, source_id);
        // cerr << "overwriting node id " << temp_snarl_rightmost_id <<  " with " << sink_id << " (which is sink_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_rightmost_id)) << endl;
        new_rightmost_handle = overwrite_node_id(temp_snarl_rightmost_id, sink_id);
    }
    else
    {
        // cerr << "backwards" << endl;
        // cerr << "overwriting node id " << temp_snarl_leftmost_id <<  " with " << sink_id << " (which is sink_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_leftmost_id)) << endl;
        new_leftmost_handle = overwrite_node_id(temp_snarl_leftmost_id, sink_id);
        // cerr << "overwriting node id " << temp_snarl_rightmost_id <<  " with " << source_id << " (which is source_id)." << " has sequence " << _graph.get_sequence(_graph.get_handle(temp_snarl_rightmost_id)) << endl;
        new_rightmost_handle = overwrite_node_id(temp_snarl_rightmost_id, source_id);
    }  

    pair<handle_t, handle_t> new_left_right = make_pair(new_leftmost_handle, new_rightmost_handle);

    return new_left_right;
}

/**
 * Deletes the given node, replacing it with a new node that has the desired
 * new_node_sequence.
 *
 * @param  {id_t} node_id     : The old node id, to be replaced with a new node id.
 * @param  {string} node_sequence    : The node id for the new node. Cannot be currently
 *                              in use in the graph.
 * @return {handle_t}        : The new handle, in the same position as the original handle
 *                              in the graph, but with the new node_sequence.
 */
handle_t SnarlNormalizer::replace_node_using_sequence(const id_t old_node_id, const string new_node_sequence, MutablePathDeletableHandleGraph& graph)
{
    handle_t old_handle = graph.get_handle(old_node_id);
    if (new_node_sequence.size() == 0) //special case of creating an empty node.
    {
        pair<handle_t, handle_t> split = graph.divide_handle(old_handle, graph.get_sequence(old_handle).size());
        graph.follow_edges(split.first, true, [&](const handle_t prev_handle) 
        {
            graph.create_edge(prev_handle, split.second);
        });
        graph.destroy_handle(split.first);
        return split.second;
    }
    handle_t new_handle = graph.create_handle(new_node_sequence);

    // move the edges:
    graph.follow_edges(old_handle, true, [&](const handle_t prev_handle) 
    {
        graph.create_edge(prev_handle, new_handle);
    });
    graph.follow_edges(old_handle, false, [&](const handle_t next_handle)
    {
        graph.create_edge(new_handle, next_handle);
    });

    // move the paths:
    graph.for_each_step_on_handle(old_handle, [&](step_handle_t step) 
    {
        handle_t properly_oriented_old_handle = graph.get_handle_of_step(step); 
        if (graph.get_is_reverse(properly_oriented_old_handle) != graph.get_is_reverse(new_handle))
        {
            new_handle = graph.flip(new_handle);
        }
        graph.rewrite_segment(step, graph.get_next_step(step), vector<handle_t>{new_handle});
    });

    // delete the old_handle:
    graph.destroy_handle(old_handle);
    return new_handle;
}

/**
 * Deletes the given handle's underlying node, and returns a new handle to a new node 
 * with the desired node_id
 * 
 * @param  {id_t} handle     : The old node id, to be replaced with a new node id.
 * @param  {id_t} node_id    : The node id for the new node. Cannot be currently in use in
 *                              the graph.
 * @return {handle_t}        : The new handle, in the same position as the original handle
 *                              in the graph, but with the new node_id.
 */
handle_t SnarlNormalizer::overwrite_node_id(const id_t old_node_id, const id_t new_node_id)
{
    handle_t old_handle = _graph.get_handle(old_node_id);
    handle_t new_handle = _graph.create_handle(_graph.get_sequence(old_handle), new_node_id);

    // move the edges:
    _graph.follow_edges(old_handle, true, [&](const handle_t prev_handle) 
    {
        _graph.create_edge(prev_handle, new_handle);
    });
    _graph.follow_edges(old_handle, false, [&](const handle_t next_handle)
    {
        _graph.create_edge(new_handle, next_handle);
    });

    // move the paths:
    _graph.for_each_step_on_handle(old_handle, [&](step_handle_t step) 
    {
        handle_t properly_oriented_old_handle = _graph.get_handle_of_step(step); 
        if (_graph.get_is_reverse(properly_oriented_old_handle) != _graph.get_is_reverse(new_handle))
        {
            new_handle = _graph.flip(new_handle);
        }
        _graph.rewrite_segment(step, _graph.get_next_step(step), vector<handle_t>{new_handle});
    });

    // delete the old_handle:
    _graph.destroy_handle(old_handle);
    return new_handle;
}

/**
 * Updates the changes that need making to the gbwt after the graph is finished being
 * normalized, so that an updated gbwt can be made.
 * @param  {list<string>} old_paths : the paths in the gbwt that need moving to the new
 * graph.
 * @param  {HandleGraph} new_snarl  : the normalized portion of the graph. Probably a 
 * subhandlegraph.
 */
void SnarlNormalizer::log_gbwt_changes(const vector<pair<gbwt::vector_type, string>>& source_to_sink_gbwt_paths, const HandleGraph &new_snarl){
    //todo: move Aligner to initialization of object, since I'm not supposed to make a new one each time I do alignments.
    Aligner aligner = Aligner();
    // cerr << "in log_gbwt_changes" << endl;
    // cerr << old_paths.size() << endl;
    for (auto path : source_to_sink_gbwt_paths)
    {
        Alignment alignment;
        alignment.set_sequence(path.second);
        aligner.align_global_banded(alignment, new_snarl,0, false);
        // cerr << "gbwt path being sent to new graph: " << endl;
        // cerr << "ALIGNMENT PATH FOR " << path << ":" << endl;
        gbwt::vector_type alignment_full_path;
        for (auto mapping : alignment.path().mapping())
        {
            // gbwt::Node::encode(id, is_reverse)
            // mapping.position().
            alignment_full_path.emplace_back(gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse()));
            // cerr << "mapping.position().node_id() " << mapping.position().node_id() << _graph.get_sequence(_graph.get_handle( mapping.position().node_id() )) << endl;
        }
        _gbwt_changelog.emplace_back(path.first, alignment_full_path);
        // cerr << pb2json(alignment.path()) << endl << alignment.query_position() << endl << alignment.path().mapping().begin() << endl << endl;
        // alignment.path().mapping()
    }
    
    // use banded global aligner. optimizations for finidng one perfect match from source to sink.

}




/** Used to help move_path_to_snarl map paths from an old snarl to its newly
 * normalized counterpart. In particular, ensures that any paths which touch the
 * source and/or sink of the old snarl still do so in the new snarl (which is
 * important to ensure that we don't break any paths partway through the snarl.)
 *
 * @param  {HandleGraph} _graph         : the _graph that contains the old and new snarl
 * nodes.
 * @param  {id_t} new_source_id        : the node id of the newly created source.
 * @param  {id_t} new_sink_id          : the node id of the newly created sink.
 * @param  {bool} touching_source      : true if the path is connected to the old
 * source.
 * @param  {bool} touching_sink        : true if the path is connected to the old
 * sink.
 * @param  {handle_t} path_start : proposed source for the path in the new snarl.
 * @param  {handle_t} path_end   : proposed sink for the path in the new snarl.
 * @return {bool}                      : true if the path satisfies the requirement
 * that, if the original path covered the old source or sink, the new path also covers
 * the same respective nodes in the new snarl.
 */
bool SnarlNormalizer::source_and_sink_handles_map_properly(
    const HandleGraph &graph, const id_t new_source_id, const id_t new_sink_id,
    const bool touching_source, const bool touching_sink, const handle_t path_start,
    const handle_t path_end) {

    bool path_map = false;
    // cerr << "touching source? " << touching_source << "touching_sink" << touching_sink
    //      << "source is source?" << (graph.get_id(path_start) == new_source_id)
    //      << " sink is sink: " << (graph.get_id(path_end) == new_sink_id) << endl;
    if (touching_source && touching_sink) {
        path_map = ((graph.get_id(path_start) == new_source_id) &&
                    (graph.get_id(path_end) == new_sink_id));
    } else if (touching_source) {
        path_map = (graph.get_id(path_start) == new_source_id);
    } else if (touching_sink) {
        path_map = (graph.get_id(path_end) == new_sink_id);
    } else {
        path_map = true;
    }
    // cerr << "path_map " << path_map << endl;
    return path_map;
}

// Determines whether some subsequence in a handle satisfies the condition of being
// the beginning of a path.
//      If the path_seq is longer than the handle_seq, only checks subsequences that
//      reach from the beginning/middle of the handle_seq to the end. If path_seq is
//      shorter than handle_seq, checks for any substring of length path_seq within
//      the handle_seq, as well as substrings smaller than length path_seq that extend
//      beyond the current handle.
// Arguments:
//      handle_seq: the sequence in the handle we're trying to identify as a
//      start_of_path_seq. path_seq: the sequence in the path we're trying to find
//      starting points for in handle_seq
// Return: a vector of all potential starting index of the subsequence in the
// handle_seq.
vector<int> SnarlNormalizer::check_handle_as_start_of_path_seq(const string &handle_seq,
                                                               const string &path_seq) {
    vector<int> possible_start_indices;
    // If the handle_seq.size <= path_seq.size, look for subsequences reaching from
    // beginning/middle of handle_seq to the end - where path_seq may run off the end
    // of this handle to the next in the snarl.
    if (handle_seq.size() <= path_seq.size()) {
        // iterate through all possible starting positions in the handle_seq.
        for (int handle_start_i = 0; handle_start_i < handle_seq.size();
             handle_start_i++) {
            int subseq_size = handle_seq.size() - handle_start_i;
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
    }
    // if handle_seq.size > path_seq.size, look for any subsequence within handle_seq
    // of path_seq.size, as well as any subsequence smaller than path_seq reaching
    // from middle of handle_seq to the end of handle_seq.
    else {
        // first, search through all handle_seq for any comparable subsequence of
        // path_seq.size. Note: only differences between this for loop and above for
        // loop is that handle_start_i stops at (<= path_seq.size() -
        // handle_seq.size()), and subseq.size() = path_seq.size()
        for (int handle_start_i = 0;
             handle_start_i <= (handle_seq.size() - path_seq.size()); handle_start_i++) {
            int subseq_size = path_seq.size();
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
        // second, search through the last few bases of handle_seq for the beginning
        // of path_seq. Note: nearly identical for loop to the one in "if
        // (handle_seq.size()
        // <= path_seq.size())"
        for (int handle_start_i = (handle_seq.size() - path_seq.size() + 1);
             handle_start_i < handle_seq.size(); handle_start_i++) {
            int subseq_size = handle_seq.size() - handle_start_i;
            // The path_seq subsequence of interest is from 0 to subseq_size;
            // The handle_seq subsequence of interest starts at handle_start_i
            // and ends at the end of the handle_seq (len subseq_size).
            // if compare returns 0, the substring matches.
            if (path_seq.compare(0, subseq_size, handle_seq, handle_start_i,
                                 subseq_size) == 0) {
                possible_start_indices.push_back(handle_start_i);
            }
        }
    }
    // Note: if we passed through the above check without returning anything, then
    // there isn't any satisfactory subsequence and we'll return an empty vector.
    return possible_start_indices;
}

// ------------------------------ DEBUG CODE BELOW:
// ------------------------------------------

// Returns pair where pair.first is a vector of all sources of the given _graph and
// path.second is all the sinks of the given _graph. If _graph is a subhandlegraph of a
// snarl, there should only be one source and sink each.
pair<vector<handle_t>, vector<handle_t>>
SnarlNormalizer::debug_get_sources_and_sinks(const HandleGraph &graph) {
    // cerr << "debug_get_source_and_sinks" << endl;
    vector<handle_t> sink;
    vector<handle_t> source;

    // identify sources and sinks
    graph.for_each_handle([&](const handle_t handle) {
        //todo: debug_statements in code below:
        // cerr << "identifying if " << graph.get_id(handle) << "is a source/sink." <<endl;
        bool is_source = true, is_sink = true;
        // cerr << "handles to the left: ";
        graph.follow_edges(handle, true, [&](const handle_t prev) {
            // cerr << graph.get_id(prev) << endl;
            is_source = false;
            return false;
        });
        // cerr << "handles to the right: ";
        graph.follow_edges(handle, false, [&](const handle_t next) {
            // cerr << graph.get_id(next) << endl;
            is_sink = false;
            return false;
        });

        if (is_source) {
            // cerr<< "determined is_source" << endl;
            source.push_back(handle);
        }
        if (is_sink) {
            // cerr<< "determined is_sink" << endl;
            sink.emplace_back(handle);
        }
    });
    return pair<vector<handle_t>, vector<handle_t>>(source, sink);
}

void SnarlNormalizer::get_all_gbwt_sequences(id_t source_id, id_t sink_id, bool backwards)
{
    cerr << "in get_all_gbwt_sequences" << endl;
    id_t leftmost_id;
    id_t rightmost_id;
    if (backwards) {
        leftmost_id = sink_id;
        rightmost_id = source_id;
    }
    else {
        leftmost_id = source_id;
        rightmost_id = sink_id;
    }

    
    SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);

    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, source_id, sink_id, backwards);

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<id_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();
    // unordered_set<string> hap_strs = format_handle_haplotypes_to_strings(get<0>(gbwt_haplotypes));
    for (auto hap : get<0>(gbwt_haplotypes))
    {
        cerr << "new hap:" << endl;
        for (auto handle : hap)
        {
            cerr << _gbwt_graph.get_id(handle) << endl;
        }
        // cerr << hap_str << endl;
    }
}

void SnarlNormalizer::make_one_edit(id_t leftmost_id, id_t rightmost_id) 
{
    cerr << "this is running. Source: " << leftmost_id << " sink: " << rightmost_id << endl;
    ///first, find a path through the snarl. In gbwt_graph.
    vector<handle_t> path;
    path.push_back(_gbwt_graph.get_handle(leftmost_id));

    gbwt::SearchState first_state = _gbwt_graph.get_state(path.back());

    cerr << "testing out a path example." << endl;
    gbwt::SearchState next_state = first_state;
    cerr << "next state: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) << endl;
    while (_gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) != rightmost_id)
    {
        _gbwt_graph.follow_paths(next_state, [&](const gbwt::SearchState one_next_state) -> bool {
                                     next_state = one_next_state;
                                     return false;
                                 });
        path.push_back(_gbwt_graph.node_to_handle(next_state.node));
        cerr << "next state: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_state.node)) << endl;
    }
    cerr << "finished a path example." << endl;

    ///store the old path as a gbwt::vector_type. and save the path sequence as a string, for the replacement_handle:
    string replacement_string;
    gbwt::vector_type old_path;
    for (auto old_handle : path)
    {
        if (_gbwt_graph.get_id(old_handle)!= leftmost_id && _gbwt_graph.get_id(old_handle)!= rightmost_id)
        {
           replacement_string.append(_gbwt_graph.get_sequence(old_handle)); 
        }
        old_path.push_back(gbwt::Node::encode(_gbwt_graph.get_id(old_handle), _gbwt_graph.get_is_reverse(old_handle)));
    }


    ///remove that path, except for the source/sink.
    for (handle_t gbwt_handle : path)
    {
        if (_gbwt_graph.get_id(gbwt_handle) != leftmost_id && _gbwt_graph.get_id(gbwt_handle) != rightmost_id )
        {
            handle_t normal_handle = _graph.get_handle(_gbwt_graph.get_id(gbwt_handle), _gbwt_graph.get_is_reverse(gbwt_handle));
            _graph.destroy_handle(normal_handle);
        }
    }

    ///replace it with another
    handle_t replacement_handle = _graph.create_handle(replacement_string);
    _graph.create_edge(_graph.get_handle(leftmost_id), replacement_handle);
    _graph.create_edge(replacement_handle, _graph.get_handle(rightmost_id));
    
    ///log that change in the gbwt_changelog.

    gbwt::vector_type new_path;
    new_path.push_back(gbwt::Node::encode(leftmost_id, _gbwt_graph.get_is_reverse(path.front())));
    new_path.push_back(gbwt::Node::encode(_graph.get_id(replacement_handle), false));
    new_path.push_back(gbwt::Node::encode(rightmost_id, _gbwt_graph.get_is_reverse(path.back())));

    _gbwt_changelog.push_back(make_pair(old_path, new_path));
}

void SnarlNormalizer::output_msa(const id_t leftmost_id, const id_t rightmost_id)
{
    SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);

    SnarlSequenceFinder sequence_finder = SnarlSequenceFinder(_graph, snarl, _gbwt_graph, leftmost_id, rightmost_id, false);

    tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<id_t>>
            gbwt_haplotypes = sequence_finder.find_gbwt_haps();

    unordered_set<string> haplotypes = format_handle_haplotypes_to_strings(_graph, _gbwt_graph, get<0>(gbwt_haplotypes));

    //true, because I want to print the msa to cout.
    poa_source_to_sink_haplotypes(haplotypes, true);
}

// // we want the sizes here to be in bases, not in handles/number of snarls.
// int _max_region_size = 1000;
// int _max_snarl_spacing = 1000;
/*
 *Goal of snarl stats: give a sense of how big the average snarl is, as well as whatever else I'd like to know about snarls. Used to figure out the recommended clustering gap, cap on alignment size, etc.
*/
// void SnarlNormalizer::snarl_stats(const vector<const Snarl *> &snarl_roots) {    
//     //first, I just want to know what the distribution of snarl sizes are.
//     vector<int> snarl_sizes;
//     int sum_size = 0;
//     int avg_size = 0;
//     int size_more_than_thousand = 0;
//     int size_more_than_three_thousand = 0;
//     for (auto snarl : snarl_roots)
//     {
//         SubHandleGraph snarl_graph = SnarlNormalizer::extract_subgraph(_graph, region.first, region.second);
//         int total_size = 0;
//         snarl_graph.for_each_handle([&](const handle_t handle){
//             total_size += _graph.get_sequence(handle).size();
//         });
//         sum_size += total_size;
//         if (total_size > 1000)
//         {
//             size_more_than_thousand += 1;
//         }
//         if (total_size > 3000)
//         {
//             size_more_than_three_thousand += 1;
//         }
//         snarl_sizes.push_back(total_size);

//     }
//     avg_size = sum_size/snarl_sizes.size(); 

//     cerr << "sum_size " << sum_size << endl;
//     cerr << "avg_size " << avg_size << endl;
//     cerr << "size_more_than_thousand " << size_more_than_thousand << endl;
//     cerr << "size_more_than_three_thousand " << size_more_than_three_thousand << endl;
    
    
// }


}
}
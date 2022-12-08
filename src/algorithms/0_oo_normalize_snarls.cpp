#pragma once // TODO: remove this, to avoid warnings + maybe bad coding practice?

#include "0_oo_normalize_snarls.hpp"
#include "0_snarl_sequence_finder.hpp"
#include <string>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

#include <abpoa/abpoa.h>


#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"

#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../msa_converter.hpp"
#include "../snarls.hpp"
#include "../vg.hpp"

#include "../types.hpp"
#include "extract_containing_graph.hpp"

#include "multipath_mapper.hpp"

// #define USE_CALLGRIND
// #ifdef USE_CALLGRIND
// #include <valgrind/callgrind.h>
// #endif

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
                                 const int max_snarl_spacing,
                                 const int threads,
                                 const int max_alignment_size, /*= MAX_INT*/
                                 const string path_finder, /*= "GBWT"*/
                                 const string alignment_algorithm, /*= "sPOA"*/
                                 const bool disable_gbwt_update, /*= false*/
                                 const bool debug_print /*= false*/)
: _graph(graph), _gbwt(gbwt), _max_alignment_size(max_alignment_size),
      _max_handle_size(max_handle_size), _max_region_size(max_region_size), _max_snarl_spacing(max_snarl_spacing), _threads(threads), _path_finder(path_finder), _gbwt_graph(gbwt_graph),
       _alignment_algorithm(alignment_algorithm), _disable_gbwt_update(disable_gbwt_update), _debug_print(debug_print){}

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


/**
 * Iterates over all top-level snarls in _graph, and normalizes them.
 * @param snarl_stream file stream from .snarl.pb output of vg snarls
*/
tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> SnarlNormalizer::normalize_snarls(const vector<const Snarl *>& snarl_roots) {
    // snarl_stats(snarl_roots);
    //Extend each of the normalize_regions to encompass multiple snarls from snarl_roots, if max_region_size > 1.
    //normalize_regions is a pair indicating: leftmost_id, rightmost_id. Note this is equivalent to a snarl of source_id=leftmost_id, sink_id=rightmost_id, backward=false.
    // select a subset of snarls for debugging purposes.
    // cerr << "snarl_roots.size()" << snarl_roots.size() << endl;
    // auto first = snarl_roots.begin() + 475064;
    // auto last = snarl_roots.begin() + 475078;
    // vector<const Snarl *> small_test(first, last);
    // cerr << "snarls in the sample: " << endl;
    // for (auto root : small_test) 
    // {
    //     cerr << root->start().node_id() << " " << root->end().node_id() << endl;
    // }
    // vector<pair<id_t, id_t>> normalize_regions = get_normalize_regions(small_test);

    // cerr << "before making normalize_regions" << endl;
    vector<pair<id_t, id_t>> normalize_regions = get_normalize_regions(snarl_roots);
    // cerr << "after making normalize_regions" << endl;
    // vector<pair<id_t, id_t>> debug_normalize_regions;
    
    // int target_node = 803849;
    // cerr << "normalize regions containing: " << target_node << endl;
    // for (int i = 0 ; i != normalize_regions.size() ; i++)
    // {
    //     if (normalize_regions[i].first <= target_node && target_node <= normalize_regions[i].second)
    //     {
    //         debug_normalize_regions.push_back(normalize_regions[i]);
    //         cerr << "iter " << i << ": " << normalize_regions[i].first << " " << normalize_regions[i].second << endl;
    //     }
    // }
    // normalize_regions = debug_normalize_regions;
    // cerr << "size of normalize_regions: " << normalize_regions.size() << endl;
    // cerr << "debug exit" << endl;
    // exit(0);

    int num_snarls_normalized = 0;
    int total_num_snarls_skipped = 0;
    
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
     * Further error records:
     *      6) snarl is trivial (either one or two nodes only), so we skipped normalizing them.
     *      7) snarl has handles not represented in the gbwt, and so would be dropped if normalized.
    */ 
    int error_record_size = 8;
    vector<int> one_snarl_error_record(error_record_size, 0);
    vector<int> full_error_record(error_record_size, 0);

    // pair<int, int> snarl_sequence_change;


    // todo: debug_code
    // int stop_size = 10000;
    // int num_snarls_touched = 0;

    // int skip_first_few = 1;
    // int skipped = 0;

    // get_all_gbwt_sequences(1, 15, false);
    // cerr << "has node before? " << _gbwt_graph.has_node(gbwt::Node::encode(7405162, false)) << endl;

    int snarl_num = 0;
    // for (auto roots : snarl_roots) 
    for (auto region : normalize_regions) 
    {
        // if (region.first != 996838)
        // {
        //     continue;
        // }
        snarl_num++;
        // if (snarl_num > start_snarl_num)
        // {
        //     cerr << "in for loop" << endl;
        // }

        // if (start_snarl_num != 0 || end_snarl_num != 0)
        // {
        //     if (start_snarl_num > end_snarl_num)
        //     {
        //         cerr << "ERROR: start_snarl_num >= end_snarl_num" << endl;
        //         exit(1);
        //     }
        //     if (snarl_num < start_snarl_num)
        //     {
        //         continue;
        //     }
        //     if (snarl_num > end_snarl_num)
        //     {
        //         cerr << "reached limit of snarl_num" << endl;
        //         break;
        //     }
        // }
        

        // if (snarl_num==1 || snarl_num==2 || snarl_num==3)
        // {
        //     continue;
        //     cerr << "not running" << endl;
        // }
        
        // if (skipped < skip_first_few){
        //     skipped++;
        //     continue;
        // }
        
        // if (num_snarls_touched == stop_size){
        //     cerr << "breakpoint here" << endl;
        //     break;
        // } else {
        //     num_snarls_touched++;
        // }

        // //todo: debug_print:
        // // if (roots->start().node_id()!=3775521)
        // if (snarl_num!=22530 && snarl_num!=22529)
        // {
        //     continue;
        // }
        // else
        // {
        //     cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        // }200: " << region.first << " and sink at: " << region.second << endl;
        if (_debug_print)
        {
            // cerr << "normalizing region number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
            cerr << "normalizing region number " << snarl_num << " with source at: " << region.first << " and sink at: " << region.second << endl;
        }
        else if (snarl_num==1 || snarl_num%10000 == 0)
        {
            // cerr << "normalizing region number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
            cerr << "normalizing region number " << snarl_num << " with source at: " << region.first << " and sink at: " << region.second << endl;
        }

        // if (!roots->start().backward())
        // {
        //     cerr << "not backward" << endl;
        //     // make_one_edit(1, 15);
        //     // get_all_gbwt_sequences(roots->start().node_id(), roots->end().node_id(), roots->start().backward());
        //     make_one_edit(roots->start().node_id(), roots->end().node_id());
        // }
        // else
        // {
        //     cerr << "backward" << endl;
        //     make_one_edit(roots->end().node_id(), roots->start().node_id());
        // }  

        // }
        // cerr << "normalizing snarl number " << snarl_num << " with source at: " << roots->start().node_id() << " and sink at: " << roots->end().node_id() << endl;
        // cerr << "seq from snarl number " << snarl_num << " with source at: " << _graph.get_sequence(_graph.get_handle(roots->start().node_id())) << " and sink at: " << _graph.get_sequence(_graph.get_handle(roots->end().node_id())) << endl;
        // cerr << "graph range of original graph " << _gbwt_graph.min_node_id() << " " <<  _gbwt_graph.max_node_id() <<endl ;
        // // cerr << "gbwt graph investigation: " << _gbwt_graph.has_node(test) << endl;
        // cerr << _gbwt_graph.get_length(_gbwt_graph.get_handle(roots->start().node_id())) << endl;
        // cerr << "seq from snarl number (using gbwt) " << snarl_num << " with source at: " << _gbwt_graph.get_sequence(_gbwt_graph.get_handle(roots->start().node_id())) << " and sink at: " << _gbwt_graph.get_sequence(_gbwt_graph.get_handle(roots->end().node_id())) << endl;

        // cerr << "backwards value? " << roots->start().backward() << endl;
        // if (roots->start().node_id() == 3881494) {
            // cerr << "root backwards?" << roots->start().backward() << endl;
            // cerr << "disambiguating snarl #"
            //         << (num_snarls_normalized + total_num_snarls_skipped)
            //         << " source: " << roots->start().node_id()
            //         << " sink: " << roots->end().node_id() << endl;

            // one_snarl_error_record = normalize_snarl(roots->start().node_id(), roots->end().node_id(), roots->start().backward(), snarl_num);
            one_snarl_error_record = normalize_snarl(region.first, region.second, false, snarl_num);
            // one_snarl_error_record = normalize_snarl(region.second, region.first, false, snarl_num);
            if (!(one_snarl_error_record[0] || one_snarl_error_record[1] ||
                    one_snarl_error_record[2] || one_snarl_error_record[3] ||
                    one_snarl_error_record[6] || one_snarl_error_record[7])) {
                // if there are no errors, then we've successfully normalized a snarl.
                num_snarls_normalized += 1;
                // track the change in size of the snarl.
                // snarl_sequence_change.first += one_snarl_error_record[4];
                // snarl_sequence_change.second += one_snarl_error_record[5];
                // cerr << "normalized snarl starting at: " << roots->start().node_id() << endl;
            } else {
                // else, there was an error. Track which errors caused the snarl to not
                // normalize.
                // note: the ints 4 and 5 are ignored here b/c they're for
                // recording the changing size of snarls that are successfully normalized.
                for (int i = 0; i < error_record_size; i++) {
                    if ( i != 4 && i != 5)
                    {
                        full_error_record[i] += one_snarl_error_record[i];
                    }
                }
                total_num_snarls_skipped += 1;
            }
            
            // //todo: debug_statement for extracting snarl of interest.
            // VG outGraph;
            // pos_t source_pos = make_pos_t(roots->start().node_id(), false, 0);
            // vector<pos_t> pos_vec;
            // pos_vec.push_back(source_pos);
            // algorithms::extract_containing_graph(&_graph, &outGraph, pos_vec, roots->end().node_id() - roots->start().node_id() + 2);
            // outGraph.serialize_to_ostream(cout);
            // break;
        // }

    }
    
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
        // if (!handlealgs::is_acyclic(&region_graph)) {
        //     // skip cyclic snarls, just like in actual normalization.
        //     cerr << "ERROR: I shouldn't have a cyclic region at this point?!" << endl;
        //     continue;
        // }
        // bool skip_this_subgraph = false;
        // region_graph.for_each_handle([&](handle_t handle){
        //     const gbwt::BidirectionalState debug_state = _gbwt_graph.get_bd_state(_gbwt_graph.get_handle(region_graph.get_id(handle)));
        //     // const gbwt::SearchState debug_state = _gbwt_graph.get_state(_gbwt_graph.get_handle(region_graph.get_id(handle))); //todo: consider trying this, in both orientations?
        //     if (debug_state.empty()) //then we skipped this graph during normalization, so we should skip it here as well.
        //     {
        //         // cerr << "skipping" << endl;
        //         skip_this_subgraph = true;
        //         return;
        //     }
        //     // else
        //     // {
        //     //     return false;
        //     // }
        // }, true);
        // if (skip_this_subgraph)
        // {
        //     continue;
        // }


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

    // //todo: remove this debug code:
    // int unskipped_snarls_not_found_after_norm = 0;
    // for (auto unskipped_snarl : _unskipped_snarls)
    // {
    //     if (snarls_found_after_norm.find(unskipped_snarl) == snarls_found_after_norm.end())
    //     {
    //         unskipped_snarls_not_found_after_norm += 1;
    //     }
    // }
    // cerr << "unskipped_snarls_not_found_after_norm: " << unskipped_snarls_not_found_after_norm << endl;
    // cerr << "_unskipped_snarls.size(): " << _unskipped_snarls.size() << endl;

    // int snarls_found_after_norm_that_were_skipped = 0;
    // for (auto snarl_found_after_norm : snarls_found_after_norm)
    // {
    //     if (_unskipped_snarls.find(snarl_found_after_norm) == _unskipped_snarls.end())
    //     {
    //         snarls_found_after_norm_that_were_skipped += 1;
    //     }
    // }
    // cerr << "snarls_found_after_norm_that_were_skipped: " << snarls_found_after_norm_that_were_skipped << endl;
    // cerr << "snarls_found_after_norm.size(): " << snarls_found_after_norm.size() << endl;

    int num_top_snarls_tracked = 10; // hardcoded for debugging convenience. //todo: change?
    // partial_sort(_snarl_size_changes.begin(), _snarl_size_changes.begin() + num_top_snarls_tracked,_snarl_size_changes.end(), []());

    // Detect the top best snarl changes:
    auto best_compare = [](const pair<pair<id_t, id_t>, pair<int, int>> left, const pair<pair<id_t, id_t>, pair<int, int>> right) 
    {
        return (right.second.second - right.second.first) >= (left.second.second - left.second.first);
    };
    vector<pair<pair<id_t, id_t>, pair<int, int>>> best_snarl_changes;


    for (pair<pair<id_t, id_t>, pair<int, int>> region : _snarl_size_changes)
    {
        // cerr << endl << "all best_snarl_changes: " << endl;
        // for (auto change : best_snarl_changes)
        // {
        //     cerr << change.second.first << " , " << change.second.second << endl;
        // }
        // cerr << "deciding whether or not to insert following value: " << region.second.first << " , " << region.second.second << endl;
        if (best_snarl_changes.size() < num_top_snarls_tracked)
        {
            // cerr << "inserted because the best_snarl_changes isn't big enough." << endl;
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
        }
        else if ((region.second.second - region.second.first) < (best_snarl_changes.back().second.second - best_snarl_changes.back().second.first)) //todo: check that this is proper comparison.
        {
            // cerr << "inserted and replaced an item. Item replaced: " << best_snarl_changes.back().second.first << " , " << best_snarl_changes.back().second.second << endl;
            best_snarl_changes.insert(upper_bound(best_snarl_changes.begin(), best_snarl_changes.end(), region, best_compare), region);
            best_snarl_changes.pop_back();
        }
        // else
        // {
        //     cerr << "didn't insert an item" << endl;
        // }
        // debug_count++;
        // if (debug_count > 30)
        // {
        //     break;
        // }
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

    // // find the snarls with the biggest changes in size before/after normalization: //note: doesn't work because I'm trying to sort a map.
    // sort(_snarl_size_changes.begin(), _snarl_size_changes.end(), [](const pair<pair<id_t, id_t>, pair<int, int>> &left, const pair<pair<id_t, id_t>, pair<int, int>> &right) {
    //     // return (5) >= (3);
    //     return (left.second.second - left.second.first) >= (right.second.second - right.second.first);
    // });




    // // NOTE: MY SILLY (but more efficient?) CODE THAT ONLY TRACKS TOP FEW SNARL CHANGES BELOW, SEE ABOVE FOR FULL-SORT:
    // // the snarls with the biggest changes in size before/after normalization:
    // vector<pair<pair<id_t, id_t>, pair<int, int>>> best_snarl_changes; // 0 is snarl with biggest change; 1 the snarl with 2nd biggest change, etc.
    // int num_top_snarls_tracked = 3; // hardcoded for debugging convenience. //todo: change?
    // for (auto snarl_size_change : _snarl_size_changes)
    // {
    //     int i = 0; // my position iterating through the best_snarl_changes list.
    //     // int best_replacement = -1; // replacements will only be valid if >= 0.
    //     while (i < num_top_snarls_tracked)
    //     {
    //         if (i < best_snarl_changes.size()) // if best_snarl_changes is not yet filled, and we haven't found a better spot to fill, fill this spot.
    //         {
    //             best_snarl_changes.push_back(snarl_size_change);
    //             break;
    //         }
    //         else if ((best_snarl_changes[i].second.second - best_snarl_changes[i].second.first) < (snarl_size_change.second.second - snarl_size_change.second.first))
    //         {
    //             best_snarl_changes.insert(best_snarl_changes.begin() + i, snarl_size_change);
    //             if (best_snarl_changes.size() > num_top_snarls_tracked)
    //             {
    //                 best_snarl_changes.pop_back();
    //             }
    //         }
    //         i++; 
    //     }
    //     if (best_snarl_changes.size() < num_top_snarls_tracked)
    //     {
    //         best_snarl_changes.push_back(snarl_size_change);
    //     }
    //     else
    //     {
            
    //     }
    // }

    //todo: replace error_record system with object-wide variables that are edited whenever needed in normalize_snarl(). If I feel like making the code easier to read and edit.
    // float percent_snarl_sequence_reduced = static_cast<float>((post_norm_net_snarl_size/_pre_norm_net_snarl_size)*100.0);
    float percent_snarl_sequence_change = ((static_cast<float>(post_norm_net_snarl_size)-static_cast<float>(_pre_norm_net_snarl_size))/static_cast<float>(_pre_norm_net_snarl_size))*100.0;
    // cerr.precision(2);
    cerr << endl
         << "normalization arguments:" << endl
         << "aligner (-A): " << _alignment_algorithm << endl
         << "max normalization region size (-k): " << _max_region_size << " snarls" << endl //todo: change to "bases" if/when I change _max_region_size's metric.
         << "max snarl spacing (-i): " << _max_snarl_spacing << endl //todo: add units. Handles? bases? I don't recall.
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

    tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> gbwt_update_items = make_tuple(_gbwt_graph, _gbwt_changelog, _gbwt);
    return gbwt_update_items;
    // cerr << "_gbwt_changelog size: " << _gbwt_changelog.size() << endl;

    // for (auto& entry : _gbwt_changelog)
    // {
    //     // if (entry.first.size() - entry.second.size() == 2) //export a good snarl for debugging.
    //     // {
    //     cerr << "old_nodes: " << endl;
    //     for (auto node : entry.first)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    
    //     cerr << "new_nodes: " << endl;
    //     for (auto node : entry.second)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    //     // }
        
    // //     if (false)
    // //     {

    // //         cerr << "hi" << endl;
    // //     }
    // //     for (auto new_node : entry.second)
    // //     {
    // //         cerr << "graph contains the new nodes? " << _graph.contains(_gbwt_graph.get_id(_gbwt_graph.node_to_handle(new_node))) << endl;

    // //     }
    //     // cerr << "size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
    //     // if (_gbwt.find(entry.first.begin(), entry.first.end()).empty())
    //     // {
    //     //     cerr << "found empty path!" << endl;
    //     //     cerr << "size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
            
    //     //     for (auto item : entry.first)
    //     //     {
    //     //         cerr << "seq in node: " << _gbwt_graph.get_sequence(_gbwt_graph.node_to_handle(item)) << endl;
    //     //         entry.second.push_back(item);
    //     //     }
    //     // }
    //     // cerr << _gbwt.find(entry.first.begin(), entry.first.end()) << endl;
    //     // for (auto item : entry.first)
    //     // {
    //     //     cerr << "item" << endl;
    //         // cerr << "seq in node: " << _gbwt_graph.get_sequence(_gbwt_graph.node_to_handle(item)) << endl;
    //         // entry.second.push_back(item);
    //     // }
    // }
    // for (auto& entry : _gbwt_changelog)
    // {
    //     cerr << "2nd size of entries: " << entry.first.size() << " " << entry.second.size() << endl;
    // }
    // cerr << _gbwt.empty() << _gbwt_changelog.empty() << endl;
    // _gbwt_changelog.push_back()
    
    // cerr << "contents of gbwt changelog: " << endl;
    // for (auto& entry : _gbwt_changelog)
    // {
    //     // if (entry.first.size() - entry.second.size() == 2) //export a good snarl for debugging.
    //     // {
    //     cerr << "old_nodes: " << endl;
    //     for (auto node : entry.first)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    
    //     cerr << "new_nodes: " << endl;
    //     for (auto node : entry.second)
    //     {
    //         cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(node)) << " " << _gbwt_graph.get_is_reverse(_gbwt_graph.node_to_handle(node)) << endl;
    //     }
    // }
    
    // _gbwt_graph.weakly_connected_components(_graph);
    //todo: use weakly_connected below for multithreading
    // handlegraph::algorithms::weakly_connected_components(_graph)
    // gbwtgraph::weakly_connected_components(_graph);


    
    // #ifdef USE_CALLGRIND
    //     // We want to profile the alignment, not the loading.
    //     CALLGRIND_START_INSTRUMENTATION;
    // #endif
    // if (!_disable_gbwt_update)
    // {
    //      gbwt::GBWT output_gbwt = apply_gbwt_changelog();
    //     // gbwt::GBWT output_gbwt = rebuild_gbwt(_gbwt, _gbwt_changelog);
    //     cerr << "finished generating gbwt." << endl;
    //     return output_gbwt;
    // }
    // else
    // {
    //     gbwt::GBWT empty_gbwt;
    //     return empty_gbwt;
    // }
    
    // //todo: debug-code for checking that I can build the gbwt_graph:
    // cerr << "making new gbwt graph." << endl;
    // gbwtgraph::GBWTGraph output_gbwt_graph = gbwtgraph::GBWTGraph(output_gbwt, _graph);
    // cerr << "new gbwt graph created." << endl;

    
    // cerr << "output have second path?" << endl;
    // for (auto& entry : _gbwt_changelog)
    // {
    //     cerr << _gbwt.find(entry.second.begin(), entry.second.end()).empty() << endl;
    // }


    //todo: debug_statement for extracting snarl of interest.
    // VG outGraph;
    // pos_t source_pos = make_pos_t(269695, false, 0);
    // vector<pos_t> pos_vec;
    // pos_vec.push_back(source_pos);
    // algorithms::extract_containing_graph(&_graph, &outGraph, pos_vec, 1000);
    // _graph = outGraph;
    // vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph *>(outGraph.get()), cout);
    // outGraph.serialize_to_ostream(cout);



}

////////////////////////////////////////////////////////////
/// Use the gbwt_changelog to perform paralellized rebuild_gbwt():
////////////////////////////////////////////////////////////


gbwt::GBWT SnarlNormalizer::apply_gbwt_changelog()
{
    // is the changelog in _graph ids? Or gbwt ids? If gbwt ids, use the gbwt weakly_connected_components. If it's the graph, I need to compute the connected components before changing it.
    // Uh, but actually the gbwt algorithm is gonna be expecting gbwt ids. So basically, use the gbwt info.
    
    // vector<unordered_set<nid_t>> components = gbwtgraph::weakly_connected_components(&_gbwt_graph);
    vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&_gbwt_graph);
    // vector<unordered_set<nid_t>> components = handlegraph::algorithms::weakly_connected_components(&_graph); // this was using the handlegraph algorithm, but I think we want the gbwt's view of the connected components.

    std::unordered_map<nid_t, size_t> node_to_job = get_node_to_job(components);

    std::vector<RebuildJob> jobs = divide_changelog_into_jobs(node_to_job, components); 

    RebuildParameters rebuild_parameters = set_parameters();
    
    gbwt::GBWT output_gbwt = rebuild_gbwt(_gbwt, jobs, node_to_job, rebuild_parameters);
    //todo: remove temporary non-parallelized rebuild_gbwt.
    // gbwt::GBWT output_gbwt = rebuild_gbwt(_gbwt, _gbwt_changelog);
    return output_gbwt;
}

std::unordered_map<nid_t, size_t> SnarlNormalizer::get_node_to_job(const vector<unordered_set<nid_t>>& weakly_connected_components)
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

std::vector<RebuildJob> SnarlNormalizer::divide_changelog_into_jobs(const std::unordered_map<nid_t, size_t>& node_to_job, const vector<unordered_set<nid_t>>& weakly_connected_components)
{
    std::vector<RebuildJob> jobs(weakly_connected_components.size());
    // cerr << "jobs.size()" << jobs.size() << endl;   
    // cerr << "jobs[0].size() " << jobs[0].mappings.size() << endl;

    for (pair<gbwt::vector_type, gbwt::vector_type> change : _gbwt_changelog)
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

RebuildParameters SnarlNormalizer::set_parameters()
{
    RebuildParameters parameters;
    parameters.num_jobs = _threads; //todo: add option for this to constructor.
    parameters.show_progress = _debug_print;
    return parameters;
    //todo: implement Jouni's advice for the max_region_size and sample_interval. Advice:
    /*
    * From Jouni:
    * The defaults should be fine in most cases.
    * If the threads are particularly long and memory usage is not a problem, you may want to increase the batch size to ~20x the length of the longest threads.
    */
    // parameters.max_region_size = ???;
    // parameters.sample_interval = ???;
    // /// Maximum number of parallel construction jobs.
    // size_t num_jobs = 1;

    // /// Print progress information to stderr.
    // bool show_progress = false;

    // /// Size of the GBWT construction buffer in nodes.
    // gbwt::size_type max_region_size = gbwt::DynamicGBWT::INSERT_max_region_size;

    // /// Sample interval for locate queries.
    // gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
}


////////////////////////////////////////////////////////////
/// Data preprocessing, as encapsulated in get_normalize_regions():
////////////////////////////////////////////////////////////

//snarls_adjacent used to identify if two snarls overlap at one of their boundary nodes.
bool SnarlNormalizer::snarls_adjacent(const Snarl& snarl_1, const Snarl& snarl_2) 
{
    //does snarl_2 overlap at the start handle of snarl_1?
    handle_t start_h = _graph.get_handle(snarl_1.start().node_id());
    bool overlap = _graph.follow_edges(start_h, snarl_1.start().backward(), [&](handle_t potential_end){
        if (snarl_2.start().node_id() == snarl_1.start().node_id() || snarl_2.end().node_id() == snarl_1.start().node_id())
        {
            return true;
        }
        return false;
    });
    if (overlap) // to pass the result of the follow_edges to the return value of the snarls_adjacent:
    {
        return true;
    }
    //does snarl_2 overlap at the end handle of snarl_1?
    handle_t end_h = _graph.get_handle(snarl_1.end().node_id());
    overlap = _graph.follow_edges(end_h, !snarl_1.end().backward(), [&](handle_t potential_end){
        if (snarl_2.start().node_id() == snarl_1.end().node_id() || snarl_2.end().node_id() == snarl_1.end().node_id())
        {
            return true;
        }
        return false;
    });
    if (overlap) // to pass the result of the follow_edges to the return value of the snarls_adjacent:
    {
        return true;
    }

    return false;
}

bool SnarlNormalizer::is_trivial(const Snarl& snarl) {
    vector<id_t> nodes_adjacent_to_source;
    handle_t first_h = _graph.get_handle(snarl.start().node_id());

    _graph.follow_edges(first_h, snarl.start().backward(), [&](handle_t next)
    {
        nodes_adjacent_to_source.push_back(_graph.get_id(next));
    });

    if (nodes_adjacent_to_source.size() == 1)
    {
        if (nodes_adjacent_to_source.back() == snarl.end().node_id())
        {
            return true;
        }
        else
        {
            cerr << "error:[vg normalize] snarl with start" << snarl.start().node_id();
            cerr << " and end " << snarl.end().node_id();
            cerr << " is performing unexpectedly. The start has only one internal-to-snarl";
            cerr << " neighbor, and it is not the end node. Malformed snarl input?" << endl;
            exit(1);
        }
    }
    else
    {
        return false;
    }
}

vector<pair<id_t, id_t>> SnarlNormalizer::convert_snarl_clusters_to_regions(const vector<vector<const Snarl *> >& clusters) {
    // cerr << "convert_snarl_clusters_to_regions" << endl;
    vector<pair<id_t, id_t>> normalize_regions;
    // int debug_count = 0;
    // bool debug = false;
    for (auto cluster : clusters )
    {
        // if (debug_count == clusters.size() - 1)
        // {
        //     debug = true;
        //     cerr << "I found the important cluster in convert_snarl_clusters_to_regions! Here 'tis: " << endl;
        //     for (auto snarl : cluster)
        //     {
        //         cerr << "leftmost: " << (*snarl).start().node_id() << " rightmost: " << (*snarl).end().node_id() << endl;
        //     }
        // }
        
        pair<id_t, id_t> cur_region; 
        if (cluster.front()->start().backward())
        {
            cur_region.first = cluster.front()->end().node_id();
            cur_region.second = cluster.front()->start().node_id();
        } 
        else
        {
            cur_region.first = cluster.front()->start().node_id();
            cur_region.second = cluster.front()->end().node_id();
        }
        for (int i = 1; i != cluster.size(); i++)
        {
            //find which part of the cur_region to replace with this snarl:
            if (cur_region.first == cluster[i]->end().node_id())
            {
                cur_region.first = cluster[i]->start().node_id();
                // if (debug) 
                // {
                //     cerr << 1 << endl;
                //     cerr << cur_region.first << " " << cur_region.second << endl;
                // }
            }
            else if (cur_region.first == cluster[i]->start().node_id())
            {
                cur_region.first = cluster[i]->end().node_id();
                // if (debug) 
                // {
                //     cerr << 2 << endl;
                //     cerr << cur_region.first << " " << cur_region.second << endl;
                // }
            }
            else if (cur_region.second == cluster[i]->end().node_id())
            {
                cur_region.second = cluster[i]->start().node_id();
                // if (debug) 
                // {
                //     cerr << 3 << endl;
                //     cerr << cur_region.first << " " << cur_region.second << endl;
                // }
            }
            else if (cur_region.second == cluster[i]->start().node_id())
            {
                cur_region.second = cluster[i]->end().node_id();
                // if (debug) 
                // {
                //     cerr << 4 << endl;
                //     cerr << cur_region.first << " " << cur_region.second << endl;
                // }
            }
            else
            {
                // the current snarl is not adjacent to the previous snarl; we have left 
                // the snarl chain. 
                // This is supposed to be guaranteed not to happen, because of the way 
                // get_normalize_regions generates clusters. Raise an error.
                cerr << "error:[vg normalize] snarl cluster is not within only a single"<< 
                " connected component. There is likely a bug in whatever passed clusters"<<
                " to convert_snarl_clusters_to_regions." << endl;
                exit(1);
            }
            
        }
        // auto test = cur_region.first;
        // cur_region.first = cur_region.second;
        // cur_region.second = test;
        normalize_regions.push_back(cur_region);
        // debug_count ++;

    }
    return normalize_regions;
}

vector<vector<const Snarl *> > SnarlNormalizer::cluster_snarls(const vector<const Snarl *> &snarl_roots) {    
    vector<vector<const Snarl *> > snarl_clusters;
    auto cur_snarl = snarl_roots.begin();
    int trivial_count = 0;
    vector<const Snarl*> first_cluster;
    snarl_clusters.push_back(first_cluster);
    // cerr << "snarl_clusters.size()" << snarl_clusters.size() << endl;

    // int debug_count = 0;
    // bool debug = false;
    while (cur_snarl != snarl_roots.end())
    {
        // if (debug_count>53000)
        // {
        //     cerr << "Inserting snarl number " << debug_count << " into a cluster" <<endl;
        // }
        //todo: Find out if this graph extraction is highly inefficient. If so, remove the graph_extraction I perform here, or else integrate it into the 
        //todo:     overall normalizer in a way that prevents me from having to extract it twice (if that's efficient).
        auto deptr_snarl = **cur_snarl;
        id_t leftmost_id, rightmost_id;
        if (deptr_snarl.start().backward())
        {
            leftmost_id = deptr_snarl.end().node_id();
            rightmost_id = deptr_snarl.start().node_id();
        }
        else
        {
            leftmost_id = deptr_snarl.start().node_id();
            rightmost_id = deptr_snarl.end().node_id();
        }            
        SubHandleGraph snarl_graph = extract_subgraph(_graph, leftmost_id, rightmost_id);
        bool cyclic = !handlealgs::is_acyclic(&snarl_graph);
        // if (!handlealgs::is_acyclic(&snarl_graph))
        // {
        //     //This snarl is cyclic. Start new, empty cluster.
        //     vector<const Snarl*> new_cluster;
        //     snarl_clusters.push_back(new_cluster);
        // }

        // if (debug_count >= 54815)
        // {
        //     debug = true;
        // }
        // if (debug_count >= 54815){
        //     debug = true;
        //     auto frog = **cur_snarl;
        //     id_t leftmost_id, rightmost_id;
        //     cerr << "is the snarl backward? Start node: " << frog.start().backward() << " end node: " << frog.end().backward() << endl;
        //     if (frog.start().backward())
        //     {
        //         leftmost_id = frog.end().node_id();
        //         rightmost_id = frog.start().node_id();
        //     }
        //     else
        //     {
        //         leftmost_id = frog.start().node_id();
        //         rightmost_id = frog.end().node_id();
        //     }            
        //     cerr << "extracting snarl: " << endl;
        //     SubHandleGraph snarl = extract_subgraph(_graph, leftmost_id, rightmost_id);
        //     cerr << "extracting snarl_reverse: " << endl;
        //     SubHandleGraph snarl_reverse = extract_subgraph(_graph, rightmost_id, leftmost_id);

        //     // only normalize non-trivial snarls (i.e. not composed of just a source and sink.):
        //     int num_handles_in_snarl = 0;
        //     int num_seq_in_snarl = 0;
        //     snarl.for_each_handle([&](const handle_t handle){
        //         num_handles_in_snarl++;
        //         num_seq_in_snarl+= snarl.get_sequence(handle).size();
        //     });

        //     // only normalize non-trivial snarls (i.e. not composed of just a source and sink.):
        //     int num_handles_in_snarl_rev = 0;
        //     int num_seq_in_snarl_rev = 0;
        //     snarl_reverse.for_each_handle([&](const handle_t handle){
        //         num_handles_in_snarl_rev++;
        //         num_seq_in_snarl_rev+= snarl_reverse.get_sequence(handle).size();
        //     });

        //     cerr << "debug count: " << debug_count << " info on snarl: " << " with start " << frog.start().node_id() << " and end " << frog.end().node_id() << " and num_handles: " << num_handles_in_snarl << " and num_seq_in_snarls: " << num_seq_in_snarl << endl;
        //     cerr << " and num_handles from reverse-extracted snarl: " << num_handles_in_snarl_rev << " and num_seq_in_snarls_rev: " << num_seq_in_snarl_rev << endl;
        // }
        // if (debug)
        // {
        //     cerr << "debug: BEFORE while iter: all nodes in current snarl cluster:" << endl;
        //     for (auto snarl : snarl_clusters.back()) 
        //     {
        //         auto frog = *snarl;
        //         cerr << frog.start().node_id() << " " << frog.end().node_id() << endl;;

        //     }
        //     cerr << "end of cluster" << endl;
        // }
        // cerr << "snarl_clusters.size()" << snarl_clusters.size() << endl;
        // if (debug_count>=10 && snarl_clusters.size()>2)
        // {
        //     break;
        // }
        // debug_count++;
        // cerr << "cur_snarl" << (*cur_snarl)->start() << endl;


        // If we are considering extending a snarl cluster, first make sure that we 
        // haven't reached the end conditions of the current cluster.
        // batch size exceeded? Or snarls aren't part of the same connected component?
        // start new cluster, and trim the previous one of any trailing trivial snarls.
        if (snarl_clusters.back().size() != 0)
        {
            // if (debug){
            //     cerr << "first if" << endl;
            // }
            const Snarl prev_snarl = *snarl_clusters.back().back();
            if (snarl_clusters.back().size() == _max_region_size || !snarls_adjacent(prev_snarl, **cur_snarl) || trivial_count > _max_snarl_spacing || cyclic)
            {
                // if (cyclic)
                // {
                //     cerr << "cyclic in the nested if" << endl;
                // }
                //If we're here, we've reached the end condition for this cluster.
                //Trim the tirivial snarls from the end of the cluster:
                for (int i = 0; i < trivial_count; i++)
                {
                    // if (debug){
                    //     cerr << "trim" << endl;
                    // }
                    // cerr << " before pop_back, trivial count: " << trivial_count << ", i" << i << " snarl_clusters.size()" << snarl_clusters.size() << endl;
                    // snarl_clusters.pop_back();
                    snarl_clusters.back().pop_back();
                    // cerr << " after pop_back " << endl;
                }
                trivial_count = 0;

                // if (debug)
                // {
                //     cerr << " make new cluster" << endl;
                // }
                //start new cluster
                vector<const Snarl*> new_cluster;
                snarl_clusters.push_back(new_cluster);
            }

        }
        
        // special case: if we find a cyclic snarl and its the first snarl of a new cluster:
        if (snarl_clusters.back().size() == 0 && cyclic)
        {
            // cerr << "skipping a cyclic snarl." << endl;
            // Then just skip the cyclic snarl.
            cur_snarl++;
            // debug_count++;
            continue;
        }
        
        bool trivial = is_trivial(**cur_snarl);
        // if we have a snarl that isn't trivial, add it to the latest cluster.
        if (!trivial)
        {
            // if (debug)
            // {
            //     cerr << "not trivial" << endl;
            // }
            // cerr << "second if" << endl;
            snarl_clusters.back().push_back(*cur_snarl);
            //reset the trivial_count
            trivial_count=0;
        }
            
        // if we have a snarl that is trivial, but we have a cluster in progress, add it to the cluster.
        else if (snarl_clusters.back().size() > 0)
        {
        //     if (debug)
        //     {
        //         cerr << "is trivial, but belongs in cluster" << endl;
        //     }
            // cerr << "third if" << endl;
            snarl_clusters.back().push_back(*cur_snarl);
            trivial_count++;
        }
        // if we have a snarl that is trivial, and we have no cluster in progress, skip it. (no code needed.)
        // if (debug)
        // {
        //     cerr << "end of loop" << endl;
        // }
        cur_snarl++;
        // debug_count++;

        // if (debug)
        // {
        //     cerr << "debug: AFTER while iter: all nodes in current snarl cluster:" << endl;
        //     for (auto snarl : snarl_clusters.back()) 
        //     {
        //         auto frog = *snarl;
        //         cerr << frog.start().node_id() << " " << frog.end().node_id() << endl;;

        //     }
        //     cerr << "end of cluster" << endl;
        // }


    }
    //check to see if the snarl_cluster at the end is empty. If so, discard it.
    if (snarl_clusters.back().size() == 0)
    {
        // if (debug)
        // {    
        //     cerr << " popping the snarl cluster at the end. Size is: " << snarl_clusters.back().size() << endl;
        // }
        snarl_clusters.pop_back();
    }
    return snarl_clusters;
}

vector<pair<id_t, id_t>> SnarlNormalizer::get_single_snarl_normalize_regions(const vector<const Snarl *> &snarl_roots)
{
    vector<pair<id_t, id_t>> simple_normalize_regions;
    ////todo: debug_code: 
    // for (int i = 1; i != snarl_roots.size(); i++)
    // {
    //     cerr << "snarls adjacent between " << i-1 << " " << i << "?" << endl;
    //     cerr << snarls_adjacent(*snarl_roots[i-1], *snarl_roots[i]) << endl;
    // }
    for (auto roots : snarl_roots)
    {
        if (is_trivial(*roots))
        {
            continue;
        }
        id_t leftmost_id;
        id_t rightmost_id;
        if (roots->start().backward())
        {
            leftmost_id = roots->end().node_id();
            rightmost_id = roots->start().node_id();
        }
        else
        {
            leftmost_id = roots->start().node_id();
            rightmost_id = roots->end().node_id();
        }
        pair<id_t, id_t> normalize_region = make_pair(leftmost_id, rightmost_id);
        simple_normalize_regions.push_back(normalize_region);
    }
    return simple_normalize_regions;
}

vector<pair<id_t, id_t>> SnarlNormalizer::get_normalize_regions(const vector<const Snarl *> &snarl_roots) {
    // cerr << "get_normalize_regions" << endl;

    
    //I don't need max_interstitial trivials, I just need to trim off trivials if I hit the _max_region_size limit.
    
// vector<pair<id_t, id_t>> SnarlNormalizer::get_normalize_regions(vector<const Snarl *> snarl_roots) {

    //if batch size is 1, just return the nontrivial roots from snarl_roots.
    //todo: uncomment below if I want? Should be equivalent, but slightly more efficient than the general case. I'm currently using it for testing.
    // if (_max_region_size==1)
    // {
    //     return get_single_snarl_normalize_regions(snarl_roots);
    // }
    //otherwise, cluster snarls.
    // vector<pair<id_t, id_t>> simple_normalize_regions = get_single_snarl_normalize_regions(snarl_roots);
    // cerr << "debug-test: what is simple_normalize regions size? " << simple_normalize_regions.size() << endl;
    
    vector<vector<const Snarl *> > snarl_clusters = cluster_snarls(snarl_roots);


    // cerr << "debug: after snarl_clusters formed: all nodes in last snarl cluster:" << endl;
    // cerr << "note: are all the snarls in the cluster abut back-to-back? They should. Otherwise, the clusterer is removing important interstitial trivial snarls." << endl;
    // for (auto snarl : snarl_clusters.back()) 
    // {
    //     auto frog = *snarl;
    //     cerr << frog.start().node_id() << " " << frog.end().node_id() << " backwards: " << frog.start().backward() << " " << frog.end().backward() << " size: " << endl;;

    // }
    // cerr << "end of cluster" << endl;


    // cerr << "debug-test: is every snarl root inside a snarl cluster?" << endl;
    
    // set<const Snarl *> snarls_in_snarl_clusters;
    // for (auto cluster : snarl_clusters)
    // {
    //     for (auto snarl : cluster)
    //     {
    //         snarls_in_snarl_clusters.emplace(snarl);
    //     }
    // }

    // cerr << "number of snarls in snarl clusters: " << snarls_in_snarl_clusters.size() << endl;
    // cerr << "number of snarls in snarl_roots: " << snarl_roots.size() << endl;
    // cerr << "(if equal, then yes)." << endl;

    vector<pair<id_t, id_t>> regions = convert_snarl_clusters_to_regions(snarl_clusters);

    //For comparison across different clustering variables, an export of all handles in the graph:
    //DEBUG step one: save all nodes in (single-snarl-mode) regions to a file.
    // I will later read these regions back in and see if they are still contained in all the larger clusters.
    // ofstream all_region_nodes_file;
    // string filename = "nodes_in_normalizable_regions_in_k_" + to_string(_max_region_size) + "_i_" + to_string(_max_snarl_spacing) + "_m_" + to_string(_max_alignment_size) + ".txt"; 
    // string directory = "~/paten_lab/vg-team/vg/robin-graphs/yeast-subset/"; //todo: make dynamic if desired.
    // cerr << "writing to file named: " << filename << endl;
    // all_region_nodes_file.open(filename);
    // for (auto region : regions)
    // {
    //     SubHandleGraph region_graph = extract_subgraph(_graph, region.first, region.second);
        //     region_graph.for_each_handle([&](const handle_t handle){
        //         all_region_nodes_file << _graph.get_id(handle) << endl;
        //     });
    // }
    // all_region_nodes_file.close();

    // cerr << "last region: " << regions.back().first << " " << regions.back().second << endl;
    // cerr << "does the last region share first/last node id in common with cluster? " << endl;

    // cerr << "debug-test: is every snarl root contained within a region?" << endl; 
    // //todo: note that this debug-test only functions if node ids can be used to indicate 
    // //todo:     relative positions in the graph. If I'd like to implement this as a full 
    // //todo:     test, I'd need to use position of the nodes in a reference path in the graph.
    // pair<id_t, id_t> debug_last_region = regions.back();
    // SubHandleGraph debug_snarl = extract_subgraph(_graph, debug_last_region.first, debug_last_region.second);
    // for (auto snarl : snarl_clusters.back()) 
    // {
    //     cerr << debug_snarl.has_node(snarl->start().node_id()) << " " << debug_snarl.has_node(snarl->start().node_id()) << " ";
    // }
    // cerr << endl;
    
    // set<id_t> all_border_nodes;
    // for (auto region : regions)
    // {

    // }



    // int uncontained_count = 0;
    // int snarl_count = 0;
    // int trivial_count = 0;
    // for (auto snarl : snarl_roots)
    // {
    //     if (is_trivial(*snarl))
    //     {
    //         trivial_count++;
    //         continue;
    //     }
    //     bool contained = false;
    //     for (auto region : regions)
    //     {
    //         if (region.first <= region.second)
    //         {
    //             if (snarl->start().node_id() >= region.first && snarl->end().node_id() <= region.second)
    //             {
    //                 contained = true;
    //             }
    //         }
    //         if (region.first > region.second)
    //         {
    //             if (snarl->start().node_id() >= region.second && snarl->end().node_id() <= region.first)
    //             {
    //                 contained = true;
    //             }
    //         }

    //     }
    //     if (!contained)
    //     {
    //         uncontained_count++;
    //         if (!is_trivial(*snarl))
    //         {
    //             cerr << "error: snarl number " << snarl_count << " with start " << snarl->start().node_id() << " and end " << snarl->end().node_id() << " not contained in the normalize regions." << endl;
    //             // cerr << "snarl trivial? " << is_trivial(*snarl) << endl;
    //         }
    //     }
    //     snarl_count += 1;
    // }
    // cerr << "trivial count: " << trivial_count << endl;
    // int contained_count = 0;
    // cerr << "number of snarls not contained in a region: " << uncontained_count << endl;

    // cerr << "size of clusters:" << endl;
    // for (auto clust : snarl_clusters)
    // {
    //     cerr << clust.size() << endl;
    //     cerr << "contents of cluster: " << endl;
    //     for (auto snarl : clust)
    //     {
    //         cerr << snarl->start().node_id() << " " << snarl->end().node_id() << endl;
    //     }
    // }

    // cerr << "contents of regions" << endl;  
    // vector<pair<id_t, id_t>> regions = convert_snarl_clusters_to_regions(snarl_clusters);
    // for (auto reg : regions)
    // {
    //     cerr << reg.first << " " << reg.second << endl; 
    // }
    return regions;
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
    // if (backwards){
    //     // swap the source and sink ids. Essentially, this guarantees I treat the leftmost node in snarl as "source".
    //     // (although some adjustments for paths need be made)
    //     id_t swap_source = sink_id; //temp storage of sink_id value. 
    //     sink_id = source_id;
    //     source_id = swap_source; 
    // }

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
    */ 
    vector<int> error_record(8, 0);
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
        
        // cerr << "haps in haplotypes: " << endl;
        // for (string hap : get<0>(haplotypes))
        // {
        //     cerr << hap << endl;
        // }
        // Align the new snarl:
        VG new_snarl;
        if (_alignment_algorithm == "TCoffee")
        {
            new_snarl = align_source_to_sink_haplotypes(get<0>(haplotypes));
        }
        else if (_alignment_algorithm == "sPOA")
        {
            new_snarl = poa_source_to_sink_haplotypes(get<0>(haplotypes), snarl_num);
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
        bool single_stranded = handlealgs::is_single_stranded(&new_snarl);
        VG* single_stranded_snarl;
        if (!single_stranded) 
        {
            handlealgs::split_strands(&new_snarl, single_stranded_snarl);
            // handlealgs::SplitStrandOverlay(new_snarl)
        }
        else
        {
            single_stranded_snarl=&new_snarl;
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
        new_snarl.for_each_handle([&](const handle_t handle) {
            error_record[5] += new_snarl.get_sequence(handle).size();
            _snarl_size_changes[make_pair(leftmost_id, rightmost_id)].second += new_snarl.get_sequence(handle).size(); //todo: robin-review-and-fix this so that I can have correct change by initializing this here as zero, otherwise set it as the original size of the snarl when setting pair.first value.. 
        });
        force_maximum_handle_size(new_snarl);
        
        // integrate the new_snarl into the _graph, removing the old snarl as you go.
        // //todo: debug_statement
        // integrate_snarl(new_snarl, embedded_paths, sink_id, source_id);
        pair<handle_t, handle_t> new_left_right = integrate_snarl(snarl, new_snarl, embedded_paths, source_id, sink_id, backwards);
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
VG SnarlNormalizer::align_source_to_sink_haplotypes(
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
    VG snarl = myMSAConverter.make_graph();

    snarl.clear_paths();

    pair<vector<handle_t>, vector<handle_t>> source_and_sink =
        debug_get_sources_and_sinks(snarl);


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

// TODO: change the arguments to handles, which contain orientation within themselves.
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
    // cerr << "extract_subgraph has source and sink: " << source_id << " " << sink_id << endl; 
    // because algorithm moves left to right, determine leftmost and rightmost nodes.
    // id_t leftmost_id;
    // id_t rightmost_id;
    //// if snarl's "backwards," source is rightmost node, sink is leftmost.
    // if (backwards) 
    // {
    //     leftmost_id = sink_id;
    //     rightmost_id = source_id;
    // }
    // else 
    // {
    //     leftmost_id = source_id;
    //     rightmost_id = sink_id;
    // }
    // cerr << "extract_subgraph" << endl;
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
    const id_t source_id, const id_t sink_id, const bool backwards) {
    // cerr << "integrate_snarl" << endl;

    //todo: debug_statement
    // cerr << "\nhandles in to_insert_snarl:" << endl;
    // to_insert_snarl.for_each_handle([&](const handle_t handle) {
    //     cerr << to_insert_snarl.get_id(handle) << " "
    //          << to_insert_snarl.get_sequence(handle) << " ";
    //     cerr << "neighbors: ";
    //     to_insert_snarl.follow_edges(handle, false, [&](const handle_t next) {
    //         cerr << "     " << to_insert_snarl.get_id(next) << endl;
    //     });
    //     cerr << " \n";
    // });
    // cerr << endl;
    // Get old _graph snarl
    // SubHandleGraph old_snarl = extract_subgraph(_graph, source_id, sink_id, backwards);

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
        // //todo: debug_statement:
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
    // For each path of interest, move it onto the new_snarl.
    for (int i = 0; i != embedded_paths.size(); i++)
    {
        // //todo: debug_statement
        // cerr << "the new sink id: " << temp_snarl_rightmost_id << endl;
        // //todo: debug_statement
        // move_path_to_snarl(path, new_snarl_topo_order, temp_snarl_rightmost_id,
        //                    temp_snarl_leftmost_id, sink_id, source_id);
        // move_path_to_snarl(path, new_snarl_topo_order, temp_snarl_leftmost_id,
        //                    temp_snarl_rightmost_id, source_id, sink_id, backwards);
        // cerr << "is path backwards? " << backwards << endl;
        // cerr << "path first: " << _graph.get_id(_graph.get_handle_of_step(path.first)) << " step after path first: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_next_step(path.first))) << " path second: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << endl;
        // cerr << "source: " << source_id << " sink: " << sink_id << endl;
        // pair<bool, bool> path_spans_left_right;
        // path_spans_left_right.first = (!backwards && _graph.get_id(_graph.get_handle_of_step(path.first)) == source_id) || (backwards && _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == source_id);
        // path_spans_left_right.second = (!backwards && _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == sink_id) || (backwards && _graph.get_id(_graph.get_handle_of_step(path.first)) == sink_id);
        // cerr << "first: " << path_spans_left_right.first << "second: " << path_spans_left_right.second << endl;
        pair<bool, bool> path_spans_left_right;
        path_spans_left_right.first = (_graph.get_id(_graph.get_handle_of_step(embedded_paths[i].first)) == source_id);
        path_spans_left_right.second = (_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(embedded_paths[i].second))) == sink_id);

        embedded_paths[i] = move_path_to_new_snarl(embedded_paths[i], temp_snarl_leftmost_id, temp_snarl_rightmost_id, path_spans_left_right, !backwards, make_pair(source_id, sink_id));
    }

    // Destroy the old snarl.
    old_snarl.for_each_handle([&](const handle_t handle) 
    {
        // //todo: debug_statement these are the handles in old_snarl:
        // cerr << "destroying old_snarl handle: " << old_snarl.get_id(handle) << " with sequence: " << old_snarl.get_sequence(handle) << endl;
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
    poa_source_to_sink_haplotypes(haplotypes, 0, true);
}

}
}
#include "0_snarl_sequence_finder.hpp"

// #include "0_oo_normalize_snarls.hpp"
// #include "0_snarl_sequence_finder.hpp"
// #include <algorithm>
// #include <string>

// #include <seqan/align.h>
// #include <seqan/graph_align.h>
// #include <seqan/graph_msa.h>

// #include <gbwtgraph/gbwtgraph.h>

// #include "../gbwt_helper.hpp"
// #include "../handle.hpp"
// #include "../msa_converter.hpp"
// #include "../snarls.hpp"
// #include "../vg.hpp"
// #include "is_acyclic.hpp"

// #include "../types.hpp"
// #include "extract_containing_graph.hpp"

// #include <algorithm>

// #include "../msa_converter.hpp"
// #include "vg.hpp"

// #include "topological_sort.hpp"


#include <gbwtgraph/gbwtgraph.h>
#include "../handle.hpp"
#include "../subgraph.hpp"


namespace vg {
namespace algorithms{
SnarlSequenceFinder::SnarlSequenceFinder(const PathHandleGraph & graph, 
                               const SubHandleGraph &snarl,
                               const gbwtgraph::GBWTGraph &gbwt_graph, 
                               const id_t source_id, const id_t sink_id, const bool backwards)
    : _graph(graph), _gbwt_graph(gbwt_graph), _snarl(snarl), _source_id(source_id), _sink_id(sink_id), _backwards(backwards) {} //todo: since find_gbwt_haps uses neither _snarl nor _graph, maybe move various arguments from the object declaration to individual fxns?

// TODO: test that it successfully extracts any haplotypes that start/end in the middle of
// TODO:    the snarl.
/**
 * Finds all haplotypes in gbwt associated with snarl.
 * @param snarl The subhandlegraph of the snarl to be normalized.
 * @param source_id The source of the snarl.
 * @param sink_id The sink of the snarl.
 * @return A 3-tuple containing 1) a vector of haps stretching from source to sink, in 
 * vector<handle_t> format; 2) a second vector containing all other haps in snarl; 
 * 3) a vector of all handles oberved by the method.
*/
tuple<vector<vector<handle_t>>, vector<vector<handle_t>>, unordered_set<handle_t>>
SnarlSequenceFinder::find_gbwt_haps() {
    // cerr << "checking that the nodes of interest exist: " << endl;
    // for (int node : {803859, 803860, 803852, 803851, 803850})
    // {
    //     cerr << "node_id, has_node, node min, node max: " << node << " " << _gbwt_graph.has_node(gbwt::Node::encode(node, false)) << " " << _gbwt_graph.min_node_id() << " " << _gbwt_graph.max_node_id() << endl;
    // }
    // cerr << "803856 -> 803852 edge exist? " << _gbwt_graph.has_edge(_gbwt_graph.get_handle(803856), _gbwt_graph.get_handle(803852)) << endl;
    // If snarl has been fed to us backwards, run the algorithm with righmost_id as source
    // and vice-versa. Otherwise, keep source as leftmost_id.
    id_t leftmost_id = _source_id;
    id_t rightmost_id = _sink_id;
    if (_backwards) {
        leftmost_id = _sink_id;
        rightmost_id = _source_id;
    }
    
    /** 
     * haplotype_queue contains all started exon_haplotypes not completed yet.
     * Every time we encounter a branch in the paths, the next node down the path
     * Is stored here, along with the vector of handles that represents the path up
     * to the SearchState.
    */
    vector<pair<vector<handle_t>, gbwt::SearchState>> haplotype_queue;

    // source and sink handle for _gbwt_graph:
    handle_t source_handle = _gbwt_graph.get_handle(leftmost_id);
    handle_t sink_handle = _gbwt_graph.get_handle(rightmost_id);


    // place source in haplotype_queue.
    vector<handle_t> source_handle_vec(1, source_handle);
    gbwt::SearchState source_state = _gbwt_graph.get_state(source_handle);
    haplotype_queue.push_back(make_pair(source_handle_vec, source_state));

    //todo: debug statement:
    // cerr << "all handles to the right of source, accroding to _graph:" << endl;

    // touched_handles contains all handles that have been touched by the
    // depth first search below, for later use in other_haplotypes_to_strings, which
    // identifies paths that didn't stretch from source to sink in the snarl.
    unordered_set<handle_t> touched_handles{source_handle, sink_handle};

    // these haplotype vecs contains all "finished" haplotypes - those that were either 
    // walked to their conclusion, or until they reached the sink.
    vector<vector<handle_t>> haplotypes_from_source_to_sink;
    vector<vector<handle_t>> other_haplotypes;

    // sometimes a gbwt thread will indicate a connection between two handles that doesn't
    // actually exist in the _graph. These connections need to be ignored.
    unordered_set<edge_t> incorrect_connections;

    // for every partly-extracted thread, extend the thread until it either reaches
    // the sink of the snarl or the end of the thread.
    while (!haplotype_queue.empty()) {
        // todo: debug_statement
        // cerr << "haplotype queue: ";
        // cerr << "size of queue:" << haplotype_queue.size() << " " << endl;
        // for (auto hap : haplotype_queue) {
        //     cerr << "size: " << hap.first.size() << endl << "handle_ids: ";
        //     for (handle_t handle : hap.first) {
        //         cerr << _gbwt_graph.get_id(handle) << " ";
        //     }
        //     cerr << endl;
        // }

        // get a haplotype out of haplotype_queue to extend -
        // a tuple of (handles_traversed_so_far, last_touched_SearchState)
        pair<vector<handle_t>, gbwt::SearchState> cur_haplotype = haplotype_queue.back();
        haplotype_queue.pop_back();

        // get all the subsequent search_states that immediately follow the searchstate
        // from cur_haplotype.
        vector<gbwt::SearchState> next_searches;
        _gbwt_graph.follow_paths(cur_haplotype.second,
                                 [&](const gbwt::SearchState next_search) -> bool {
                                     next_searches.push_back(next_search);
                                    //  cerr << "adjacent handles to state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(cur_haplotype.second.node)) << " is " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
                                     return true;
                                 });
        // const gbwt::SearchState debug_state = _gbwt_graph.get_state(_gbwt_graph.node_to_handle(cur_haplotype.second.node));
        // _gbwt_graph.follow_paths(debug_state,
        //                         [&](const gbwt::SearchState debug_next) -> bool {
        //                             cerr << "and if I make a new search state, I get adjacent: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(debug_next.node)) << endl;
        //                             return true;
        //                         });
        // const gbwt::BidirectionalState debug_state_2 = _gbwt_graph.get_bd_state(_gbwt_graph.node_to_handle(cur_haplotype.second.node));
        // cerr << "bidirectionalstate follow paths: " << endl;
        // cerr << "follow paths false: " << endl;
        // _gbwt_graph.follow_paths(debug_state_2, false, [&](const gbwt::BidirectionalState& previous) -> bool
        // {
        //     cerr << "previous states forwards: " <<  _gbwt_graph.get_id(_gbwt_graph.node_to_handle(previous.forward.node))<< endl;
        //     cerr << "previous states backwards: " <<  _gbwt_graph.get_id(_gbwt_graph.node_to_handle(previous.backward.node))<< endl;
        //     return true;
        // });
        // cerr << "follow paths true: " << endl;
        // _gbwt_graph.follow_paths(debug_state_2, true, [&](const gbwt::BidirectionalState& previous) -> bool
        // {
        //     cerr << "previous states forwards: " <<  _gbwt_graph.get_id(_gbwt_graph.node_to_handle(previous.forward.node))<< endl;
        //     cerr << "previous states backwards: " <<  _gbwt_graph.get_id(_gbwt_graph.node_to_handle(previous.backward.node))<< endl;
        //     return true;
        // });
        // cerr << "Here is search state made using prefix at node: " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(cur_haplotype.second.node)) << endl;
        // gbwt::SearchState new_search =
        //     _gbwt_graph.index->prefix(cur_haplotype.second.node);
        
        // cerr << "what does the search state look like?" << endl;
        // cerr << "adjacent handles to state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(new_search.node)) << " is: " << endl;
        // _gbwt_graph.follow_paths(new_search,
        //                          [&](const gbwt::SearchState next_search) -> bool {
        //                              cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
        //                              return true;
        //                          });
        // cerr << endl;

        // if next_searches > 1, then we need to make multiple new haplotypes to be
        // recorded in haplotype_queue or one of the finished haplotype_handle_vectors.
        if (next_searches.size() > 1) {
            // for every next_search in next_searches, either create a new, extended
            // cur_haplotype to push into haplotype queue, or place in the
            // haplotypes_from_source_to_sink if haplotype extends to sink, or place in
            // the other_haplotypes if haplotype ends before reaching sink.
            for (gbwt::SearchState next_search : next_searches) {
                handle_t next_handle = _gbwt_graph.node_to_handle(next_search.node);
                // if (!_snarl.has_edge(cur_haplotype.first.back(), next_handle)) {
                //     if (incorrect_connections.find(
                //             _snarl.edge_handle(cur_haplotype.first.back(), next_handle)) ==
                //         incorrect_connections.end()) {
                //         cerr << "_snarl with source " << _source_id
                //             << " and sink " << _sink_id
                //             << " has a thread that incorrectly connects two nodes that "
                //                "don't have any edge connecting them. These two nodes are "
                //             << _gbwt_graph.get_id(cur_haplotype.first.back()) << " and "
                //             << _gbwt_graph.get_id(next_handle)
                //             << ". This thread connection will be ignored." << endl;
                //         incorrect_connections.emplace(
                //             _snarl.edge_handle(cur_haplotype.first.back(), next_handle));

                //         // todo: debug_statement
                //         // cerr << "next handle(s) of handle "
                //         //      << _snarl.get_id(cur_haplotype.first.back())
                //         //      << " according to _snarl:" << endl;
                //         // _snarl.follow_edges(cur_haplotype.first.back(), false,
                //         //                    [&](const handle_t handle) {
                //         //                        cerr << "\t" << _snarl.get_id(handle);
                //         //                    });
                //         // cerr << endl;
                //     }
                //     continue;
                // }
                // copy over the vector<handle_t> of cur_haplotype:
                vector<handle_t> next_handle_vec(cur_haplotype.first);

                // add the new handle to the vec:
                next_handle_vec.push_back(next_handle);

                // if new_handle is the sink, put in haplotypes_from_source_to_sink
                if (_gbwt_graph.get_id(next_handle) == rightmost_id) {
                    haplotypes_from_source_to_sink.push_back(next_handle_vec);
                } else // keep extending the haplotype!
                {
                    pair<vector<handle_t>, gbwt::SearchState> next_haplotype =
                        make_pair(next_handle_vec, next_search);
                    haplotype_queue.push_back(next_haplotype);
                }
                // next_handle will be touched.
                touched_handles.emplace(next_handle);
            }
        }
        // if next_searches is empty, the path has ended but not reached sink.
        else if (next_searches.empty()) {
            // We have reached the end of the path, but it doesn't reach the sink.
            // we need to add cur_haplotype to other_haplotypes.
            other_haplotypes.push_back(cur_haplotype.first);

        }
        // if next_handle is the "sink"/rightmost_id, put in haplotypes_from_source_to_sink
        else if (_gbwt_graph.get_id(
                     _gbwt_graph.node_to_handle(next_searches.back().node)) == rightmost_id) {
            // Then we need to add cur_haplotype + next_search to
            // haplotypes_from_source_to_sink.
            handle_t next_handle = _gbwt_graph.node_to_handle(next_searches.back().node);
            cur_haplotype.first.push_back(next_handle);
            haplotypes_from_source_to_sink.push_back(cur_haplotype.first);

            // touched next_search's handle
            touched_handles.emplace(next_handle);
        }
        // else, there is just one next_search, and it's not the end of the path.
        // just extend the search by adding (cur_haplotype + next_search to
        // haplotype_queue.
        else {
            // get the next_handle from the one next_search.
            handle_t next_handle = _gbwt_graph.node_to_handle(next_searches.back().node);

            // modify cur_haplotype with next_handle and next_search.
            cur_haplotype.first.push_back(next_handle);
            cur_haplotype.second =
                next_searches.back(); // there's only one next_search in next_searches.

            // put cur_haplotype back in haplotype_queue.
            haplotype_queue.push_back(cur_haplotype);
            touched_handles.emplace(next_handle);
        }
    }

    // Find any haplotypes starting from handles not starting at the source, but which
    // still start somewhere inside the snarl.
    vector<vector<handle_t>> haplotypes_not_starting_at_source =
        find_haplotypes_not_at_source();

    //todo: debug_statement:
    // cerr << "number of haplotypes starting at source, but ending before source: " << other_haplotypes.size() << endl;
    // cerr << "number of haplotypes starting after source: " << haplotypes_not_starting_at_source.size() << endl;
        

    // move haplotypes_not_starting_at_source into other_haplotypes:
    other_haplotypes.reserve(other_haplotypes.size() +
                             haplotypes_not_starting_at_source.size());
    move(haplotypes_not_starting_at_source.begin(),
         haplotypes_not_starting_at_source.end(), back_inserter(other_haplotypes));

    // //todo: debug_statement
    // cerr << "lets look through all the haplotypes after extraction:" << endl;
    // for (vector<handle_t> hap_vec : haplotypes_from_source_to_sink) {
    //     cerr << "new hap:" << endl;
    //     for (handle_t handle : hap_vec){
    //         // cerr << _gbwt_graph.get_id(handle) << " " << _graph.get_sequence(_graph.get_handle(_gbwt_graph.get_id(handle))) << endl;
    //         cerr << _gbwt_graph.get_id(handle) << " " << _gbwt_graph.get_sequence(handle) << endl; //note: this print *must* work if you want a working gbwt_graph.
    //     }
    // }

    return tuple<vector<vector<handle_t>>, vector<vector<handle_t>>,
                 unordered_set<handle_t>>{haplotypes_from_source_to_sink,
                                          other_haplotypes, touched_handles};
}


// Function that completes the traversal of a snarl along its haplotype threads, when there are
// handles connected to the snarl by threads that start after the source handle. (Threads
// that merely end before the sink handle are addressed in extract_gbwt_haplotypes).
// Arguments:
//      touched_handles: any handles found in the snarl so far.
// Returns:
//      a vector of haplotypes in vector<handle_t> format that start in the middle of the
//      snarl.
vector<vector<handle_t>>
SnarlSequenceFinder::find_haplotypes_not_at_source() {
    // cerr << "find_haplotypes_not_at_source begin" << endl;
    // If snarl has been fed to us backwards, run the algorithm with righmost_id as source
    // and vice-versa. Otherwise, keep source as leftmost_id.
    id_t leftmost_id = _source_id;
    id_t rightmost_id = _sink_id;
    if (_backwards) {
        leftmost_id = _sink_id;
        rightmost_id = _source_id;
    }
    // cerr << "leftmost_id: " << leftmost_id << endl;
    // cerr << "rightmost_id: " << rightmost_id << endl;
    // //todo: debug_statement
    // for (handle_t handle : touched_handles){
    //     cerr << "touched handles find_gbwt_haps: " << _graph.get_id(handle) << endl;
    // }
    // cerr << "find_haplotypes_not_at_source" << endl;

    /// Search every handle in touched handles for haplotypes starting at that point.
    // Any new haplotypes will be added to haplotype_queue.
    vector<pair<vector<handle_t>, gbwt::SearchState>> haplotype_queue;

    // Fully extended haplotypes (or haplotypes extended to the snarl's sink)
    // will be added to finished_haplotypes.
    vector<vector<handle_t>> finished_haplotypes;

    // We don't need to ever check the sink handle, since paths from the sink handle
    // extend beyond snarl.
    handle_t sink_handle = _gbwt_graph.get_handle(rightmost_id);
    // touched_handles.erase(sink_handle);

    // Nested function for making a new_search. Identifies threads starting at a given
    // handle and
    //      either adds them as a full haplotype (if the haplotype is one handle long) or
    //      makes a new entry to haplotype_queue.
    auto make_new_search = [&](handle_t handle) {
        // Are there any new threads starting at this handle?
        // if (_gbwt_graph.get_id(handle) == )
        // cerr << "make_new_search at handle: " << _gbwt_graph.get_id(handle) << endl;
        gbwt::SearchState new_search =
            _gbwt_graph.index->prefix(_gbwt_graph.handle_to_node(handle));
        
        // cerr << "what does the search state look like?" << endl;
        // cerr << "adjacent handles to state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(new_search.node)) << " is: " << endl;
        // _gbwt_graph.follow_paths(new_search,
        //                          [&](const gbwt::SearchState next_search) -> bool {
        //                              cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
        //                              return true;
        //                          });
        // cerr << endl;
        // gbwt::SearchState new_search_2 = _gbwt_graph.get_state(handle);
        // cerr << "adjacent handles to *standard* state of " << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(new_search_2.node)) << " is: " << endl;
        // _gbwt_graph.follow_paths(new_search_2,
        //                          [&](const gbwt::SearchState next_search) -> bool {
        //                              cerr << _gbwt_graph.get_id(_gbwt_graph.node_to_handle(next_search.node)) << endl;
        //                              return true;
        //                          });
        // cerr << endl;

        if (!new_search.empty()) {
            // Then add them to haplotype_queue.
            _gbwt_graph.follow_paths(
                new_search, [&](const gbwt::SearchState &next_search) -> bool {
                    handle_t next_handle = _gbwt_graph.node_to_handle(next_search.node);

                    /// check to make sure that the thread isn't already finished:
                    // if next_handle is the sink, or if this thread is only one handle
                    // long, then there isn't any need for this function to make a new 
                    // search from this. 
                    // This is because I already have a separate method 
                    // for storing haplotypes that have reached the sink in the while loop
                    // below.
                    if (next_handle != sink_handle ||
                        next_search == gbwt::SearchState()) {
                        // establish a new thread to walk along.
                        vector<handle_t> new_path;
                        new_path.push_back(handle);
                        new_path.push_back(next_handle);

                        pair<vector<handle_t>, gbwt::SearchState> mypair =
                            make_pair(new_path, next_search);

                        // add the new path to haplotype_queue to be extended.
                        haplotype_queue.push_back(make_pair(new_path, next_search));
                    }
                    return true;
                });
        }
    };

    _snarl.for_each_handle([&](const handle_t handle)
    {
        //todo: add quality to make_new_search that allows for handles that lead directly to sink next. 
        //todo:     This is important because I haven't dealt with the edge case where I add 
        //todo:     that path outside the while loop below; and the path might cover an edge that is uninvestigated in source_to_sink_haps.
        if (_snarl.get_id(handle) != _source_id && _snarl.get_id(handle) != _sink_id )
        {
            make_new_search(_gbwt_graph.get_handle(_snarl.get_id(handle), _snarl.get_is_reverse(handle)));
        }
    });

    /// Extend any paths in haplotype_queue, and add any newly found handles to to_search.
    /// Then, check to see if there are any new threads on handles in to_search.
    /// Extend those threads, and add any newly found handles to to_search,
    /// then search for threads again in to_search again... repeat until to_search remains
    /// emptied of new handles.

    // cerr << "haplotype_queue.size() " << haplotype_queue.size() << endl;
    // for tracking whether the haplotype thread is still extending:
    bool still_extending;
    while (!haplotype_queue.empty()) {
        // get a haplotype to extend out of haplotype_queue - a tuple of
        // (handles_traversed_so_far, last_touched_SearchState)
        pair<vector<handle_t>, gbwt::SearchState> cur_haplotype =
            haplotype_queue.back();
        haplotype_queue.pop_back();

        // get all the subsequent search_states that immediately follow the
        // searchstate from cur_haplotype.
        vector<gbwt::SearchState> next_searches;
        _gbwt_graph.follow_paths(cur_haplotype.second,
                                    [&](const gbwt::SearchState &next_search) -> bool {
                                        next_searches.push_back(next_search);
                                        return true;
                                    });

        for (gbwt::SearchState next_search : next_searches) {
            handle_t next_handle = _gbwt_graph.node_to_handle(next_search.node);

            // if next_search is empty, then we've fallen off the thread,
            // and cur_haplotype can be placed in finished_haplotypes as is for this
            // thread.
            if (next_search == gbwt::SearchState()) {
                finished_haplotypes.push_back(cur_haplotype.first);
            }

            // if next_search is on the sink_handle,
            // then cur_haplotype.first + next_search goes to finished_haplotypes.
            else if (_gbwt_graph.get_id(next_handle) == rightmost_id) {

                // copy over the vector<handle_t> of cur_haplotype:
                vector<handle_t> next_handle_vec(cur_haplotype.first);
                // add next_handle
                next_handle_vec.push_back(next_handle);
                // place in finished_haplotypes
                finished_haplotypes.push_back(next_handle_vec);


            }
            // otherwise, just place an extended cur_haplotype in haplotype_queue.
            else {
                // copy over cur_haplotype:
                pair<vector<handle_t>, gbwt::SearchState> cur_haplotype_copy =
                    cur_haplotype;
                // modify with next_handle/search
                cur_haplotype_copy.first.push_back(next_handle);
                cur_haplotype_copy.second = next_search;
                // place back in haplotype_queue for further extension.
                haplotype_queue.push_back(cur_haplotype_copy);
            }
        }
    }
    // cerr << "find_haplotypes_not_at_source end" << endl;
    return finished_haplotypes;
}

//////////////////////////////////////////////////////////////////////////////////////////
// embedded path finding:
//////////////////////////////////////////////////////////////////////////////////////////

// Finds all embedded paths that either start or end in a snarl (or both) defined by
// source_id, sink_id.
//      returns a vector of the embedded paths, where each entry in the vector is defined
//      by the pair of step_handles closest to the beginning and end of the path. If the
//      path is fully contained within the snarl, these step_handles will the be the
//      leftmost and rightmost handles in the path.
// Arguments:
//      _graph: a pathhandlegraph containing the snarl with embedded paths.
//      source_id: the source of the snarl of interest.
//      sink_id: the sink of the snarl of interest.
// Returns:
//      a vector containing all the embedded paths in the snarl, in pair< step_handle_t,
//      step_handle_t > > format. Pair.first is the first step in the path's range of
//      interest, and pair.second is the step *after* the last step in the path's range of
//      interest (can be the null step at end of path).
vector<pair<step_handle_t, step_handle_t>>
SnarlSequenceFinder::find_embedded_paths() {

    // key is path_handle, value is a step in that path from which to extend.
    unordered_map<path_handle_t, step_handle_t> paths_found;
    // cerr << "1: getting id of sink: " << (_graph).get_id((_graph).get_handle(_sink_id)) << endl;

    // look for handles with paths we haven't touched yet.
    _snarl.for_each_handle([&](const handle_t &handle) {
        // cerr << "looking for paths at handle " << _graph.get_id(handle) << endl; 
        vector<step_handle_t> steps = _graph.steps_of_handle(handle);
        // do any of these steps belong to a path not in paths_found?
        for (step_handle_t &step : steps) {
            path_handle_t path = _graph.get_path_handle_of_step(step);
            // cerr << "found a path. Is it new?" << endl;
            // If it's a step along a new path, save the first step to that path we find.
            // (The avoidance
            // of source and sink here is to ensure that we can properly check to see if
            // we've reached the end of an embedded path walking in any arbitrary
            // direction (i.e. source towards sink or sink towards source).
            //todo: should the following if statement only contain the first conditional? The other two conditionals don't do what the comment says, and also don't seem to make sense.
            if (paths_found.find(path) == paths_found.end() ||
                _graph.get_id(_graph.get_handle_of_step(paths_found.at(path))) == _source_id ||
                _graph.get_id(_graph.get_handle_of_step(paths_found.at(path))) == _sink_id) {
                // cerr << "found a new path." << endl;
                // then we need to mark it as found and save the step.
                paths_found.emplace(path, step);
            }
        }
    });
    // cerr << "2: getting id of sink: " << (_graph).get_id((_graph).get_handle(_sink_id)) << endl;

    /// for each step_handle_t corresponding to a unique path, we want to get the steps
    /// closest to both the end and beginning step that still remains in the snarl.
    // TODO: Note copy paste of code here. In python I'd do "for fxn in [fxn1, fxn2]:",
    // TODO      so that I could iterate over the fxn. That sounds template-messy in C++
    // tho'. Should I?
    int debug_counter = 0;
    // cerr << "2.3: getting sink directly: " << _sink_id << endl;
    vector<pair<step_handle_t, step_handle_t>> paths_in_snarl;
    // cerr << "2.3: getting sink directly: " << _sink_id << endl;
    for (auto &it : paths_found) {
        
        if ((_graph).get_id((_graph).get_handle(_sink_id)) != _sink_id)
        {
            cerr << "2.6: getting id of sink from handle isn't the same as sink: " << (_graph).get_id((_graph).get_handle(_sink_id)) << endl;
        }
        if (_sink_id != 997029)
        {
            cerr << "2) ";
            cerr << "sink id has changed after " << debug_counter << " iterations." << endl;
        }
        step_handle_t step = it.second;
        // path_in_snarl describes the start and end steps in the path,
        // as constrained by the snarl.
        pair<step_handle_t, step_handle_t> path_in_snarl;

        // Look for the step closest to the beginning of the path, as constrained by the
        // snarl.
        step_handle_t begin_in_snarl_step = step;
        id_t begin_in_snarl_id =
            _graph.get_id(_graph.get_handle_of_step(begin_in_snarl_step));

        while (((begin_in_snarl_id != _source_id)) &&
               _graph.has_previous_step(begin_in_snarl_step)) {
            begin_in_snarl_step = _graph.get_previous_step(begin_in_snarl_step);
            begin_in_snarl_id =
                _graph.get_id(_graph.get_handle_of_step(begin_in_snarl_step));
        }
        path_in_snarl.first = begin_in_snarl_step;

        // Look for the step closest to the end of the path, as constrained by the snarl.
        step_handle_t end_in_snarl_step = step;
        id_t end_in_snarl_id = _graph.get_id(_graph.get_handle_of_step(end_in_snarl_step));

        // while (end_in_snarl_id != source_id and end_in_snarl_id != sink_id and
        //        _graph.has_next_step(end_in_snarl_step)) {
        while ((end_in_snarl_id != _sink_id) and _graph.has_next_step(end_in_snarl_step)) {
            end_in_snarl_step = _graph.get_next_step(end_in_snarl_step);
            end_in_snarl_id = _graph.get_id(_graph.get_handle_of_step(end_in_snarl_step));
        }
        // Note: when adding the end step, path notation convention requires that we add
        // the null step at the end of the path (or the next arbitrary step, in the case
        // of a path that extends beyond our snarl.)
        path_in_snarl.second = _graph.get_next_step(end_in_snarl_step);

        paths_in_snarl.push_back(path_in_snarl);

        debug_counter++;
    }
    // cerr << "3: getting id of sink: " << (_graph).get_id((_graph).get_handle(_sink_id)) << endl;

    //todo: move the following to unit tests:
    unordered_set<string> path_names;
    for (auto path : paths_in_snarl) {
        // if (!(_graph.get_id(_graph.get_handle_of_step(path.first)) == _source_id)) {
        //     cerr << "************in UNIT_TEST for find_embedded_paths************" << endl;
        //     cerr << "in snarl with source: " << _source_id << " and sink " << _sink_id << ":" << endl;
        //     cerr << "path " << _graph.get_path_name(_graph.get_path_handle_of_step(path.first)) << " doesn't start at source of snarl. " << " source: " << _source_id << "; start of path: " << _graph.get_id(_graph.get_handle_of_step(path.first)) << endl;
        // }
        // if (!(_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == _sink_id)) {
        //     cerr << "************in UNIT_TEST for find_embedded_paths************" << endl;
        //     cerr << "in snarl with source: " << _source_id << " and sink " << _sink_id << ":" << endl;
        //     cerr << "path " << _graph.get_path_name(_graph.get_path_handle_of_step(path.second)) << " doesn't end at sink of snarl. " << " source: " << _sink_id << "; end of path: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) << endl;
        //     cerr << "note that the 'true' end of the path is one step further than the sink. Print statement above corrects for that convention." << endl;
        // }
        if (!(path_names.find(_graph.get_path_name(_graph.get_path_handle_of_step(path.first))) == path_names.end())) {
            cerr << "************in UNIT_TEST for find_embedded_paths************" << endl;
            cerr << "in snarl with source: " << _source_id << " and sink " << _sink_id << ":" << endl;
            cerr << "path " << _graph.get_path_name(_graph.get_path_handle_of_step(path.second)) << " has been found more than once in find_embedded_paths, when it should only have been extracted once. " << endl;
        }
        path_names.emplace(_graph.get_path_name(_graph.get_path_handle_of_step(path.first)));
    } 
    if ((path_names.size() == 0)) {
        if (_full_log_print)
        {
            cerr << "************in UNIT_TEST for find_embedded_paths************" << endl;
            cerr << "in snarl with source: " << _source_id << " and sink " << _sink_id << ":" << endl;
            cerr << "no embedded paths found in find_embedded_paths." << endl;
        }
    }
    for (auto path : paths_in_snarl) {
        step_handle_t cur_step = path.first;
        step_handle_t last_step = path.second;
        int path_size = 0;
        while (true) //todo: implement this style of iteration whenever looping until I find the last step of the path? Essentially, I'm concerned that the path_handle_ts aren't guaranteed to match. I can look that up in the definition of step_handle_t equals, I guess. It may not confirm, however.
        {
            path_size += _graph.get_sequence(_graph.get_handle_of_step(cur_step)).size();
            if (_graph.get_id(_graph.get_handle_of_step(cur_step)) == _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(last_step))))
            {
                break;
            }
            cur_step = _graph.get_next_step(cur_step);
        }
        cerr << path_size << endl;
    //     cerr << "path starts at source? " << (_graph.get_id(_graph.get_handle_of_step(path.first)) == _source_id) << endl;
    //     cerr << "path ends at sink? " << (_graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(path.second))) == _sink_id) << endl;

    //     cerr << "is a path a duplicate of one we've already extracted? " << (path_names.find(_graph.get_path_name(_graph.get_path_handle_of_step(path.first))) == path_names.end()) << endl;
    //     path_names.emplace(_graph.get_path_name(_graph.get_path_handle_of_step(path.first)));
    } 
    // cerr << "tested " << path_names.size() << " paths in UNIT_TEST." << endl;
    // cerr << "************END-UNIT_TEST for find_embedded_paths. Tested " << path_names.size() << " paths in UNIT_TEST.************"<< endl;


    return paths_in_snarl;
}

//////////////////////////////////////////////////////////////////////////////////////////
// exhaustive sequence extraction:
//////////////////////////////////////////////////////////////////////////////////////////

// Runs through the whole snarl and generates all possible strings representing walks
// from source to sink. Generates a combinatorial number of possible paths with splits
// in the snarl.
//todo: for consistency, have source_to_sink_exhaustive_path_finder return paths in format
//todo:     vector<vector<handle_t>> instead of vector<string>
pair<unordered_set<string>, unordered_set<handle_t>>
SnarlSequenceFinder::find_exhaustive_paths() {
    // cerr << "debug_graph_to_strings" << endl;
    unordered_set<handle_t> touched_handles;

    unordered_map<handle_t, vector<string>> sequences;
    vector<handle_t> sinks;
    unordered_map<handle_t, size_t> count;
    count.reserve(_snarl.get_node_count());     // resize count to contain enough buckets
                                               // for size of _snarl
    sequences.reserve(_snarl.get_node_count()); // resize sequences to contain enough
                                               // buckets for size of _snarl

    // identify sources and sinks //TODO: once we've established that this fxn works,
    // we can just use start_id and sink_id.
    _snarl.for_each_handle([&](const handle_t &handle) {
        bool is_source = true, is_sink = true;
        _snarl.follow_edges(handle, true, [&](const handle_t &prev) {
            is_source = false;
            return false;
        });
        _snarl.follow_edges(handle, false, [&](const handle_t &next) {
            is_sink = false;
            return false;
        });

        // base case for dynamic programming
        if (is_source) {
            count[handle] = 1;
            sequences[handle].push_back(
                _snarl.get_sequence(handle)); // TODO: presented in the handle's local
                                             // forward orientation. An issue?
        }
        if (is_sink) {
            sinks.emplace_back(handle);
        }
    });

    // count walks by dynamic programming
    bool overflowed = false;
    for (const handle_t &handle : handlealgs::lazier_topological_order(&_snarl)) {
        touched_handles.emplace(handle);
        size_t count_here = count[handle];
        vector<string> seqs_here = sequences[handle];

        _snarl.follow_edges(handle, false, [&](const handle_t &next) {
            size_t &count_next = count[next];
            string seq_next = _snarl.get_sequence(next);

            if (numeric_limits<size_t>::max() - count_here < count_next) {
                overflowed = true;
            }

            else {
                count_next += count_here;
                for (string seq : seqs_here) {
                    sequences[next].push_back(seq + seq_next);
                }
            }
        });
        /// TODO: figure out how to deal with overflow.
        // if (overflowed) {
        //     return numeric_limits<size_t>::max();
        // }
    }

    // total up the walks at the sinks
    size_t total_count = 0;
    for (handle_t &sink : sinks) {
        total_count += count[sink];
    }

    // all the sequences at the sinks will be all the sequences in the _snarl.
    unordered_set<string> walks;
    for (handle_t &sink : sinks) {
        for (string seq : sequences[sink]) {
            walks.emplace(seq);
        }
    }

    return make_pair(walks, touched_handles);
}


}
}
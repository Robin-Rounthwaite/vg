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

#include "topological_sort.hpp"


#include <gbwtgraph/gbwtgraph.h>
#include "../handle.hpp"
#include "../subgraph.hpp"


namespace vg {
namespace algorithms{
SnarlSequenceFinder::SnarlSequenceFinder(const SubHandleGraph &snarl,
                               const gbwtgraph::GBWTGraph &haploGraph, 
                               const id_t &source_id, const id_t &sink_id)
    : _haploGraph(haploGraph), _snarl(snarl), _source_id(source_id), _sink_id(sink_id) {}

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
    /** 
     * haplotype_queue contains all started exon_haplotypes not completed yet.
     * Every time we encounter a branch in the paths, the next node down the path
     * Is stored here, along with the vector of handles that represents the path up
     * to the SearchState.
    */
    vector<pair<vector<handle_t>, gbwt::SearchState>> haplotype_queue;

    // source and sink handle for _haploGraph:
    handle_t source_handle = _haploGraph.get_handle(_source_id);
    handle_t sink_handle = _haploGraph.get_handle(_sink_id);

    // place source in haplotype_queue.
    vector<handle_t> source_handle_vec(1, source_handle);
    gbwt::SearchState source_state = _haploGraph.get_state(source_handle);
    haplotype_queue.push_back(make_pair(source_handle_vec, source_state));

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
        //         cerr << _haploGraph.get_id(handle) << " ";
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
        _haploGraph.follow_paths(cur_haplotype.second,
                                 [&](const gbwt::SearchState next_search) -> bool {
                                     next_searches.push_back(next_search);
                                     return true;
                                 });

        // if next_searches > 1, then we need to make multiple new haplotypes to be
        // recorded in haplotype_queue or one of the finished haplotype_handle_vectors.
        if (next_searches.size() > 1) {
            // for every next_search in next_searches, either create a new, extended
            // cur_haplotype to push into haplotype queue, or place in the
            // haplotypes_from_source_to_sink if haplotype extends to sink, or place in
            // the other_haplotypes if haplotype ends before reaching sink.
            for (gbwt::SearchState next_search : next_searches) {
                handle_t next_handle = _haploGraph.node_to_handle(next_search.node);
                // if (!snarl.has_node(snarl.get_id(next_handle)) &&
                // make_pair(_haploGraph.get_id(cur_haplotype.first.back()),_haploGraph.get_id(next_handle)))
                // {
                if (!_snarl.has_edge(cur_haplotype.first.back(), next_handle)) {
                    if (incorrect_connections.find(
                            _snarl.edge_handle(cur_haplotype.first.back(), next_handle)) ==
                        incorrect_connections.end()) {
                        cerr << "_snarl starting at node " << _source_id
                            << " and ending at " << _sink_id
                            << " has a thread that incorrectly connects two nodes that "
                               "don't have any edge connecting them. These two nodes are "
                            << _haploGraph.get_id(cur_haplotype.first.back()) << " and "
                            << _haploGraph.get_id(next_handle)
                            << ". This thread connection will be ignored." << endl;
                        incorrect_connections.emplace(
                            _snarl.edge_handle(cur_haplotype.first.back(), next_handle));

                        // todo: debug_statement
                        // cerr << "next handle(s) of handle "
                        //      << _snarl.get_id(cur_haplotype.first.back())
                        //      << " according to _snarl:" << endl;
                        // _snarl.follow_edges(cur_haplotype.first.back(), false,
                        //                    [&](const handle_t handle) {
                        //                        cerr << "\t" << _snarl.get_id(handle);
                        //                    });
                        // cerr << endl;
                    }
                    continue;
                }
                // copy over the vector<handle_t> of cur_haplotype:
                vector<handle_t> next_handle_vec(cur_haplotype.first);

                // add the new handle to the vec:
                next_handle_vec.push_back(next_handle);

                // if new_handle is the sink, put in haplotypes_from_source_to_sink
                if (_haploGraph.get_id(next_handle) == _sink_id) {
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
        // if next_handle is the sink, put in haplotypes_from_source_to_sink
        else if (_haploGraph.get_id(
                     _haploGraph.node_to_handle(next_searches.back().node)) == _sink_id) {
            // Then we need to add cur_haplotype + next_search to
            // haplotypes_from_source_to_sink.
            handle_t next_handle = _haploGraph.node_to_handle(next_searches.back().node);
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
            handle_t next_handle = _haploGraph.node_to_handle(next_searches.back().node);

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
        find_haplotypes_not_at_source(touched_handles);

    // move haplotypes_not_starting_at_source into other_haplotypes:
    other_haplotypes.reserve(other_haplotypes.size() +
                             haplotypes_not_starting_at_source.size());
    move(haplotypes_not_starting_at_source.begin(),
         haplotypes_not_starting_at_source.end(), back_inserter(other_haplotypes));

    //todo: debug_statement
    cerr << "lets look through all the haplotypes after extraction:" << endl;
    for (vector<handle_t> hap_vec : haplotypes_from_source_to_sink) {
        cerr << "new hap:" << endl;
        for (handle_t handle : hap_vec){
            cerr << _haploGraph.get_id(handle) << " " << _haploGraph.get_sequence(handle) << endl;
        }
    }

    return tuple<vector<vector<handle_t>>, vector<vector<handle_t>>,
                 unordered_set<handle_t>>{haplotypes_from_source_to_sink,
                                          other_haplotypes, touched_handles};
}


// Used to complete the traversal of a snarl along its haplotype threads, when there are
// handles connected to the snarl by threads that start after the source handle. (Threads
// that merely end before the sink handle are addressed in extract_gbwt_haplotypes).
// Arguments:
//      _haploGraph: the GBWTgraph containing the haplotype threads.
//      touched_handles: any handles found in the snarl so far.
//      sink_id: the id of the final handle in the snarl.
// Returns:
//      a vector of haplotypes in vector<handle_t> format that start in the middle of the
//      snarl.
vector<vector<handle_t>>
SnarlSequenceFinder::find_haplotypes_not_at_source(unordered_set<handle_t> &touched_handles) {
    // cerr << "find_haplotypes_not_at_source" << endl;

    /// Search every handle in touched handles for haplotypes starting at that point.
    // Any new haplotypes will be added to haplotype_queue.
    vector<pair<vector<handle_t>, gbwt::SearchState>> haplotype_queue;

    // Fully extended haplotypes (or haplotypes extended to the snarl's sink)
    // will be added to finished_haplotypes.
    vector<vector<handle_t>> finished_haplotypes;

    // In addition, we need to put the new handle into to_search, because a path may have
    // started on the new handle (which means we need to start a searchstate there.)
    unordered_set<handle_t> to_search;

    // We don't need to ever check the sink handle, since paths from the sink handle
    // extend beyond snarl.
    handle_t sink_handle = _haploGraph.get_handle(_sink_id);
    // touched_handles.erase(sink_handle);

    // Nested function for making a new_search. Identifies threads starting at a given
    // handle and
    //      either adds them as a full haplotype (if the haplotype is one handle long) or
    //      makes a new entry to haplotype_queue.
    auto make_new_search = [&](handle_t handle) {
        // Are there any new threads starting at this handle?
        gbwt::SearchState new_search =
            _haploGraph.index->prefix(_haploGraph.handle_to_node(handle));
        if (!new_search.empty()) {
            // Then add them to haplotype_queue.
            _haploGraph.follow_paths(
                new_search, [&](const gbwt::SearchState &next_search) -> bool {
                    handle_t next_handle = _haploGraph.node_to_handle(next_search.node);

                    /// check to make sure that the thread isn't already finished:
                    // if next_handle is the sink, or if this thread is only one handle
                    // long, then there isn't any useful string to extract from this.
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

                        // if next_handle hasn't been checked for starting threads, add to
                        // to_search.
                        if (touched_handles.find(next_handle) == touched_handles.end()) {
                            to_search.emplace(next_handle);
                        }
                    }
                    return true;
                });
        }
    };

    /// Extend any paths in haplotype_queue, and add any newly found handles to to_search.
    /// Then, check to see if there are any new threads on handles in to_search.
    /// Extend those threads, and add any newly found handles to to_search,
    /// then search for threads again in to_search again... repeat until to_search remains
    /// emptied of new handles.

    // for tracking whether the haplotype thread is still extending:
    bool still_extending;
    while (!to_search.empty() || !haplotype_queue.empty()) {
        while (!haplotype_queue.empty()) {
            // get a haplotype to extend out of haplotype_queue - a tuple of
            // (handles_traversed_so_far, last_touched_SearchState)
            pair<vector<handle_t>, gbwt::SearchState> cur_haplotype =
                haplotype_queue.back();
            haplotype_queue.pop_back();

            // get all the subsequent search_states that immediately follow the
            // searchstate from cur_haplotype.
            vector<gbwt::SearchState> next_searches;
            _haploGraph.follow_paths(cur_haplotype.second,
                                     [&](const gbwt::SearchState &next_search) -> bool {
                                         next_searches.push_back(next_search);
                                         return true;
                                     });

            for (gbwt::SearchState next_search : next_searches) {
                handle_t next_handle = _haploGraph.node_to_handle(next_search.node);

                // if next_search is empty, then we've fallen off the thread,
                // and cur_haplotype can be placed in finished_haplotypes as is for this
                // thread.
                if (next_search == gbwt::SearchState()) {
                    finished_haplotypes.push_back(cur_haplotype.first);
                }

                // if next_search is on the sink_handle,
                // then cur_haplotype.first + next_search goes to finished_haplotypes.
                else if (_haploGraph.get_id(next_handle) == _sink_id) {

                    // copy over the vector<handle_t> of cur_haplotype:
                    vector<handle_t> next_handle_vec(cur_haplotype.first);
                    // add next_handle
                    next_handle_vec.push_back(next_handle);
                    // place in finished_haplotypes
                    finished_haplotypes.push_back(next_handle_vec);

                    // also, if next_handle hasn't been checked for new threads, add to
                    // to_search.
                    if (touched_handles.find(next_handle) != touched_handles.end()) {
                        to_search.emplace(next_handle);
                    }

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

                    // also, if next_handle hasn't been checked for new threads, add to
                    // to_search.
                    if (touched_handles.find(next_handle) != touched_handles.end()) {
                        to_search.emplace(next_handle);
                    }
                }
            }
        }
        // Then, make more new_searches from the handles in to_search.
        for (handle_t handle : to_search) {
            make_new_search(handle); // will add to haplotype_queue if there's any
                                     // new_searches to be had.
        }
        to_search.clear();
    }
    return finished_haplotypes;
}

///////////////////////////////////
// exhaustive sequence extraction:
///////////////////////////////////
// Runs through the whole snarl and generates all possible strings representing walks
// from source to sink. Generates a combinatorial number of possible paths with splits
// in the snarl.
//todo: for consistency, have source_to_sink_exhaustive_path_finder return paths in format
//todo:     vector<vector<handle_t>> instead of vector<string>
pair<vector<string>, unordered_set<handle_t>>
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
    for (const handle_t &handle : lazier_topological_order(&_snarl)) {
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
    vector<string> walks;
    for (handle_t &sink : sinks) {
        for (string seq : sequences[sink]) {
            walks.push_back(seq);
        }
    }

    return make_pair(walks, touched_handles);
}


}
}
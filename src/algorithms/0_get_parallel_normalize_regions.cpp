#include "0_get_parallel_normalize_regions.hpp"


namespace vg {
namespace algorithms {
NormalizeRegionFinder::NormalizeRegionFinder(MutablePathDeletableHandleGraph &graph,
                                 const int max_region_size,
                                 const int max_snarl_spacing)
: _graph(graph), _max_region_size(max_region_size), _max_snarl_spacing(max_snarl_spacing) {}

/// @brief The main function to be called in NormalizeRegionFinder. Sets ups normalizable
/// regions that are completely independent of each other (no overlapping handles between
/// regions). This involves editing the graph so that two regions with a shared source &
/// sink have that source/sink split into two.
/// @param snarl_roots A vector of all the snarls in the graph that you want turned into
/// normalizable regions. (from the snarl manager, most likely)
/// @param parallel_normalize_regions An empty vector<pair<id_t, id_t>>, to be modified by
/// the function. This will be the regions that meet the NormalizeRegionFinder's
/// specifications.
/// @param nodes_to_remove An empty vector<id_t>, to be modified by the
/// function. Because some regions may have a shared source/sink that is of length one,
/// this function will duplicate that single node so that the leftmost of the two regions will
/// be empty. This needs to be deduplicated after the
/// parallelized normalization.
/// @return A tuple of arguments to pass to gbwt_update_items if you want an updated gbwt
/// to match the parallel normalize regions. Tuple: (_gbwt_graph, _gbwt_changelog, _gbwt)
std::vector<vg::RebuildJob::mapping_type> NormalizeRegionFinder::get_parallel_normalize_regions(const vector<const Snarl *> &snarl_roots, vector<pair<id_t, id_t>>& parallel_normalize_regions, set<id_t>& nodes_to_remove)
{
    vector<pair<id_t, id_t>> normalize_regions = get_normalize_regions(snarl_roots);
    vector<pair<id_t, id_t>> split_normalize_regions = split_sources_and_sinks(normalize_regions, nodes_to_remove);
    return _gbwt_changelog;
}


////////////////////////////////////////////////////////////
/// Functions called by get_parallel_normalize_regions:
////////////////////////////////////////////////////////////

vector<pair<id_t, id_t>> NormalizeRegionFinder::get_normalize_regions(const vector<const Snarl *> &snarl_roots) 
{
        //todo: uncomment below if I want? Should be equivalent, but slightly more efficient than the general case. I'm currently using it for testing.
        //if batch size is 1, just return the nontrivial roots from snarl_roots.
        // if (_max_region_size==1)
        // {
        //     return get_single_snarl_normalize_regions(snarl_roots);
        // }

    //otherwise, cluster snarls.
    vector<vector<const Snarl *> > snarl_clusters = cluster_snarls(snarl_roots);
    vector<pair<id_t, id_t>> regions = convert_snarl_clusters_to_regions(snarl_clusters);
    return regions;
}


/// @brief Checks to see if original source and sink of each snarl is shared with another 
///         snarl. If so, that source and/or sink is split into two. The ids of the 
///         original source and sink will always be the further of the two split nodes. 
///         (i.e., the neighboring snarl gets to keep the node with the original id, and 
///         the current snarl gets its node ids updated.)
/// @param normalize_regions 
/// @return a set of sources and sinks representing each (now isolated) snarl.
vector<pair<id_t, id_t>> NormalizeRegionFinder::split_sources_and_sinks(vector<pair<id_t, id_t>> normalize_regions, set<id_t>& nodes_to_remove){
    vector<pair<id_t, id_t>> new_normalize_regions;
    //todo: remove debug comment:
    // snarl 1883644 1883647 righttmost is of size 1
    // normalize_regions.clear();
    // pair<id_t, id_t> debug_region = make_pair(1883644, 1883647);
    // normalize_regions.push_back(debug_region);
    for (auto region : normalize_regions){
        cerr << "region: " << region.first << " " << region.second << endl;
        handle_t leftmost_handle = _graph.get_handle(region.first);
        handle_t rightmost_handle = _graph.get_handle(region.second);

        // look only to the right of rightmost_handle
        int right_of_rightmost = 0;
        int left_of_leftmost = 0;
        _graph.follow_edges(rightmost_handle, false, [&](const handle_t handle) 
        {
            right_of_rightmost += 1;
        });
        // look only to the left of leftmost_handle
        _graph.follow_edges(leftmost_handle, true, [&](const handle_t handle) 
        {
            left_of_leftmost += 1;
        });

        //the leftmost and rightmost node ids will only be updated if they need to be changed.
        id_t new_leftmost = region.first;
        id_t new_rightmost = region.second;
        if (left_of_leftmost >1)
        {
            int original_handle_seq_len = _graph.get_sequence(leftmost_handle).size();
            // if (_graph.get_sequence(leftmost_handle).size() == 1)
            // {

            //     // cerr << "snarl " << region.first << " " << region.second << " leftmost is of size 1" << endl;
            // }
            // divide handle always gives the original node id to the leftmost of the two 
            // new handles. In this case, that's what we want.
            pair<handle_t, handle_t> new_leftmosts = _graph.divide_handle(leftmost_handle, _graph.get_sequence(leftmost_handle).size()/2);
            cerr << "new leftmosts: 1: " << _graph.get_id(new_leftmosts.first) << " " << _graph.get_sequence(new_leftmosts.first) << " 2: " << _graph.get_id(new_leftmosts.second) << " " << _graph.get_sequence(new_leftmosts.second) << endl; 
            new_leftmost = _graph.get_id(new_leftmosts.second); 
            // if seq is length one, mark the empty handle as to-remove after the separated regions are done with.
            if (original_handle_seq_len == 1)
            {
                nodes_to_remove.emplace(_graph.get_id(new_leftmosts.first));
            }
            // encode the change in gbwt path in the gbwt.
            gbwt::vector_type original_gbwt_path;
            original_gbwt_path.emplace_back(gbwt::Node::encode(region.first, false));
            gbwt::vector_type new_gbwt_path;
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_leftmosts.first), false));
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_leftmosts.second), false));
            _gbwt_changelog.emplace_back(original_gbwt_path, new_gbwt_path);
            // cerr << "contents of _gbwt_changelog" << endl; 
            // for (auto region : _gbwt_changelog)
            // {
            //     cerr << "region original: ";
            //     for (auto node : region.first)
            //     {
            //         cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
            //     }
            //     cerr << endl;
            //     cerr << "new region: ";
            //     for (auto node : region.second)
            //     {
            //         cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
            //         // cerr << node << " ";
            //     }
            //     cerr << endl;
            // // cerr << "region: " << region.first.back() << " " << region.second.back() << endl; 
            // }

        }
        if (right_of_rightmost >1)
        {
            // cerr << "_gbwt_changelog ids: " << endl;;
            // for (auto region : _gbwt_changelog)
            // {
            //     cerr << "original node series:" << endl;
            //     for (auto node : region.first)
            //     {
            //         cerr << node ;
            //     }
            //     cerr << endl;
            //     cerr << "updated node series:" << endl;
            //     for (auto node : region.second)
            //     {
            //         cerr << node ;
            //     }
            //     cerr << endl;
            // }
            // cerr << "end _gbwt_changelog ids. " << endl;;

            int original_handle_seq_len = _graph.get_sequence(leftmost_handle).size();

            // cerr << "before adding changes:" << endl;
            // for (auto region : _gbwt_changelog)
            // {
            //     cerr << "region original: ";
            //     for (auto node : region.first)
            //     {
            //         cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
            //     }
            //     cerr << endl;
            //     cerr << "new region: ";
            //     for (auto node : region.second)
            //     {
            //         cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
            //         // cerr << node << " ";
            //     }
            //     cerr << endl;
            //     // cerr << "region: " << region.first.back() << " " << region.second.back() << endl; 
            // }
            // if (_graph.get_sequence(rightmost_handle).size() == 1)
            // {
            //     cerr << "snarl " << region.first << " " << region.second << " righttmost is of size 1" << endl;
            // }

            // cerr << "original node id: " << _graph.get_id(rightmost_handle) << endl;
            pair<handle_t, handle_t> new_rightmosts = _graph.divide_handle(rightmost_handle, _graph.get_sequence(rightmost_handle).size()/2);
            // cerr << "0" << endl;
            // cerr << "new rightmosts left id: " << _graph.get_id(new_rightmosts.first) << " " << _graph.get_sequence(new_rightmosts.first) << endl;
            // cerr << "new rightmosts right id: " << _graph.get_id(new_rightmosts.second) << " " << _graph.get_sequence(new_rightmosts.second) << endl;
            // gotta move the original node id to the rightmost of the divided handles, 
            // rather than the leftmost:
            handle_t dummy_handle = _graph.create_handle("A");
            id_t new_node_id = _graph.get_id(dummy_handle);
            _graph.destroy_handle(dummy_handle);
            overwrite_node_id(_graph.get_id(new_rightmosts.first), new_node_id);
            overwrite_node_id(_graph.get_id(new_rightmosts.second), region.second);
            // cerr << "1" << endl;
            // cerr << "new rightmosts left id: " << _graph.get_id(_graph.get_handle(new_node_id)) << " " << _graph.get_sequence(_graph.get_handle(new_node_id)) << endl;
            // cerr << "new rightmosts right id: " << _graph.get_id(_graph.get_handle(region.second)) << " " << _graph.get_sequence(_graph.get_handle(region.second)) << endl;
            new_rightmost = new_node_id; 
            // if seq is length one, mark the empty handle as to-remove after the separated regions are done with.
            if (original_handle_seq_len == 1)
            {
                nodes_to_remove.emplace(new_node_id);
            }
            // encode the change in gbwt path in the gbwt.
            gbwt::vector_type original_gbwt_path;
            original_gbwt_path.emplace_back(gbwt::Node::encode(region.second, false));
            gbwt::vector_type new_gbwt_path;
            // cerr << "_graph.get_id(new_rightmosts.second): " << _graph.get_id(new_rightmosts.second) << endl;
            // cerr << "gbwt encoding of _graph.get_id(new_rightmosts.second): " << gbwt::Node::encode(_graph.get_id(new_rightmosts.second), false) << endl;
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_rightmosts.first), false));
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_rightmosts.second), false));
            // cerr << "original_gbwt_path: " << endl;
        //     for (auto id : original_gbwt_path)
        //     {
        //         cerr << id << " ";
        //     }
        //     cerr << endl;
        //     cerr << "end original_gbwt_path: " << endl;

        //     cerr << "new_gbwt_path: " << endl;
        //     for (auto id : new_gbwt_path)
        //     {
        //         cerr << id << " ";
        //     }
        //     cerr << endl;
        //     cerr << "end new_gbwt_path: " << endl;
        //     _gbwt_changelog.emplace_back(original_gbwt_path, new_gbwt_path);
            
        //     cerr << "_gbwt_changelog ids: " << endl;;
        //     for (auto region : _gbwt_changelog)
        //     {
        //         cerr << "original node series:" << endl;
        //         for (auto node : region.first)
        //         {
        //             cerr << node << " ";
        //         }
        //         cerr << endl;
        //         cerr << "updated node series:" << endl;
        //         for (auto node : region.second)
        //         {
        //             cerr << node << " ";
        //         }
        //         cerr << endl;
        //     }
        //     cerr << "end _gbwt_changelog ids. " << endl;;


        //     cerr << "contents of _gbwt_changelog" << endl; 
        //     for (auto region : _gbwt_changelog)
        //     {
        //         cerr << "region original: ";
        //         for (auto node : region.first)
        //         {
        //             cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
        //         }
        //         cerr << endl;
        //         cerr << "new region: ";
        //         for (auto node : region.second)
        //         {
        //             cerr << node << " "  << _graph.get_sequence(_graph.get_handle(node)) << endl;
        //             // cerr << node << " ";
        //         }
        //         cerr << endl;
        //         // cerr << "region: " << region.first.back() << " " << region.second.back() << endl; 
        //     }
        }


        // cerr << "here is the sequences neighboring and in the new leftmost and rightmost " << endl;
        // cerr << "new rightmost id: " << new_rightmost << endl;
        // cerr << "seq in new_rightmost " << _graph.get_sequence(_graph.get_handle(new_rightmost)) << endl; 
        // // look around the rightmost_handle
        // cerr << "left of rightmost" << endl;
        // _graph.follow_edges(_graph.get_handle(new_rightmost), true, [&](const handle_t handle) 
        // {
        //     cerr << _graph.get_id(handle) << _graph.get_sequence(handle) << endl;
        // });
        // cerr << "right of rightmost" << endl;
        // _graph.follow_edges(_graph.get_handle(new_rightmost), false, [&](const handle_t handle) 
        // {
        //     cerr << _graph.get_id(handle) << _graph.get_sequence(handle) << endl;
        // });
        // cerr << "seq in new_leftmost " << _graph.get_sequence(_graph.get_handle(new_leftmost)) << endl; 
        // cerr << "new leftmost id: " << new_rightmost << endl;
        // // look around the leftmost_handle
        // cerr << "left of leftmost" << endl;
        // _graph.follow_edges(_graph.get_handle(new_leftmost), true, [&](const handle_t handle) 
        // {
        //     cerr << _graph.get_id(handle) << _graph.get_sequence(handle) << endl;
        // });
        // cerr << "right of leftmost" << endl;
        // _graph.follow_edges(_graph.get_handle(new_leftmost), false, [&](const handle_t handle) 
        // {
        //     cerr << _graph.get_id(handle) << _graph.get_sequence(handle) << endl;
        // });

        new_normalize_regions.push_back(make_pair(new_leftmost, new_rightmost));

    }
    return new_normalize_regions;
}

////////////////////////////////////////////////////////////
/// Functions called by get_normalize_regions:
////////////////////////////////////////////////////////////


vector<vector<const Snarl *> > NormalizeRegionFinder::cluster_snarls(const vector<const Snarl *> &snarl_roots) {    
    vector<vector<const Snarl *> > snarl_clusters;
    auto cur_snarl = snarl_roots.begin();
    int trivial_count = 0;
    vector<const Snarl*> first_cluster;
    snarl_clusters.push_back(first_cluster);

    while (cur_snarl != snarl_roots.end())
    {
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

        // If we are considering extending a snarl cluster, first make sure that we 
        // haven't reached the end conditions of the current cluster.
        // batch size exceeded? Or snarls aren't part of the same connected component?
        // start new cluster, and trim the previous one of any trailing trivial snarls.
        if (snarl_clusters.back().size() != 0)
        {
            const Snarl prev_snarl = *snarl_clusters.back().back();
            if (snarl_clusters.back().size() == _max_region_size || !snarls_adjacent(prev_snarl, **cur_snarl) || trivial_count > _max_snarl_spacing || cyclic)
            {
                //If we're here, we've reached the end condition for this cluster.
                //Trim the tirivial snarls from the end of the cluster:
                for (int i = 0; i < trivial_count; i++)
                {
                    snarl_clusters.back().pop_back();
                }
                trivial_count = 0;

                //start new cluster
                vector<const Snarl*> new_cluster;
                snarl_clusters.push_back(new_cluster);
            }

        }
        
        // special case: if we find a cyclic snarl and its the first snarl of a new cluster:
        if (snarl_clusters.back().size() == 0 && cyclic)
        {
            // Then just skip the cyclic snarl.
            cur_snarl++;
            continue;
        }
        
        bool trivial = is_trivial(**cur_snarl);
        // if we have a snarl that isn't trivial, add it to the latest cluster.
        if (!trivial)
        {
            snarl_clusters.back().push_back(*cur_snarl);
            //reset the trivial_count
            trivial_count=0;
        }
            
        // if we have a snarl that is trivial, but we have a cluster in progress, add it to the cluster.
        else if (snarl_clusters.back().size() > 0)
        {
            snarl_clusters.back().push_back(*cur_snarl);
            trivial_count++;
        }
        // if we have a snarl that is trivial, and we have no cluster in progress, skip it. (no code needed.)
        cur_snarl++;
    }
    //check to see if the snarl_cluster at the end is empty. If so, discard it.
    if (snarl_clusters.back().size() == 0)
    {
        snarl_clusters.pop_back();
    }
    return snarl_clusters;
}


vector<pair<id_t, id_t>> NormalizeRegionFinder::convert_snarl_clusters_to_regions(const vector<vector<const Snarl *> >& clusters) {
    vector<pair<id_t, id_t>> normalize_regions;
    for (auto cluster : clusters )
    {
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
            }
            else if (cur_region.first == cluster[i]->start().node_id())
            {
                cur_region.first = cluster[i]->end().node_id();
            }
            else if (cur_region.second == cluster[i]->end().node_id())
            {
                cur_region.second = cluster[i]->start().node_id();
            }
            else if (cur_region.second == cluster[i]->start().node_id())
            {
                cur_region.second = cluster[i]->end().node_id();
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
        normalize_regions.push_back(cur_region);

    }
    return normalize_regions;
}

////////////////////////////////////////////////////////////
/// Functions called by cluster_snarls:
////////////////////////////////////////////////////////////

// TODO: change the arguments to handles, which contain orientation within themselves.
// Given a start and end node id, construct an extract subgraph between the two nodes
// (inclusive). Arguments:
//      graph: a pathhandlegraph containing the snarl with embedded paths.
//      source_id: the source of the snarl of interest.
//      sink_id: the sink of the snarl of interest.
// Returns:
//      a SubHandleGraph containing only the handles in _graph that are between start_id
//      and sink_id.
SubHandleGraph NormalizeRegionFinder::extract_subgraph(const HandleGraph &graph,
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

//snarls_adjacent used to identify if two snarls overlap at one of their boundary nodes.
bool NormalizeRegionFinder::snarls_adjacent(const Snarl& snarl_1, const Snarl& snarl_2) 
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

bool NormalizeRegionFinder::is_trivial(const Snarl& snarl) {
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

////////////////////////////////////////////////////////////
/// Functions called by split_sources_and_sinks:
////////////////////////////////////////////////////////////


/**
 * Deletes the given handle's underlying node, and returns a new handle to a new node 
 * with the desired node_id (copied from normalize_snarls)//todo: move to an appropiate 
 *                                                          todo: place in handle library?
 * 
 * @param  {id_t} handle     : The old node id, to be replaced with a new node id.
 * @param  {id_t} node_id    : The node id for the new node. Cannot be currently in use in
 *                              the graph.
 * @return {handle_t}        : The new handle, in the same position as the original handle
 *                              in the graph, but with the new node_id.
 */
handle_t NormalizeRegionFinder::overwrite_node_id(const id_t old_node_id, const id_t new_node_id)
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

// /// @brief Specifically for replacing node six in tiny's example with a one base character for testing purposes.
// std::vector<vg::RebuildJob::mapping_type> NormalizeRegionFinder::debug_replace_node_six()
// {
//     handle_t old_six = _graph.get_handle(6);
//     handle_t new_six = _graph.create_handle("A");
//     _graph.destroy_handle(old_six);
//     overwrite_node_id(_graph.get_id(new_six), 6);
    
//     return _gbwt_changelog;
// }
}
}

#include "0_get_parallel_normalize_regions.hpp"
// #include <bdsg/snarl_distance_index.hpp>


namespace vg {
namespace algorithms {
NormalizeRegionFinder::NormalizeRegionFinder(MutablePathDeletableHandleGraph &graph,
                                 const int max_region_size,
                                 const int max_region_gap) //todo: remove snarl spacing, or adapt it for when I want to divide clusters, esp once abpoa allows for larger clusters.
: _graph(graph), _max_region_size(max_region_size), _max_region_gap(max_region_gap){}

/// @brief The main function to be called in NormalizeRegionFinder. Sets ups normalizable
/// regions that are completely independent of each other (no overlapping handles between
/// regions). This involves editing the graph so that two regions with a shared source &
/// sink have that source/sink split into two.
/// @param snarl_roots A vector of all the snarls in the graph that you want turned into
/// normalizable regions. (from the snarl manager, most likely)
/// @param parallel_normalize_regions An empty vector<pair<id_t, id_t>>, to be modified by
/// the function. This will be the regions that meet the NormalizeRegionFinder's
/// specifications.
/// @param desegregation_candidates a vector of pairs. Each pair's first item is the two
/// new_nodes created for the parallelization process. The second item is the original
/// node id, which will be reinstated after normalization using desegregated_nodes (via a
/// fresh NormalizeRegionFinder object).
/// @param segregated_node_to_parent this is a map for tracking which of the new_node ids
/// correspond to which original_node id. This allows normalize to record the update_gbwt
/// process without ever actually mentioning any of the segregated_nodes, since we'll be
/// removing them from the graph via fxn desegregated_nodes before updating the gbwt for
/// the 2nd and final time.
/// @return A tuple of arguments to pass to gbwt_update_items if you want an updated gbwt
/// to match the parallel normalize regions. Tuple: (_gbwt_graph, _gbwt_changelog, _gbwt)
std::vector<vg::RebuildJob::mapping_type> NormalizeRegionFinder::get_parallel_normalize_regions(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots, const SnarlDistanceIndex& distance_index, vector<pair<id_t, id_t>>& parallel_normalize_regions, vector< pair< pair< id_t, id_t >, id_t > >& desegregation_candidates, unordered_map<id_t, id_t>& segregated_node_to_parent)
{
    vector<pair<vg::id_t, vg::id_t>> clustered_snarls = cluster_snarls(snarl_roots, distance_index);
    // vector<pair<id_t, id_t>> clustered_regions = convert_snarl_clusters_to_regions(clustered_snarls);

    // parallel_normalize_regions = split_sources_and_sinks(snarl_roots, desegregation_candidates);
    parallel_normalize_regions = split_sources_and_sinks(clustered_snarls, desegregation_candidates, segregated_node_to_parent);

    cerr << "clusters were reset " << _reset_cluster_because_too_big << " times." << endl;

    // cerr << "segregated_node_to_parent output: " << endl;
    // for (auto item : segregated_node_to_parent)
    // {
    //     cerr << item.first << " " << item.second << endl;
    // }

    return _gbwt_changelog;
}



////////////////////////////////////////////////////////////
/// Functions called by get_parallel_normalize_regions:
////////////////////////////////////////////////////////////

// vector<pair<id_t, id_t>> NormalizeRegionFinder::get_normalize_regions(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots) 
// {
//         //todo: uncomment below if I want? Should be equivalent, but slightly more efficient than the general case. I'm currently using it for testing.
//         //if batch size is 1, just return the nontrivial roots from snarl_roots.
//         // if (_max_region_size==1)
//         // {
//         //     return get_single_snarl_normalize_regions(snarl_roots);
//         // }

//     //otherwise, cluster snarls.
//     vector<vector<pair<vg::id_t, vg::id_t>> > snarl_clusters = cluster_snarls(snarl_roots);
//     // cerr << "snarl clusters: " << endl;
//     // for (auto cluster : snarl_clusters)
//     // {
//     //     cerr << "one cluster: " << endl;
//     //     for (auto root : cluster)
//     //     {
//     //         cerr << root.first << " " << root.second << endl;
//     //     }
//     // }
//     vector<pair<id_t, id_t>> regions = convert_snarl_clusters_to_regions(snarl_clusters); //todo: MAKE THIS NOT BROKEN! right now, it assumes the snarls are ordered in the cluster from left to right.
//     return regions;
// }


/// @brief Checks to see if original source and sink of each snarl is shared with another 
///         snarl. If so, that source and/or sink is split into two. The ids of the 
///         original source and sink will always be the further of the two split nodes. 
///         (i.e., the neighboring snarl gets to keep the node with the original id, and 
///         the current snarl gets its node ids updated.)
/// @param normalize_regions 
/// @param desegregation_candidates a vector of pairs. The first element of the pair is a
/// pair of the two graph node ids created by handle.split(). The second element in the pair is the
/// id_t that was the original graph node id of the handle before splitting. This wil lbe used for
/// de-splitting the snarl later.
/// @return a set of sources and sinks representing each (now isolated) snarl.
vector<pair<id_t, id_t>> NormalizeRegionFinder::split_sources_and_sinks(vector<pair<id_t, id_t>> normalize_regions, vector< pair< pair< id_t, id_t >, id_t > >& desegregation_candidates, unordered_map<id_t, id_t>& segregated_node_to_parent){
    vector<pair<id_t, id_t>> new_normalize_regions;

    //todo: remove debug comment:
    // snarl 1883644 1883647 righttmost is of size 1
    // normalize_regions.clear();
    // pair<id_t, id_t> debug_region = make_pair(1883644, 1883647);
    // normalize_regions.push_back(debug_region);
    for (auto region : normalize_regions){
        // cerr << "region: " << region.first << " " << region.second << endl;
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
            id_t original_leftmost = _graph.get_id(leftmost_handle);
            //get context of original leftmost for _gbwt_changelog:
            id_t original_left_context;
            _graph.follow_edges(leftmost_handle, true, [&](const handle_t handle) 
            {
                original_left_context = _graph.get_id(handle);
                return false;
            });
            id_t original_right_context;
            _graph.follow_edges(leftmost_handle, false, [&](const handle_t handle) 
            {
                original_right_context = _graph.get_id(handle);
                return false;
            });

            //todo: begin debug_code
            std::cerr << "looking at original leftmost_handle: " << _graph.get_id(leftmost_handle) << endl; 
            vector<step_handle_t> steps_0 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(leftmost_handle)));
            for (step_handle_t step : steps_0)
            {
                std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            }
            //todo: end debug_code

            pair<handle_t, handle_t> new_leftmosts = _graph.divide_handle(leftmost_handle, _graph.get_sequence(leftmost_handle).size()/2);
            //todo: begin debug_code
            std::cerr << "looking at handle new_leftmosts.first: " << _graph.get_id(new_leftmosts.first) << endl; 
            vector<step_handle_t> steps_1 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(new_leftmosts.first)));
            for (step_handle_t step : steps_1)
            {
                std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            }
            std::cerr << "looking at handle new_leftmosts.second: " << _graph.get_id(new_leftmosts.second) << endl; 
            vector<step_handle_t> steps_2 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(new_leftmosts.second)));
            for (step_handle_t step : steps_2)
            {
                std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            }

            //todo: end debug_code

            
            desegregation_candidates.push_back(make_pair(make_pair(_graph.get_id(new_leftmosts.first), _graph.get_id(new_leftmosts.second)), original_leftmost));

            if (_graph.get_id(new_leftmosts.first) != original_leftmost)
            {
                segregated_node_to_parent[_graph.get_id(new_leftmosts.first)] = original_leftmost;
            }
            if (_graph.get_id(new_leftmosts.second) != original_leftmost)
            {
                segregated_node_to_parent[_graph.get_id(new_leftmosts.second)] = original_leftmost;
            }
            // cerr << "new leftmosts: 1: " << _graph.get_id(new_leftmosts.first) << " " << _graph.get_sequence(new_leftmosts.first) << " 2: " << _graph.get_id(new_leftmosts.second) << " " << _graph.get_sequence(new_leftmosts.second) << endl; 
            new_leftmost = _graph.get_id(new_leftmosts.second); 
            //// if seq is length one, mark the empty handle as to-remove after the separated regions are done with.
            // if (original_handle_seq_len == 1)
            // {
            //     nodes_to_remove.emplace(_graph.get_id(new_leftmosts.first));
            // }
            // encode the change in gbwt path in the gbwt.
            gbwt::vector_type original_gbwt_path;
            // original_gbwt_path.emplace_back(gbwt::Node::encode(original_left_context, false));
            original_gbwt_path.emplace_back(gbwt::Node::encode(region.first, false));
            // original_gbwt_path.emplace_back(gbwt::Node::encode(original_right_context, false));
            gbwt::vector_type new_gbwt_path;
            // new_gbwt_path.emplace_back(gbwt::Node::encode(original_left_context, false));
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_leftmosts.first), false));
            new_gbwt_path.emplace_back(gbwt::Node::encode(_graph.get_id(new_leftmosts.second), false));
            // new_gbwt_path.emplace_back(gbwt::Node::encode(original_right_context, false));

            // cerr << "replacing original leftmost " << region.first << " with " << _graph.get_id(new_leftmosts.first) << " " << _graph.get_id(new_leftmosts.second) << endl;
            _gbwt_changelog.emplace_back(original_gbwt_path, new_gbwt_path);
            // cerr << "contents of _gbwt_changelog" << endl; 
            // for (auto region : _gbwt_changelog)
            // {
            //     cerr << "region original: ";
            //     for (auto node : region.first)
            //     {
            //         cerr << node << " "; //<< _graph.get_sequence(_graph.get_handle(node)) << endl;
            //     }
            //     cerr << endl;
            //     cerr << "new region: ";
            //     for (auto node : region.second)
            //     {
            //         cerr << node << " "; //<< _graph.get_sequence(_graph.get_handle(node)) << endl;
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

            id_t original_rightmost = _graph.get_id(rightmost_handle);
            //todo: begin debug_code
            // if (original_rightmost==170311)
            // {

            //     std::cerr << "looking at original rightmost_handle: " << _graph.get_id(rightmost_handle) << endl; 
            //     vector<step_handle_t> steps_0 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(rightmost_handle)));
            //     for (step_handle_t step : steps_0)
            //     {
            //         std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            //     }
            // }
            //todo: end debug_code
            pair<handle_t, handle_t> new_rightmosts = _graph.divide_handle(rightmost_handle, _graph.get_sequence(rightmost_handle).size()/2);

            //todo: begin debug_code
            // if (original_rightmost==170311)
            // {

            //     std::cerr << "looking at handle new_rightmosts.first: " << _graph.get_id(new_rightmosts.first) << endl; 
            //     vector<step_handle_t> steps_1 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(new_rightmosts.first)));
            //     for (step_handle_t step : steps_1)
            //     {
            //         std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            //     }
            //     std::cerr << "looking at handle new_rightmosts.second: " << _graph.get_id(new_rightmosts.second) << endl; 
            //     vector<step_handle_t> steps_2 = _graph.steps_of_handle(_graph.get_handle(_graph.get_id(new_rightmosts.second)));
            //     for (step_handle_t step : steps_2)
            //     {
            //         std::cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " on " << _graph.get_id(_graph.get_handle_of_step(step)) << endl; 
            //     }
            // }
            //todo: end debug_code

            // gotta move the original node id to the rightmost of the divided handles, 
            // rather than the leftmost:
            handle_t dummy_handle = _graph.create_handle("A");
            id_t new_node_id = _graph.get_id(dummy_handle);
            _graph.destroy_handle(dummy_handle);

            overwrite_node_id(_graph.get_id(new_rightmosts.first), new_node_id);
            overwrite_node_id(_graph.get_id(new_rightmosts.second), region.second);
            desegregation_candidates.push_back(make_pair(make_pair(new_node_id, region.second), original_rightmost));
            if (_graph.get_id(new_rightmosts.first) != original_rightmost)
            {
                segregated_node_to_parent[_graph.get_id(new_rightmosts.first)] = original_rightmost;
            }
            if (_graph.get_id(new_rightmosts.second) != original_rightmost)
            {
                segregated_node_to_parent[_graph.get_id(new_rightmosts.second)] = original_rightmost;
            }

            

            // cerr << "1" << endl;
            // cerr << "new rightmosts left id: " << _graph.get_id(_graph.get_handle(new_node_id)) << " " << _graph.get_sequence(_graph.get_handle(new_node_id)) << endl;
            // cerr << "new rightmosts right id: " << _graph.get_id(_graph.get_handle(region.second)) << " " << _graph.get_sequence(_graph.get_handle(region.second)) << endl;
            new_rightmost = new_node_id; 
            // if seq is length one, mark the empty handle as to-remove after the separated regions are done with.
            // if (original_handle_seq_len == 1)
            // {
            //     nodes_to_remove.emplace(new_node_id);
            // }
            // encode the change in gbwt path in the gbwt.
            gbwt::vector_type original_gbwt_path;
            original_gbwt_path.emplace_back(gbwt::Node::encode(region.second, false));
            gbwt::vector_type new_gbwt_path;
            // cerr << "_graph.get_id(new_rightmosts.second): " << _graph.get_id(new_rightmosts.second) << endl;
            // cerr << "gbwt encoding of _graph.get_id(new_rightmosts.second): " << gbwt::Node::encode(_graph.get_id(new_rightmosts.second), false) << endl;
            new_gbwt_path.emplace_back(gbwt::Node::encode(new_node_id, false));
            new_gbwt_path.emplace_back(gbwt::Node::encode(region.second, false));
            _gbwt_changelog.emplace_back(original_gbwt_path, new_gbwt_path);

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
        // break;

    }

    return new_normalize_regions;
}

// ////////////////////////////////////////////////////////////
// /// Functions called by get_normalize_regions:
// ////////////////////////////////////////////////////////////


// vector< vector<pair<vg::id_t, vg::id_t>> > NormalizeRegionFinder::cluster_snarls(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots) {    
//     vector< vector<pair<vg::id_t, vg::id_t>> > snarl_clusters;
//     auto cur_snarl = snarl_roots.begin();
//     int trivial_count = 0;
//     vector<pair<vg::id_t, vg::id_t>> first_cluster;
//     snarl_clusters.push_back(first_cluster);

//     while (cur_snarl != snarl_roots.end())
//     {
//         //todo: Find out if this graph extraction is highly inefficient. If so, remove the graph_extraction I perform here, or else integrate it into the 
//         //todo:     overall normalizer in a way that prevents me from having to extract it twice (if that's efficient).
        
//         SubHandleGraph snarl_graph = extract_subgraph(_graph, cur_snarl->first, cur_snarl->second);
//         bool cyclic = !handlealgs::is_acyclic(&snarl_graph);

//         // If we are considering extending a snarl cluster, first make sure that we 
//         // haven't reached the end conditions of the current cluster.
//         // batch size exceeded? Or snarls aren't part of the same connected component?
//         // start new cluster, and trim the previous one of any trailing trivial snarls.
//         if (snarl_clusters.back().size() != 0)
//         {
//             pair<id_t, id_t> prev_snarl = snarl_clusters.back().back();
//             if (snarl_clusters.back().size() == _max_region_size || !snarls_adjacent(prev_snarl, *cur_snarl) || trivial_count > _max_snarl_spacing || cyclic)
//             {
//                 //If we're here, we've reached the end condition for this cluster.
//                 //Trim the trivial snarls from the end of the cluster:
//                 for (int i = 0; i < trivial_count; i++)
//                 {
//                     snarl_clusters.back().pop_back();
//                 }
//                 trivial_count = 0;

//                 //start new cluster
//                 vector<pair<id_t, id_t>> new_cluster;
//                 snarl_clusters.push_back(new_cluster);
//             }

//         }
        
//         // special case: if we find a cyclic snarl and its the first snarl of a new cluster:
//         if (snarl_clusters.back().size() == 0 && cyclic)
//         {
//             // Then just skip the cyclic snarl.
//             cur_snarl++;
//             continue;
//         }
        
//         bool trivial = is_trivial(*cur_snarl);
//         // if we have a snarl that isn't trivial, add it to the latest cluster.
//         if (!trivial)
//         {
//             snarl_clusters.back().push_back(*cur_snarl);
//             //reset the trivial_count
//             trivial_count=0;
//         }
            
//         // if we have a snarl that is trivial, but we have a cluster in progress, add it to the cluster.
//         else if (snarl_clusters.back().size() > 0)
//         {
//             snarl_clusters.back().push_back(*cur_snarl);
//             trivial_count++;
//         }
//         // if we have a snarl that is trivial, and we have no cluster in progress, skip it. (no code needed.)
//         cur_snarl++;
//     }
//     //check to see if the snarl_cluster at the end is empty. If so, discard it.
//     if (snarl_clusters.back().size() == 0)
//     {
//         snarl_clusters.pop_back();
//     }
//     return snarl_clusters;
// }

// //todo: instead of either of these problematic hacky-ways of turning (poor) clusters to regions (which in the most recent implementation doesn't take into account snarls put in the lcuster in reverse order), try doing something based on the distance_index itself.
// vector<pair<id_t, id_t>> NormalizeRegionFinder::convert_snarl_clusters_to_regions(const vector<vector<pair<vg::id_t, vg::id_t>> >& clusters) {
//     vector<pair<id_t, id_t>> normalize_regions;
//     for (auto cluster : clusters)
//     {
//         pair<id_t, id_t> normalize_region;
//         normalize_region.first = cluster.front().first;
//         normalize_region.second = cluster.back().second;
//         normalize_regions.push_back(normalize_region);
//     }
//     return normalize_regions;

    // for (auto cluster : clusters )
    // {
    //     pair<id_t, id_t> cur_region = cluster.front(); 
    //     // if (cluster.front()->start().backward())
    //     // {
    //     //     cur_region.first = cluster.front()->end().node_id();
    //     //     cur_region.second = cluster.front()->start().node_id();
    //     // } 
    //     // else
    //     // {
    //     //     cur_region.first = cluster.front()->start().node_id();
    //     //     cur_region.second = cluster.front()->end().node_id();
    //     // }
    //     for (int i = 1; i != cluster.size(); i++)
    //     {
    //         //find which part of the cur_region to replace with this snarl:
    //         if (cur_region.first == cluster[i].second())
    //         {
    //             cur_region.first = cluster[i].first();
    //         }
    //         else if (cur_region.first == cluster[i].first())
    //         {
    //             cur_region.first = cluster[i].second();
    //         }
    //         else if (cur_region.second == cluster[i].second())
    //         {
    //             cur_region.second = cluster[i].first();
    //         }
    //         else if (cur_region.second == cluster[i].first())
    //         {
    //             cur_region.second = cluster[i].second();
    //         }
    //         else
    //         {
    //             // the current snarl is not adjacent to the previous snarl; we have left 
    //             // the snarl chain. 
    //             // This is supposed to be guaranteed not to happen, because of the way 
    //             // get_normalize_regions generates clusters. Raise an error.
    //             cerr << "error:[vg normalize] snarl cluster is not within only a single"<< 
    //             " connected component. There is likely a bug in whatever passed clusters"<<
    //             " to convert_snarl_clusters_to_regions." << endl;
    //             exit(1);
    //         }
            
    //     }
    //     normalize_regions.push_back(cur_region);

    // }
// }

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

// //snarls_adjacent used to identify if two snarls overlap at one of their boundary nodes.
// bool NormalizeRegionFinder::snarls_adjacent(const pair<id_t, id_t>& snarl_1, const pair<id_t, id_t>& snarl_2) 
// {
//     return (snarl_1.second == snarl_2.first || snarl_2.second == snarl_1.first);
// }


// bool NormalizeRegionFinder::is_trivial(const pair<id_t, id_t>& snarl) {
//     vector<id_t> nodes_adjacent_to_source;
//     handle_t first_h = _graph.get_handle(snarl.first);

//     _graph.follow_edges(first_h, false, [&](handle_t next)
//     {
//         nodes_adjacent_to_source.push_back(_graph.get_id(next));
//     });

//     if (nodes_adjacent_to_source.size() == 1)
//     {
//         if (nodes_adjacent_to_source.back() == snarl.second)
//         {
//             return true;
//         }
//         else
//         {
//             cerr << "error:[vg normalize] snarl with start" << snarl.first;
//             cerr << " and end " << snarl.second;
//             cerr << " is performing unexpectedly. The start has only one internal-to-snarl";
//             cerr << " neighbor, and it is not the end node. Malformed snarl input?" << endl;
//             exit(1);
//         }
//     }
//     else
//     {
//         return false;
//     }
// }

////////////////////////////////////////////////////////////
/// Functions called by split_sources_and_sinks:
////////////////////////////////////////////////////////////


/**
 * Copied from SnarlNormalizer.
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

/// @brief Function used to remove the splits introduced in split_sources_and_sinks.
///         (NOTE: this will overwrite the current _gbwt_changelog.
/// @param desegregation_candidates every candidate in desegregation_candidates has a pair
///         of ids that need to be merged into one id. That id is the other id in the outer pair.
///         i.e. pair<to_merge_pair, original_id>.
/// @return whatever needs to be updated in the gbwt.
std::vector<vg::RebuildJob::mapping_type> NormalizeRegionFinder::desegregate_nodes(vector<pair<pair<id_t, id_t>, id_t>> desegregation_candidates)
{
    _gbwt_changelog.clear();
    for (pair<pair<id_t, id_t>, id_t> candidate : desegregation_candidates)
    {
        pair<id_t, id_t> to_merge_pair = candidate.first;
        id_t original_id = candidate.second;
        //get edges pointing to the left of the to_merge_pair
        vector<handle_t> left_handles;
        _graph.follow_edges(_graph.get_handle(to_merge_pair.first), true, [&](handle_t left_handle){
            left_handles.push_back(left_handle);
            // cerr << "tracking edge between left_handle " << _graph.get_id(left_handle) << " and current id " << original_id << endl;
        });

        //get edges pointing to the right of the to_merge_pair.
        vector<handle_t> right_handles;
        _graph.follow_edges(_graph.get_handle(to_merge_pair.second), false, [&](handle_t right_handle){
            right_handles.push_back(right_handle);
        });

        //get combined sequence of the to_merge_pair.
        string to_merge_seq = _graph.get_sequence(_graph.get_handle(to_merge_pair.first)) + _graph.get_sequence(_graph.get_handle(to_merge_pair.second));

        //get gbwt nodes of the to_merge_pair.
        gbwt::vector_type before_merge;
        before_merge.push_back(gbwt::Node::encode(to_merge_pair.first, false));
        before_merge.push_back(gbwt::Node::encode(to_merge_pair.second, false));

        //move paths moving over the to_merge_pair to a temporary handle:
        handle_t temp_handle = _graph.create_handle(to_merge_seq);
        for (handle_t left : left_handles)
        {
            _graph.create_edge(left, temp_handle);
        }
        for (handle_t right : right_handles)
        {
            _graph.create_edge(temp_handle, right);
        }
        _graph.for_each_step_on_handle(_graph.get_handle(to_merge_pair.first), [&](step_handle_t step)
        {
            path_handle_t path = _graph.get_path_handle_of_step(step);
            if (!(step == _graph.path_end(path) || _graph.get_next_step(step) == _graph.path_end(path) || step == _graph.path_begin(path)))
            {
                //only move the path if it doesn't end partway through the node, because otherwise there's no way to anchor the end of the path to the new, bigger snarl.
                step_handle_t prev = _graph.get_previous_step(step);
                step_handle_t next = _graph.get_next_step(_graph.get_next_step(step));
                vector<handle_t> new_path_location;
                new_path_location.push_back(_graph.get_handle_of_step(prev));
                new_path_location.push_back(temp_handle);

                _graph.rewrite_segment(prev, next, new_path_location);
            }
        });

        //delete the to_merge_pair. (note: this should delete a node with the original_id.)
        //todo: comment out this if check if I'm confident it functions as expected.
        if (!(to_merge_pair.first == original_id || to_merge_pair.second == original_id))
        {
            cerr << "original id not in to_merge_pair." << endl;
        }
        _graph.destroy_handle(_graph.get_handle(to_merge_pair.first));
        _graph.destroy_handle(_graph.get_handle(to_merge_pair.second));

        //create a new node with id original_id. And give it the sequence and edges of the to_merge_pair.
        overwrite_node_id(_graph.get_id(temp_handle), original_id);

        //add gbwt nodes of the merged pair
        gbwt::vector_type after_merge;
        after_merge.push_back(gbwt::Node::encode(original_id, false));

        //record the graph edit in the _gbwt_changelog.
        _gbwt_changelog.push_back(make_pair(before_merge, after_merge));
    }

    return _gbwt_changelog;

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



vector<pair<id_t, id_t>> NormalizeRegionFinder::cluster_snarls(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots, const SnarlDistanceIndex& distance_index)
{
    vector<pair<id_t, id_t>> clustered_snarls;
    vector<pair<id_t, id_t>> skipped_snarls;
    int cur_cluster_size = 0;
    pair<id_t, id_t> cur_cluster = make_pair(-1, -1);
    pair<id_t, id_t> default_cluster = make_pair(-1, -1); //for checking if cluster has been reinitialized as empty again.
    auto cur_snarl_i = snarl_roots.begin();
    while (cur_snarl_i != snarl_roots.end())
    {
        SubHandleGraph snarl_graph = extract_subgraph(_graph, cur_snarl_i->first, cur_snarl_i->second);

        // //todo: begin debug_code:
        // cerr << "cur_cluster: " << cur_cluster.first << " " << cur_cluster.second << endl;
        // if (cur_cluster.first != -1)
        // {
        //     cerr << " reverse? " << _graph.get_is_reverse(_graph.get_handle(cur_cluster.first)) << " " << _graph.get_is_reverse(_graph.get_handle(cur_cluster.second)) << endl;
        // }
        // cerr << "cur_snarl: " << cur_snarl_i->first << " " << cur_snarl_i->second << endl;
        // if (cur_snarl_i->first != -1)
        // {
        //     cerr << " reverse? " << _graph.get_is_reverse(_graph.get_handle(cur_snarl_i->first)) << " " << _graph.get_is_reverse(_graph.get_handle(cur_snarl_i->second)) << endl;
        // }
        // for (auto node_id : {cur_cluster.first, cur_cluster.second, cur_snarl_i->first, cur_snarl_i->second})
        // {
        //     cerr << "follow edges of " << node_id << endl;
        //     if (node_id != -1)
        //     {
        //         cerr << "go left is true: " << endl;
        //         _graph.follow_edges(_graph.get_handle(node_id), true, [&](handle_t handle)
        //         {
        //             cerr << _graph.get_id(handle) << " " << _graph.get_sequence(handle) << endl;
        //         });
        //         cerr << "go left is false: " << endl;
        //         _graph.follow_edges(_graph.get_handle(node_id), false, [&](handle_t handle)
        //         {
        //             cerr << _graph.get_id(handle) << " " << _graph.get_sequence(handle) << endl;
        //         });
        //     }

            
        // }

        // if (cur_snarl_i->first == 177070) 
        // {
        //     cerr << "for each handle in cur snarl: " << endl;
        //     cerr << _graph.get_sequence(_graph.get_handle(177073)) << endl;
        //     cerr << _graph.get_sequence(_graph.get_handle(177068)) << endl;
        //     snarl_graph.for_each_handle([&](handle_t handle)
        //     {
        //         cerr << snarl_graph.get_id(handle) << " " << snarl_graph.get_sequence(handle) << endl;
        //     });
        //     cerr << "follow edges of 177072:" << endl;
        //     _graph.follow_edges(_graph.get_handle(177072), true, [&](handle_t handle)
        //     {
        //         cerr << _graph.get_id(handle) << " " << _graph.get_sequence(handle) << endl;
        //     });
        //     cerr << "go left is false: " << endl;
        //     _graph.follow_edges(_graph.get_handle(177072), false, [&](handle_t handle)
        //     {
        //         cerr << _graph.get_id(handle) << " " << _graph.get_sequence(handle) << endl;
        //     });
        // }
        //todo: end debug_code:

        //case: initialize cur_cluster:
        int cur_snarl_size = distance_index.maximum_distance(cur_snarl_i->first, false, 0, cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1);


        if (cur_cluster == default_cluster)
        {
            //then initialize the cur_cluster as the current snarl, or
            //if the current snarl fails the tests, move on to the next snarl.  
            int snarl_gap = 0; //because there is no cur_cluster to get distance to snarl
            if (!test_snarl_for_clustering(snarl_graph, cur_snarl_size))
            {
                // cerr << "initial snarl fails test" << endl;
                skipped_snarls.push_back(*cur_snarl_i);
                cur_snarl_i++;
                continue;
            }
            else
            {
                // cerr << "initial snarl continues." << endl;

                // Initialize the cur_cluster as the current snarl.
                cur_cluster = *cur_snarl_i;
                cur_cluster_size = cur_snarl_size;
                cur_snarl_i++;
                continue;
            }
        }

        // we are considering extending the cluster. So, let's get the distance between the end of the cur_cluster and the beginning of cur_snarl.
        //the smaller of these two gap calculations must not be the one containing a snarl inside it. Thus it is the true gap calculation. 
        // cerr << "distance of graph: " << distance_index.maximum_distance(1, false, 0, 18, false, _graph.get_sequence(_graph.get_handle(18)).size() - 1) << endl;
        // cerr << "distance of 15 to 18: " << distance_index.maximum_distance(15, false, 0, 18, false, _graph.get_sequence(_graph.get_handle(18)).size() - 1) << endl;
        // cerr << "distance of 12 to 18: " << distance_index.maximum_distance(12, false, 0, 18, false, _graph.get_sequence(_graph.get_handle(18)).size() - 1) << endl;
        // cerr << "distance of 12 to 12: " << distance_index.maximum_distance(12, false, 0, 12, false, _graph.get_sequence(_graph.get_handle(12)).size() - 1) << endl;
        // cerr << "reverse distance of 12 to 12: " << distance_index.maximum_distance(12, false, _graph.get_sequence(_graph.get_handle(12)).size() - 1, 12, false, 0) << endl;
        
        // cerr << "cur_cluster.second " << cur_cluster.second << endl;
        // cerr << "cur_snarl_i->first " << cur_snarl_i->first << endl;
        // cerr << "cur_snarl_i->second " << cur_snarl_i->second << endl;
        // cerr << "cur_cluster.first " << cur_cluster.first << endl;
        int right_gap = INT_MAX;
        int left_gap = INT_MAX;
        int snarl_gap = INT_MAX;
        // if (cur_cluster.second == cur_snarl_i->first)
        // {
        //     //then the right_gap is correct, and the gap is the length of the shared handle.
        //     right_gap = _graph.get_sequence(_graph.get_handle(cur_cluster.second)).size();
        //     snarl_gap = right_gap;

        // }
        // else if (cur_snarl_i->second == cur_cluster.first)
        // {
        //     //then the left_gap is correct, and the gap is the length of the shared handle.
        //     left_gap = _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size();
        //     snarl_gap = left_gap;
        // }
        // else
        // {
        //     //we have to calculate minimum distance:
        //     right_gap = distance_index.maximum_distance(cur_cluster.second, false, _graph.get_sequence(_graph.get_handle(cur_cluster.second)).size() - 1, cur_snarl_i->first, false, 0);
        //     left_gap = distance_index.maximum_distance(cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1, cur_cluster.first, false, 0);
        //     if (right_gap == -1)
        //     {
        //         snarl_gap = left_gap;
        //     }
        //     else if (left_gap == -1)
        //     {
        //         snarl_gap = right_gap;
        //     }
        //     else
        //     {
        //         snarl_gap = min(left_gap, right_gap);
        //         // cerr << "Although the snarls do not overlap, both distances are this should never happen, based on what I understand."
        //     }
        // }
        // cerr << "distance_index.maximum_distance(12, false, _graph.get_sequence(_graph.get_handle(12)).size() - 1, 12, false, 0)" << distance_index.maximum_distance(12, false, _graph.get_sequence(_graph.get_handle(12)).size() - 1, 12, false, 0) << endl;
        // cerr << "distance_index.maximum_distance(15, false, _graph.get_sequence(_graph.get_handle(15)).size() - 1, 15, false, 0)" << distance_index.maximum_distance(15, false, _graph.get_sequence(_graph.get_handle(15)).size() - 1, 15, false, 0) << endl;
        // cerr << "distance_index.maximum_distance(cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1, cur_cluster.first, false, 0) " << distance_index.maximum_distance(cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1, cur_cluster.first, false, 0) << endl;

        // int test_1 = distance_index.maximum_distance(12, false, _graph.get_sequence(_graph.get_handle(12)).size() - 1, 12, false, 0);
        // int test_2 = distance_index.maximum_distance(15, false, _graph.get_sequence(_graph.get_handle(15)).size() - 1, 15, false, 0);
        // int test_3 = distance_index.maximum_distance(cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1, cur_cluster.first, false, 0);
        // cerr << test_1 << " " << test_2 << " " << test_3 << endl;
        



        right_gap = distance_index.maximum_distance(cur_cluster.second, false, _graph.get_sequence(_graph.get_handle(cur_cluster.second)).size() - 1, cur_snarl_i->first, false, 0, true) + 1; //+1 to inclusively count the last base in the handle; .size() without -1 in maximum_dist has undefined behavior.
        left_gap = distance_index.maximum_distance(cur_snarl_i->second, false, _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() - 1, cur_cluster.first, false, 0, true) + 1; //+1 to inclusively count the last base in the handle; .size() without -1 in maximum_dist has undefined behavior.
        snarl_gap = min(left_gap, right_gap);
        // if (cur_snarl_i->first == 177070) 
        // {
        //     cerr << "cur_cluster.second " << cur_cluster.second << " _graph.get_sequence(_graph.get_handle(cur_cluster.second)).size() " << _graph.get_sequence(_graph.get_handle(cur_cluster.second)).size() << " cur_snarl_i->first " << cur_snarl_i->first << endl;
        //     cerr << "cur_snarl_i->second " << cur_snarl_i->second << " _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() " << _graph.get_sequence(_graph.get_handle(cur_snarl_i->second)).size() << " cur_cluster.first " << cur_cluster.first << endl;

        //     cerr << " left_gap " << left_gap << " right_gap " << right_gap << " snarl_gap " << snarl_gap << endl; //case: cur_snarl fails the tests, so reset the cluster.
        // }
        if (!test_snarl_for_clustering(snarl_graph, cur_snarl_size))
        {
            // 
            // cerr << "cur_snarl fails the tests, reset cluster." << endl;
            // discard the cur_snarl, since it's too large.
            skipped_snarls.push_back(*cur_snarl_i);
            // save the current cluster and make a new one.
            clustered_snarls.push_back(cur_cluster);
            cur_cluster_size = 0;
            cur_cluster = make_pair(-1, -1);
            cur_snarl_i++;

            continue;
        }
        else if (right_gap == 0 && left_gap == 0)
        {
            // cerr << "both gaps are 0. Starting a new cluster." << endl;
            //this means that the new snarl is on a separate connected component of the graph (e.g. another chromosome).
            //save the previous cluster and make this the new cluster.
            clustered_snarls.push_back(cur_cluster);
            cur_cluster_size = cur_snarl_size; //+1 to inclusively count the last base in the handle; .size() without -1 in maximum_dist has undefined behavior.
            cur_cluster = *cur_snarl_i;
            cur_snarl_i++;
            continue;

        }
        else if (snarl_gap > _max_region_gap)
        {
            //reset the cluster, but do not skip/discard the cur_snarl.
            clustered_snarls.push_back(cur_cluster);
            cur_cluster_size = 0;
            cur_cluster = make_pair(-1, -1);
            continue;
        }


        //case: the snarl can be added to the cluster. 
        int extended_cluster_size = cur_cluster_size + snarl_gap + cur_snarl_size;
        
        // if ((cur_cluster_size + max_size + cur_snarl_size) <= _max_region_size)
        if (extended_cluster_size <= _max_region_size)
        {
            //combine the cluster with the new cluster.
            if (left_gap < right_gap)
            {
                cur_cluster = make_pair(cur_snarl_i->first, cur_cluster.second);
            }
            else if (right_gap < left_gap)
            {
                cur_cluster = make_pair(cur_cluster.first, cur_snarl_i->second);
            }
            else
            {
                cerr << "in snarl " << cur_snarl_i->first << " " << cur_snarl_i->second << endl;
                cerr << "unfortunately, distance index's maximum_distance acts differently than I thought. two sizes that should have been different were treated as identical." << endl;
                cerr << "right_gap between cur_cluster and the cur_snarl is: " << right_gap << endl;
                cerr << "left_gap between cur_cluster and the cur_snarl is: " << left_gap << endl;
                cerr << "the cur_cluster in question is: " << cur_cluster.first << " " << cur_cluster.second << endl;
                //todo: debug-code
                //get path names on the nodes:
                cerr << "these are the path names running through the cur_cluster: " << endl;
                vector<step_handle_t> steps_of_cluster =  _graph.steps_of_handle(_graph.get_handle(cur_cluster.second), false);
                for (auto step : steps_of_cluster)
                {
                    cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " " << endl;

                }
                cerr << endl;
                cerr << "these are the path names running through the cur_snarl: " << endl;
                vector<step_handle_t> steps_of_snarl =  _graph.steps_of_handle(_graph.get_handle(cur_snarl_i->first), false);
                for (auto step : steps_of_snarl)
                {
                    cerr << _graph.get_path_name(_graph.get_path_handle_of_step(step)) << " " << endl;

                }
                cerr << endl;
                //print the gap sizes:
                

                //todo: end debug-code
                exit(1);
            }
            
            cur_cluster_size = extended_cluster_size;
            cur_snarl_i++;
        }
        //case: the snarl cannot be added to the cluster because it would make the cluster
        // too large. Reset the cluster with cur_snarl.
        else
        {
            // save the current cluster and initialize the new one with the new region.
            clustered_snarls.push_back(cur_cluster);
            cur_cluster_size = cur_snarl_size;
            cur_cluster = *cur_snarl_i;
            cur_snarl_i++;
            _reset_cluster_because_too_big += 1;
            continue;
        }
        
    }
    if (cur_cluster != default_cluster)
    {
        //then there was a cluster being extended when the while loop ended.
        clustered_snarls.push_back(cur_cluster);
    }
    return clustered_snarls;
}

/// @brief Determines if a snarl is fit for adding to the current cluster.
/// @param snarl_graph the SubHandleGraph containing the snarl.
/// @param snarl_size the size of the cur_snarl
/// @param snarl_gap the gap between the cur_snarl and the cluster that is being extended. (set to 0 if this is the first snarl in the cluster.) 
/// @return bool, true if passes tests; false if fails tests.
bool NormalizeRegionFinder::test_snarl_for_clustering(const HandleGraph& snarl_graph, const int snarl_size)
{
    if(snarl_size > _max_region_size)
    {
        // cerr << "snarl size (" << snarl_size << ") > " << _max_region_size << endl;
        return false;
    }
    if (!handlealgs::is_acyclic(&snarl_graph))
    {
        // cerr << "snarl is cyclic." << endl;
        return false;
    }
    return true;
}

}
}



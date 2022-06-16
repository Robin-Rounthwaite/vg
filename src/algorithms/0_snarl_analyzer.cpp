#include "0_snarl_analyzer.hpp"
#include "0_oo_normalize_snarls.hpp"
#include "../snarls.hpp"

namespace vg {
namespace algorithms{

void print_handles_in_snarl(const HandleGraph& graph, const id_t& source, const id_t& sink, const int& max_search_dist, int autostop/*=20*/)
{
    // vector<int> mapping_nodes { 803806, 803807, 803809, 803810, 803812, 803813, 803815, 803816, 803817, 803818, 803821, 803822, 803823, 803824, 803825, 803826, 803827, 803828, 803829, 803830, 803831, 803832, 803833, 803834 };
    // for (id_t node : mapping_nodes)
    // {
    //     cerr << "printing info for node " << node << ":" << endl;
    //     cerr << "node next:" << endl;
    //     graph.follow_edges(graph.get_handle(node), false, [&](const handle_t &next) {
    //         cerr << "     " << graph.get_id(next) << " ";
    //     });
    //     cerr << endl;
    //     cerr << "node prev:" << endl;
    //     graph.follow_edges(graph.get_handle(node), true, [&](const handle_t &next) {
    //         cerr << "     " << graph.get_id(next) << " ";
    //     });
    //     cerr << endl;
    //     cerr << endl;
    // }
    
    // cerr << "source size: " << graph.get_sequence(graph.get_handle(source)).size() << endl;
    // cerr << "source next:" << endl;
    // graph.follow_edges(graph.get_handle(source), false, [&](const handle_t &next) {
    //     cerr << "     " << graph.get_id(next) << " ";
    // });
    // cerr << endl;
    // cerr << "source prev:" << endl;
    // graph.follow_edges(graph.get_handle(source), true, [&](const handle_t &next) {
    //     cerr << "     " << graph.get_id(next) << " ";
    // });
    // cerr << endl;
    // cerr << "sink next:" << endl;
    // graph.follow_edges(graph.get_handle(sink), false, [&](const handle_t &next) {
    //     cerr << "     " << graph.get_id(next) << " ";
    // });
    // cerr << endl;
    // cerr << "sink prev:" << endl;
    // graph.follow_edges(graph.get_handle(sink), true, [&](const handle_t &next) {
    //     cerr << "     " << graph.get_id(next) << " ";
    // });
    // cerr << endl;
    
    if (!graph.has_node(source))
    {
        cerr << "graph does not contain node " << source << ". exiting." << endl;
        exit(1);
    }
    else if (!graph.has_node(sink))
    {
        cerr << "graph does not contain node " << sink << ". exiting." << endl;
        exit(1);
    }
    //Note: extract_subgraph here assumes that the source is leftmost node of the region.
    //todo: somehow avoid that assumption? Can I look up the snarl in snarl_roots, for example?
    cerr << "searching with leftmost handle as " << source << " and rightmost handle as " << sink << endl;
    // int max_search_dist = 500*32; // 500 standard handles.
    SubHandleGraph snarl = SnarlNormalizer::extract_subgraph(graph, source, sink);
    // SubHandleGraph snarl = extract_subgraph(graph, source, sink, false, max_search_dist, autostop);
    // if (snarl.get_node_count() == 0){
    //     cerr << "failed to extract sequence; exceeded max size. Trying opposite orientation." << endl;
    //     snarl = extract_subgraph(graph, source, sink, true, max_search_dist, autostop);
    //     if (snarl.get_node_count() == 0)
    //     {
    //         cerr << "failed to extract sequence; exceeded max size. Neither orientation had the far node within " << max_search_dist << " bases. However, both nodes exist in the graph." << endl;
    //         exit(1);
    //     }
    // } 

    
    // cout << "node ids of snarl with source " << source << " and sink " << sink << endl; 
    int cur_search_dist = 0;
    snarl.for_each_handle([&](const handle_t handle) 
    {
        cur_search_dist += graph.get_sequence(handle).size();
        cout << graph.get_id(handle) << " ";
    });
}


SnarlAnalyzer::SnarlAnalyzer(const HandleGraph& graph, ifstream &snarl_stream, bool skip_source_sink)
    :_graph(graph) 
{
    // extract the stats on each snarl in the _graph.
    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);

    vector<const Snarl *> snarl_roots = snarl_manager->top_level_snarls();
    for (auto roots : snarl_roots) {
        SubHandleGraph snarl = extract_subgraph(graph, roots->start().node_id(), roots->end().node_id(), roots->start().backward(), INT_MAX, -1); //note: autostop disabled becauses rightmost handle is guaranteed to be a sink.

        int snarl_size = 0;
        if (skip_source_sink) // skip the source and sink to avoid double-counting of seq.
        {
            snarl.for_each_handle([&](const handle_t handle) {
                if (_graph.get_id(handle) == roots->start().node_id() || _graph.get_id(handle) == roots->end().node_id())
                {
                    // count the number of bases in the snarl, except if it's the source or sink.
                    // this is an imperfect solution to avoid the double-counting problem for 
                    // snarls that share a handle at their borders.
                    snarl_size += snarl.get_sequence(handle).size();
                }
            });

        }
        else // allow potential double-counting of sequence in source and sink.
        {
            snarl.for_each_handle([&](const handle_t handle) {
                // count the number of bases in the snarl.
                snarl_size += snarl.get_sequence(handle).size();
            });

        }

        // save the stats on the snarl for later use.
        _snarl_sources.push_back(roots->start().node_id());
        _snarl_sinks.push_back(roots->end().node_id());
        _snarl_sizes.push_back(snarl_size);
    }

    // TODO: if desired, include this count of all snarls, including non-top-level.
    // int general_count = 0;
    // snarl_manager->for_each_snarl_preorder([&](const vg::Snarl * ignored){
    //     general_count++;
    // });
    // cerr << "number of total snarls in graph: " << general_count << endl;

}

void SnarlAnalyzer::output_snarl_sizes(string& file_name)
{
    std::ofstream outfile;
    outfile.open(file_name);

    for (int i=0; i != _snarl_sources.size(); i++)
    {
        outfile << _snarl_sources[i] << "\t" << _snarl_sinks[i] << "\t" << _snarl_sizes[i] << endl;
    }
}



// // TODO: Undo this terrible copy-paste from SnarlNormalizer method, and somehow use 
// // TODO:   SnarlNormalizer's code instead. (make this a non-object available function?) 
// // Given a start and end node id, construct an extract subgraph between the two nodes
// // (inclusive). Arguments:
// //      _graph: a pathhandlegraph containing the snarl with embedded paths.
// //      source_id: the source of the snarl of interest.
// //      sink_id: the sink of the snarl of interest.
// //      autostop: stop extraction of subgraph after finding the rightmost handle in the graph, and then adding up to autostop handles to the graph. (if the rightmost handle is not a sink for the graph). Set to -1 to disable.
// // Returns:
// //      a SubHandleGraph containing only the handles in _graph that are between start_id
// //      and sink_id.
// SubHandleGraph extract_subgraph(const HandleGraph &graph,
//                                                  id_t source_id,
//                                                  id_t sink_id,
//                                                  const bool backwards,
//                                                  int max_search_dist,
//                                                  int autostop) 
// {
//     // cerr << "extract_subgraph has source and sink: " << source_id << " " << sink_id << endl; 
//     // because algorithm moves left to right, determine leftmost and rightmost nodes.
//     id_t leftmost_id;
//     id_t rightmost_id;
//     // if snarl's "backwards," source is rightmost node, sink is leftmost.
//     if (backwards) 
//     {
//         leftmost_id = sink_id;
//         rightmost_id = source_id;
//     }
//     else 
//     {
//         leftmost_id = source_id;
//         rightmost_id = sink_id;
//     }
//     // cerr << "extract_subgraph" << endl;
//     /// make a subgraph containing only nodes of interest. (e.g. a snarl)
//     // make empty subgraph
//     SubHandleGraph subgraph = SubHandleGraph(&graph);

//     int cur_search_dist = 0;

//     unordered_set<id_t> visited;  // to avoid counting the same node twice.
//     unordered_set<id_t> to_visit; // nodes found that belong in the subgraph.

//     // initialize with leftmost_handle (because we move only to the right of leftmost_handle):
//     handle_t leftmost_handle = graph.get_handle(leftmost_id);
//     subgraph.add_handle(leftmost_handle);
//     visited.insert(graph.get_id(leftmost_handle));
//     cur_search_dist += graph.get_sequence(leftmost_handle).size();

//     // look only to the right of leftmost_handle
//     graph.follow_edges(leftmost_handle, false, [&](const handle_t &handle) {
//         // mark the nodes to come as to_visit
//         if (visited.find(graph.get_id(handle)) == visited.end()) {
//             to_visit.insert(graph.get_id(handle));
//         }
//     });

//     int autostop_count = -1;
//     /// explore the rest of the snarl:
//     while (to_visit.size() != 0 && autostop_count!=autostop) {
//         // remove cur_handle from to_visit
//         unordered_set<id_t>::iterator cur_index = to_visit.begin();
//         handle_t cur_handle = graph.get_handle(*cur_index);
//         if (graph.get_id(cur_handle) == rightmost_id && autostop_count == -1 && autostop != -1) // when autostop == -1, it's disabled.
//         {
//             autostop_count = 1;
//         }
//         if (autostop_count != -1 && autostop != -1)
//         {
//             if (autostop_count == autostop)
//             {
//                 break; // we have found the rightmost node, and afterwards extended the graph by autostop handles. End graph extension.
//             }
//             autostop_count++;
//         }

//         to_visit.erase(cur_index);

//         /// visit cur_handle
//         visited.insert(graph.get_id(cur_handle));

//         subgraph.add_handle(cur_handle);
//         if (cur_search_dist >= max_search_dist)
//         {
//             cerr << "exceeded max_search_dist. Returning empty SubHandleGraph" << endl;
//             return SubHandleGraph(&graph);
//         }
//         cur_search_dist += graph.get_sequence(cur_handle).size();


//         if (graph.get_id(cur_handle) != rightmost_id) { // don't iterate past rightmost node! (note: this code doesn't work if there are multiple paths past the rightmost node region. Hence the autostop feature.)
//             // look for all nodes connected to cur_handle that need to be added
//             // looking to the left,
//             graph.follow_edges(cur_handle, true, [&](const handle_t &handle) {
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//             // looking to the right,
//             graph.follow_edges(cur_handle, false, [&](const handle_t &handle) {
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//         }
//     }
//     return subgraph;
// }

}//algorithms
}//vg
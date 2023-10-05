// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce
// more efficient representations of snarls.

//todo: remove unnecessary includes.
#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include <bdsg/snarl_distance_index.hpp>
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../snarl_distance_index.hpp"

#include "subcommand.hpp"

// todo: should be able to remove '../../include/...' and replace with e.g.
// <bdsg/hash...>
#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"

#include "../io/save_handle_graph.hpp"

#include "../algorithms/0_get_parallel_normalize_regions.hpp"
#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../algorithms/0_snarl_analyzer.hpp"
#include "../algorithms/0_update_gbwt_wrapper.hpp"
#include "../algorithms/0_snarl_sequence_finder.hpp"


#include "../snarls.hpp"
#include "../clip.hpp"
#include <bdsg/overlays/path_position_overlays.hpp>

#include <chrono> // for high_resolution_clock

using namespace std;
using namespace vg;
using namespace vg::subcommand;



SubHandleGraph extract_subgraph(const HandleGraph &graph,
                                                 const vg::id_t leftmost_id,
                                                 const vg::id_t rightmost_id) {
    /// make a subgraph containing only nodes of interest. (e.g. a snarl)
    // make empty subgraph
    SubHandleGraph subgraph = SubHandleGraph(&graph);

    // this bool tracks when we find the rightmost_id while extending the snarl. If we 
    // never find it, then we must have been passed leftmost_id and rightmost_id in 
    // reverse order. In that case, we'll throw an error.
    bool found_rightmost = false;

    unordered_set<vg::id_t> visited;  // to avoid counting the same node twice.
    unordered_set<vg::id_t> to_visit; // nodes found that belong in the subgraph.

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
        unordered_set<vg::id_t>::iterator cur_index = to_visit.begin();
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


void help_normalize(char **argv) {
  cerr
        << "usage: " << argv[0] << " normalize [options] <graph.vg> >[normalized.vg]"
        << endl
        << "Modifies snarls, outputs modified on stdout." << endl
        << endl
        << "options:" << endl
        << "necessary options: " << endl
        << "    -g, --gbwt       gbwt corresponding to input graph." << endl
        << "    -d, --distance_index       distance index corresponding to input graph." << endl
        << "recommended options: " << endl
        << "    -r, --gbwt_graph       the gbwt_graph corresponding to the input graph and gbwt. If not supplied, will be temporarily generated." << endl
        << "    -o, --output_gbwt   name for the gbwt corresponding to the normalized graph. Default is normalized.gbwt."
        << "optional options: " << endl
        << "    -l, --max_handle_size       currently, default is 32, to match the default "
            "handle size of most graphs. This "
            "changes what the largest size of handle can be in normalized regions. "
        << endl
        //todo: allow the user to specify normalize and update_gbwt threads separately.
        << "    -t, --threads      The number of threads used in the normalization process. In addition to the normalization step, it affects the update_gbwt step, which can grow very long without multithreading. To prevent overusing memory, the update_gbwt step will only use up to 14 threads (even if more is specified), but the normalization step can use as many as specified, so long as there are regions that still need realigning.  Default:14."
        << endl
        << "    -h, --help      print this help info." << endl;
}

int main_normalize(int argc, char **argv) {

    if (argc == 2) {
        help_normalize(argv);
        return 1;
    }

    string gbwt_file;
    string distance_index_file;
    string gbwt_graph_file;
    string output_gbwt_file = "normalized.gbwt";
    int max_handle_size = 32;
    int threads = 14;

    //todo: make options for these variables.
    // dictates the number of snarls that may be clustered into a single region.
        //todo: make it be based on sequence length-related metric instead.
    int max_snarls_per_region = 1;
    int max_snarl_spacing = 1;
    int max_strings_per_alignment = INT_MAX; // default cutoff used to be 200 threads in a snarl.

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {{"help", no_argument, 0, 'h'},
            {"gbwt", required_argument, 0, 'g'},
            {"distance_index", required_argument, 0, 'd'},
            {"gbwt_graph", required_argument, 0, 'r'},
            {"output_gbwt", required_argument, 0, 'o'},
            {"max_handle_size", required_argument, 0, 'l'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "hg:d:r:o:l:t:", long_options,
                        &option_index);
        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {

        case 'g':
            gbwt_file = optarg;
            break;

        case 'd':
            distance_index_file = optarg;
            break;

        case 'r':
            gbwt_graph_file = optarg;
            break;

        case 'o':
            output_gbwt_file = optarg;
            break;

        case 'l':
            max_handle_size = parse<int>(optarg);
            break;

        case 't':
            threads = parse<int>(optarg);
            break;

        case 'h':
        case '?':
            help_normalize(argv);
            exit(1);
            break;

        default:
            cerr << "error:[vg normalize] abort" << endl;
            abort();
        }
    }


    // =======loading files for normalization:=======
    cerr << "=======loading files for normalization=======" << endl;
    auto start = chrono::high_resolution_clock::now();

    //graph
    cerr << "loading graph" << endl;
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    // gbwt
    cerr << "loading gbwt" << endl;
    ifstream gbwt_stream;
    gbwt_stream.open(gbwt_file);    
    unique_ptr<gbwt::GBWT> gbwt;
    gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);

    // gbwt graph 
    cerr << "loading gbwt graph" << endl;
    unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph;
    if (gbwt_graph_file.size() == 0)
    {
      string gbwt_graph_output_file = gbwt_file + ".gg";
      cerr << "gbwt_graph option is empty. Making new GBWTGraph from gbwt and graph. Saving as " << gbwt_graph_output_file << endl;
      gbwtgraph::GBWTGraph new_gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
      save_gbwtgraph(new_gbwt_graph, gbwt_graph_output_file);
      //todo: find way to load gbwtgraph's unique pointer without saving and then reloading file.
      gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_output_file);
      gbwt_graph->set_gbwt(*gbwt);
    }
    else 
    {
      gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_file);
      gbwt_graph->set_gbwt(*gbwt);
    }


    // //todo: begin debug code
    // SubHandleGraph subgraph = extract_subgraph(*graph, 1, 15);
    // vg::algorithms::SnarlSequenceFinder seq_finder = vg::algorithms::SnarlSequenceFinder(*graph, subgraph, *gbwt_graph, 1, 15, false);
    // std::tuple<std::vector<std::vector<gbwtgraph::handle_t>>, 
    // std::vector<std::vector<gbwtgraph::handle_t>>, 
    // std::unordered_set<handlegraph::nid_t>>
    //     seqs = seq_finder.find_gbwt_haps();
        
    // unordered_set<string> haplotype_strings;
    // cerr << "ids of haplotypes:" << endl;
    // for (vector<handle_t> haplotype_handles : get<0>(seqs)) {
    //     string hap;
    //     for (handle_t handle : haplotype_handles) {
    //         cerr << " " << gbwt_graph->get_id(handle);

    //         // hap += _gbwt_graph.get_sequence(handle);
    //         hap += gbwt_graph->get_sequence(handle);
    //     }
    //     cerr << endl;
    //     haplotype_strings.emplace(hap);
    //     // cerr << "hap: " << hap << endl;
    // }

    // // gbwt::SearchState cur_state = gbwt_graph->get_state(gbwt_graph->get_handle(1));
    // // gbwt::SearchState last_state = gbwt_graph->get_state(gbwt_graph->get_handle(15));
    // // vector<string> paths;
    // // paths.
    // // while (cur_state != last_state)
    // // {
    // //     cur_state.
    // // }
    // // gbwt_graph->get_sequence()
    
    // cerr << "gbwt graph contents: " << endl;
    // gbwt_graph->for_each_handle([&](handle_t handle){

    //     cerr << gbwt_graph->get_id(handle) << " " << gbwt_graph->get_sequence(handle) << " | left:";
    //     gbwt_graph->follow_edges(handle, true, [&](handle_t left)
    //     {
    //         cerr << " " << gbwt_graph->get_id(left);
    //     });
    //     cerr << " | right:";
    //     gbwt_graph->follow_edges(handle, false, [&](handle_t right)
    //     {
    //         cerr << " " << gbwt_graph->get_id(right);
    //     });
    //     cerr << endl;
    //     cerr << "does the gbwt contiain this node? " << gbwt->contains(gbwt::Node::encode(gbwt_graph->get_id(handle), false)) << endl;
    //     // cerr << "does the gbwt contiain this node? " << _gbwt.contains(gbwt::Node::encode(snarl.get_id(handle), false)) << endl;
    //     cerr << "does the gbwt_graph contian this node? " << gbwt_graph->has_node(gbwt_graph->get_id(handle)) << endl;

        
    //     // gbwt->contains(node);
    // });
    
    // //todo: end debug code

    // v2 distance index
    cerr << "loading distance index" << endl;
    auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_index_file);

    // // =======isolating normalize regions for multithreading=======
    cerr << "=======isolating normalize regions for multithreading=======" << endl;
    cerr << "getting non-isolated normalize regions." << endl;
    vector<pair<vg::id_t, vg::id_t>> snarl_roots;
    distance_index->traverse_decomposition([&] (const bdsg::net_handle_t& snarl) {
        handle_t inward_source_handle = distance_index->get_handle(distance_index->get_bound(snarl, false, true), &*graph);
        handle_t outward_sink_handle = distance_index->get_handle(distance_index->get_bound(snarl, true, false), &*graph);

        //source is "rightmost" handle if inward_source is_reverse=true. (where the
        //direction "right" is forward for the default orientation of a handle generated
        //on the node)
        if (graph->get_is_reverse(inward_source_handle) && graph->get_is_reverse(outward_sink_handle))
        {
            snarl_roots.push_back(
                make_pair(
                    graph->get_id(outward_sink_handle),
                    graph->get_id(inward_source_handle)
                    ));
        }
        else if (!graph->get_is_reverse(inward_source_handle) && !graph->get_is_reverse(outward_sink_handle))
        {
          snarl_roots.push_back(
            make_pair(
                graph->get_id(inward_source_handle),
                graph->get_id(outward_sink_handle)
                ));
        }
        else
        {
            //something unexpected is happening. There shouldn't ever be a reverse source but non-reversed sink, or vice-versa. Throw error. 
            //todo: permit this by turning snarl_roots into pairs of handles instead of id_ts. 
            cerr << "error:[vg normalize] there is a snarl with a source and sink that have mismatched orientations. This is currently unsupported. " << endl;
            exit(1);
        }
        return true;
    },
    [&](const bdsg::net_handle_t& chain) {return true;},
    [&] (const bdsg::net_handle_t& node) {return true;});


    // cerr << "calling weakly connected components before region finder" << endl;
    // vector<unordered_set<nid_t>> components_0 = handlegraph::algorithms::weakly_connected_components(&(*gbwt_graph));
    // cerr << "components_0.size() " << components_0.size() << endl;
    //todo: re-uncomment begin .

    // cerr << "segregating regions for parallelization" << endl; //so that edits in two parallel jobs can't touch the same node.
    // vg::algorithms::NormalizeRegionFinder region_finder = vg::algorithms::NormalizeRegionFinder(*graph, max_snarls_per_region, max_snarl_spacing);

    // // cerr << "calling weakly connected components after region finder" << endl;
    // // vector<unordered_set<nid_t>> components_0_5 = handlegraph::algorithms::weakly_connected_components(&(*gbwt_graph));

    // std::vector<std::pair<vg::id_t, vg::id_t>> parallel_normalize_regions;
    // vector< pair< pair< vg::id_t, vg::id_t >, vg::id_t > > desegregation_candidates;
    // std::vector<vg::RebuildJob::mapping_type> parallel_regions_gbwt_updates = region_finder.get_parallel_normalize_regions(snarl_roots, parallel_normalize_regions, desegregation_candidates);
    // // cerr << "first parallel_regions_gbwt_update: " << parallel_regions_gbwt_updates.front().first << " " << parallel_regions_gbwt_updates.front().second << endl; 
    //todo: re-uncomment end.
    handle_t debug_handle = graph->get_handle(168094);
    cerr << "gbwt->contains(graph->get_id(debug_handle)) " << gbwt->contains(graph->get_id(debug_handle)) << endl;
    cerr << "debug handle contents: " << graph->get_sequence(debug_handle) << endl;
    pair<handle_t, handle_t> new_debugs = graph->divide_handle(debug_handle, graph->get_sequence(debug_handle).size()/2);
    cerr << "new debugs seq: " << graph->get_sequence(new_debugs.first) << " " << graph->get_sequence(new_debugs.second) << endl;
    gbwt::vector_type before;
    before.push_back(gbwt::Node::encode(graph->get_id(debug_handle), false));
    gbwt::vector_type after;
    after.push_back(gbwt::Node::encode(graph->get_id(new_debugs.first), false));
    after.push_back(gbwt::Node::encode(graph->get_id(new_debugs.second), false));
    std::vector<vg::RebuildJob::mapping_type> debug_updates;
    debug_updates.push_back(make_pair(before, after));

    cerr << "updating the gbwt and gbwt_graph with the isolated regions" << endl;
    auto isolated_regions_gbwt_update_start = chrono::high_resolution_clock::now();    

    int gbwt_threads = threads;
    if (threads > 12)
    {
      gbwt_threads = 12;
    }

    // //todo: begin debug code
    // vg::RebuildJob::mapping_type first = parallel_regions_gbwt_updates.front();
    // cerr << "before update: ";
    // for (auto before : first.first)
    // {
    //     cerr << " " << gbwt::Node::id(before);
    // }
    // cerr << endl;
    // cerr << "after update: ";
    // for (auto after : first.second)
    // {
    //     cerr << " " << gbwt::Node::id(after);
    // }
    // cerr << endl;

    // // cerr << " " << first.second << endl;
    // parallel_regions_gbwt_updates.clear();
    // parallel_regions_gbwt_updates.push_back(first);
    // //todo: end debug code
    gbwt::GBWT parallel_regions_gbwt = vg::algorithms::apply_gbwt_changelog(*gbwt_graph, debug_updates, *gbwt, gbwt_threads, false);
    // gbwt::GBWT parallel_regions_gbwt = vg::algorithms::apply_gbwt_changelog(*gbwt_graph, parallel_regions_gbwt_updates, *gbwt, gbwt_threads, false);
    //todo: *******************************************************************************DECOMMENT ALL LINES BELOW THIS FOR NORM:*******************************8888
    // auto isolated_regions_gbwt_update_end = chrono::high_resolution_clock::now();    
    // chrono::duration<double> gbwt_changelog_time = isolated_regions_gbwt_update_end - isolated_regions_gbwt_update_start;
    // cerr << "Time spent generating the updated gbwt after isolating regions for parallelization: " << gbwt_changelog_time.count() << " s" << endl;

    // cerr << "generating the updated gbwt graph" << endl;
    // // make a new gbwt_graph for the parallel_regions_gbwt.
    // gbwtgraph::GBWTGraph parallel_regions_gbwt_graph = gbwtgraph::GBWTGraph(parallel_regions_gbwt, *graph);

    // //todo: begin debug code
    // cerr << "after parallel regions update." << endl;
    // cerr << "gbwt graph contents: " << endl;
    // parallel_regions_gbwt_graph.for_each_handle([&](handle_t handle){

    //     cerr << parallel_regions_gbwt_graph.get_id(handle) << " " << parallel_regions_gbwt_graph.get_sequence(handle) << " | left:";
    //     parallel_regions_gbwt_graph.follow_edges(handle, true, [&](handle_t left)
    //     {
    //         cerr << " " << parallel_regions_gbwt_graph.get_id(left);
    //     });
    //     cerr << " | right:";
    //     parallel_regions_gbwt_graph.follow_edges(handle, false, [&](handle_t right)
    //     {
    //         cerr << " " << parallel_regions_gbwt_graph.get_id(right);
    //     });
    //     cerr << endl;
    //     cerr << "does the gbwt contiain this node? " << parallel_regions_gbwt.contains(gbwt::Node::encode(parallel_regions_gbwt_graph.get_id(handle), false)) << endl;
    //     // cerr << "does the gbwt contiain this node? " << _gbwt.contains(gbwt::Node::encode(snarl.get_id(handle), false)) << endl;
    //     cerr << "does the parallel_regions_gbwt_graph contain this node? " << parallel_regions_gbwt_graph.has_node(parallel_regions_gbwt_graph.get_id(handle)) << endl;

        
    //     // gbwt->contains(node);
    // });
    
    // //todo: end debug code


    // // free some memory, since we don't need the original gbwt anymore.
    // gbwt.reset();
    // // //todo: begin debug code

    // SubHandleGraph subgraph_2 = extract_subgraph(*graph, 1, 15);
    // vg::algorithms::SnarlSequenceFinder seq_finder_2 = vg::algorithms::SnarlSequenceFinder(*graph, subgraph_2, parallel_regions_gbwt_graph, 1, 15, false);
    // std::tuple<std::vector<std::vector<gbwtgraph::handle_t>>, 
    // std::vector<std::vector<gbwtgraph::handle_t>>, 
    // std::unordered_set<handlegraph::nid_t>>
    //     seqs_2 = seq_finder_2.find_gbwt_haps();

    // // unordered_set<string> haplotype_strings;
    // cerr << "ids of haplotypes:" << endl;
    // for (vector<handle_t> haplotype_handles : get<0>(seqs_2)) {
    //     string hap;
    //     for (handle_t handle : haplotype_handles) {
    //         cerr << " " << parallel_regions_gbwt_graph.get_id(handle);

    //         // hap += _gbwt_graph.get_sequence(handle);
    //         // hap += gbwt_graph->get_sequence(handle);
    //     }
    //     cerr << endl;
    //     // haplotype_strings.emplace(hap);
    //     // cerr << "hap: " << hap << endl;
    // }

    // // //todo: end debug code

    // // //todo: delete debug code
    // // cerr << " parallel_regions_gbwt_graph.get_sequence(parallel_regions_gbwt_graph.get_handle(18, false)) " << parallel_regions_gbwt_graph.get_sequence(parallel_regions_gbwt_graph.get_handle(18, false)) << endl;

    // // =======running normalize=======
    // cerr << "=======running normalize=======" << endl;
    // string alignment_algorithm="sPOA"; bool disable_gbwt_update=false; bool debug_print=false; //todo: remove these options or else ensure they are implemented correctly.

    // //todo: note: debug code.
    // // std::set<vg::id_t> nodes_to_delete;
    // // vg::algorithms::SnarlNormalizer normalizer = vg::algorithms::SnarlNormalizer(
    // //   *graph, *gbwt, *gbwt_graph, nodes_to_delete, max_handle_size, max_snarls_per_region, max_snarl_spacing, threads, max_strings_per_alignment, "GBWT", alignment_algorithm, disable_gbwt_update, debug_print);
    // // std::vector<vg::RebuildJob::mapping_type> gbwt_normalize_updates = normalizer.parallel_normalization(snarl_roots);
    // vg::algorithms::SnarlNormalizer normalizer = vg::algorithms::SnarlNormalizer(
    //   *graph, parallel_regions_gbwt, parallel_regions_gbwt_graph, max_handle_size, max_snarls_per_region, max_snarl_spacing, threads, max_strings_per_alignment, "GBWT", alignment_algorithm, disable_gbwt_update, debug_print);

    // cerr << "calling weakly connected components before normalization" << endl;
    // vector<unordered_set<nid_t>> components_1 = handlegraph::algorithms::weakly_connected_components(&(parallel_regions_gbwt_graph));
    // cerr << "components_1.size() " << components_1.size() << endl;

    // std::vector<vg::RebuildJob::mapping_type> gbwt_normalize_updates = normalizer.parallel_normalization(parallel_normalize_regions);

    // cerr << "calling weakly connected components after normalization" << endl;
    // vector<unordered_set<nid_t>> components_2 = handlegraph::algorithms::weakly_connected_components(&(parallel_regions_gbwt_graph));
    // cerr << "components_2.size() " << components_2.size() << endl;


    // //todo: begin debug code
    // cerr << "graph contents: " << endl;
    // graph->for_each_handle([&](handle_t handle){

    //     cerr << graph->get_id(handle) << " " << graph->get_sequence(handle) << " | left:";
    //     graph->follow_edges(handle, true, [&](handle_t left)
    //     {
    //         cerr << " " << graph->get_id(left);
    //     });
    //     cerr << " | right:";
    //     graph->follow_edges(handle, false, [&](handle_t right)
    //     {
    //         cerr << " " << graph->get_id(right);
    //     });
    //     cerr << endl;
    // });
    
    // //todo: end debug code


    // cerr << "=======updating gbwt after normalization and before desegregating regions=======" << endl;

    // // //todo: begin debug code
    // // cerr << "first I'll try running update with just the first gbwt update." << endl;
    // // cerr << "size of gbwt_normalize_updates: " << gbwt_normalize_updates.size() << endl;
    // // auto saved = gbwt_normalize_updates.back();
    // // cerr << saved.first << " | " << saved.second << endl;
    // // gbwt_normalize_updates.clear();
    // // gbwt_normalize_updates.push_back(saved);
    // // cerr << saved.first << " | " << saved.second << endl;
    // // // gbwt::Node::id(node);

    // // cerr << "saving updated graph to file" << endl;
    // // //save normalized graph
    // // vg::io::save_handle_graph(graph.get(), std::cout);
    // // exit(1);

    // // //todo: end debug code

    // cerr << "updating gbwt after normalization" << endl;
    // auto _post_norm_gbwt_update_start = chrono::high_resolution_clock::now();
    // gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(parallel_regions_gbwt_graph, gbwt_normalize_updates, parallel_regions_gbwt, gbwt_threads, false);
    // auto _post_norm_gbwt_update_end = chrono::high_resolution_clock::now();    
    // chrono::duration<double> _post_norm_gbwt_update_elapsed = _post_norm_gbwt_update_end - _post_norm_gbwt_update_start;
    // cerr << "elapsed time from updating gbwt after normalization: " << _post_norm_gbwt_update_elapsed.count() << " s" << endl;


    // cerr << "generating the post-normalization gbwt graph" << endl;
    // // make a new gbwt_graph for the parallel_regions_gbwt.
    // gbwtgraph::GBWTGraph normalized_gbwt_graph = gbwtgraph::GBWTGraph(normalized_gbwt, *graph);

    // cerr << "=======desegregating normalization regions after parallelized normalization=======" << endl;

    // cerr << "desegregating the normalize regions." << endl;
    // vg::algorithms::NormalizeRegionFinder post_norm_region_finder = vg::algorithms::NormalizeRegionFinder(*graph, max_snarls_per_region, max_snarl_spacing);
    // //merges nodes, updates entries in gbwt_normalize to match those merged nodes. (note: is there a way to make this O(n), rather than O(n^2)? Maybe a reverse index... seems possibly a distraction. Could be worth just running desegregate nodes after updating the gbwt, and update the gbwt a third time.)
    // std::vector<vg::RebuildJob::mapping_type> desegregated_regions_gbwt_updates = post_norm_region_finder.desegregate_nodes(desegregation_candidates);

    // cerr << "======preparing and saving output=======" << endl;

    // cerr << "saving updated graph to file" << endl;
    // //save normalized graph
    // vg::io::save_handle_graph(graph.get(), std::cout);

    // cerr << "updating gbwt after de-isolation." << endl;
    // auto _desegregated_regions_gbwt_update_start = chrono::high_resolution_clock::now();    
    // gbwt::GBWT normalized_desegregated_gbwt = vg::algorithms::apply_gbwt_changelog(normalized_gbwt_graph, desegregated_regions_gbwt_updates, normalized_gbwt, gbwt_threads, false);

    // auto _desegregated_regions_gbwt_update_end = chrono::high_resolution_clock::now();    
    // chrono::duration<double> _desegregated_regions_gbwt_update_elapsed = _desegregated_regions_gbwt_update_end - _desegregated_regions_gbwt_update_start;
    // cerr << "elapsed time from updating gbwt after de-isolating regions: " << _desegregated_regions_gbwt_update_elapsed.count() << " s" << endl;

    // cerr << "saving updated gbwt" << endl;
    // save_gbwt(normalized_desegregated_gbwt, output_gbwt_file, true);

    // // Record end time
    // auto finish = std::chrono::high_resolution_clock::now();
    // chrono::duration<double> elapsed = finish - start;
    // cerr << "Total elapsed time: " << elapsed.count() << " s" << endl;


    return 0;
}


// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize) ;

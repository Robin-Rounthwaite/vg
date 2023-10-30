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



// SubHandleGraph extract_subgraph(const HandleGraph &graph,
//                                                  const vg::id_t leftmost_id,
//                                                  const vg::id_t rightmost_id) {
//     /// make a subgraph containing only nodes of interest. (e.g. a snarl)
//     // make empty subgraph
//     SubHandleGraph subgraph = SubHandleGraph(&graph);

//     // this bool tracks when we find the rightmost_id while extending the snarl. If we 
//     // never find it, then we must have been passed leftmost_id and rightmost_id in 
//     // reverse order. In that case, we'll throw an error.
//     bool found_rightmost = false;

//     unordered_set<vg::id_t> visited;  // to avoid counting the same node twice.
//     unordered_set<vg::id_t> to_visit; // nodes found that belong in the subgraph.

//     // initialize with leftmost_handle (because we move only to the right of leftmost_handle):
//     handle_t leftmost_handle = graph.get_handle(leftmost_id);
//     subgraph.add_handle(leftmost_handle);
//     visited.insert(graph.get_id(leftmost_handle));

//     // look only to the right of leftmost_handle
//     graph.follow_edges(leftmost_handle, false, [&](const handle_t handle) {
//         // mark the nodes to come as to_visit
//         if (visited.find(graph.get_id(handle)) == visited.end()) {
//             to_visit.insert(graph.get_id(handle));
//         }
//     });

//     /// explore the rest of the snarl:
//     while (to_visit.size() != 0) {
//         // remove cur_handle from to_visit
//         unordered_set<vg::id_t>::iterator cur_index = to_visit.begin();
//         handle_t cur_handle = graph.get_handle(*cur_index);

//         to_visit.erase(cur_index);

//         /// visit cur_handle
//         visited.insert(graph.get_id(cur_handle));

//         subgraph.add_handle(cur_handle);

//         if (graph.get_id(cur_handle) != rightmost_id) { // don't iterate past rightmost node!
//             // look for all nodes connected to cur_handle that need to be added
//             // looking to the left,
//             graph.follow_edges(cur_handle, true, [&](const handle_t handle) {
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//             // looking to the right,
//             graph.follow_edges(cur_handle, false, [&](const handle_t handle) {
//                 // mark the nodes to come as to_visit
//                 if (visited.find(graph.get_id(handle)) == visited.end()) {
//                     to_visit.insert(graph.get_id(handle));
//                 }
//             });
//         }
//         else
//         {
//             // cerr << "found the righmost node id. Here's the cur_handle: " << graph.get_id(cur_handle) << endl;
//             found_rightmost = true;
//         }
//     }
//     if (!found_rightmost)
//     {
//         // we never found rightmost! We probably were fed a rightmost_id that was actually to the left of the leftmost_id. Throw error.
//         cerr << "error:[vg normalize] in function extract_subgraph, was passed snarl with leftmost_id" << leftmost_id;
//         cerr << " and rightmost_id " << rightmost_id;
//         cerr << ". However, rightmost_id was not found to the right of leftmost_id.";
//         cerr << " Were the ids swapped?" << endl;
//         exit(1);
//     }
//     return subgraph;
// }


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
        // << "    -O, --output_gbwt_graph   name for the gbwt_graph corresponding to the normalized graph and gbwt. Default is normalized.gbwt.gg."
        << "    -l, --max_handle_size       currently, default is 32, to match the default "
            "handle size of most graphs. This "
            "changes what the largest size of handle can be in normalized regions. "
        << endl
        << "    -m, --max_region_size       This determines the maximum size of an alignment region, during cluster_snarls. Default is 750, based on sPOA's known high functionality for MSAs of ~500-1000 bases per string. (But further experimentation could be helpful.)" << endl
        << "    -n, --max_region_gap       This determines the maximum length of sequence found between snarls that can still be clustered into a single snarl. Once exceeded, a new cluster is started with the next snarl. Set to 0 if you want each snarl to be evaluated individually. (not recommended, because small snarls usually need surrounding context to assist with effective realignment.)" << endl
        //todo: allow the user to specify normalize and update_gbwt threads separately.
        << "    -t, --threads      The number of threads used in the normalization process. In addition to the normalization step, it affects the update_gbwt step, which can grow very long without multithreading. To prevent overusing memory, the update_gbwt step will only use up to 14 threads (even if more is specified), but the normalization step can use as many as specified, so long as there are regions that still need realigning.  Default:14."
        << endl
        << "    -s, --output_segregate_regions_only_file       specify the output file for the segregation of the graph into segregated "
            "regions, and generates an updated gbwt. This allows for easy tests of "
            "normalize on the segregated graph without having to update the gbwt. "
            "(debugging tool)" << endl
        << "    -S, --input_segregate_regions_only_file       input from a previous run "
            "that had an output of segregate_regions_only. This. along with the segregated "
            "regions gbwt, gg, and graph, contains all the program needs to continue a run "
            "of normalize after the first segregation of regions. (distance index "
            "may be passed to the program untouched.)." << endl
        << "    -u, --run_tests       run tests to make sure that normalize is still functioning properly." << endl
        << "    -b, --debug_print       print some information during normalization for debugging." << endl
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
    int max_region_size = 750; //todo: when abpoa added, make default region size much larger. However much that memory/time can handle.
    int max_region_gap = 32; //todo: when abpoa added, possibly make default region gap larger.
    int threads = 14;
    string output_segregate_regions_only_file;
    string input_segregate_regions_only_file;
    bool run_tests = false;
    bool debug_print = false;

    //todo: make options for these variables.
    // dictates the number of snarls that may be clustered into a single region.
        //todo: make it be based on sequence length-related metric instead.
    int gbwt_threads = 12;
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
            {"max_region_size", required_argument, 0, 'm'},
            {"max_region_gap", required_argument, 0, 'n'},
            {"threads", required_argument, 0, 't'},
            {"output_segregate_regions_only_file", required_argument, 0, 's'},
            {"input_segregate_regions_only_file", required_argument, 0, 'S'},
            {"run_tests", no_argument, 0, 'u'},
            {"debug_print", no_argument, 0, 'b'},
            {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "hg:d:r:o:l:m:n:t:s:S:ub", long_options,
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

        case 'm':
            max_region_size = parse<int>(optarg);
            break;

        case 'n':
            max_region_gap = parse<int>(optarg);
            break;

        case 't':
            threads = parse<int>(optarg);
            break;

        case 's':
            output_segregate_regions_only_file = optarg;
            break;

        case 'S':
            input_segregate_regions_only_file = optarg;
            break;

        case 'u':
            run_tests = true;
            break;

        case 'b':
            debug_print = true;
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


    std::vector<std::pair<vg::id_t, vg::id_t>> parallel_normalize_regions;
    vector< pair< pair< vg::id_t, vg::id_t >, vg::id_t > > desegregation_candidates; // all id_t are from node ids in the graph 
    std::vector<vg::RebuildJob::mapping_type> parallel_regions_gbwt_updates;
    if (input_segregate_regions_only_file.size() == 0)
    {
        // v2 distance index
        cerr << "loading distance index" << endl;
        auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_index_file);

        // // =======isolating normalize regions for multithreading=======
        cerr << "=======isolating normalize regions for multithreading=======" << endl;
        cerr << "getting non-isolated normalize regions." << endl;
        vector<pair<vg::id_t, vg::id_t>> snarl_roots;
        // distance_index->is_root_snarl
        distance_index->traverse_decomposition([&] (const bdsg::net_handle_t& snarl) {
            
            if (distance_index->get_depth(snarl) != 1) // snarl must be top-level.
            {
                return true;
            }
            handle_t inward_source_handle = distance_index->get_handle(distance_index->get_bound(snarl, false, true), &*graph);
            handle_t outward_sink_handle = distance_index->get_handle(distance_index->get_bound(snarl, true, false), &*graph);
            // cerr << graph->get_id(inward_source_handle) << " " << graph->get_id(outward_sink_handle) << endl;
            // cerr << "distance_index->get_depth(snarl): " << distance_index->get_depth(snarl) << endl;
            // if (!distance_index->is_root_snarl(snarl))
            // {
            //     cerr << "not root snarl" << graph->get_id(inward_source_handle) << " " << graph->get_id(outward_sink_handle) << endl;
            //     return true;
            // } 

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

        cerr << "segregating regions for parallelization" << endl; //so that edits in two parallel jobs can't touch the same node.
        vg::algorithms::NormalizeRegionFinder region_finder = vg::algorithms::NormalizeRegionFinder(*graph, max_region_size, max_region_gap);

        parallel_regions_gbwt_updates = region_finder.get_parallel_normalize_regions(snarl_roots, *distance_index, parallel_normalize_regions, desegregation_candidates);
    }
    vg::io::save_handle_graph(graph.get(), std::cout);
    return 0;

    //     if (run_tests)
    //     {
    //         // cerr << "non-parallel normalize regions: " << endl;
    //         // for (auto region: snarl_roots)
    //         // {
    //         //     cerr << region.first << " " << region.second << endl;
    //         // }
    //         // cerr << "parallel normalize regions: " << endl;
    //         // for (auto region : parallel_normalize_regions)
    //         // {
    //         //     cerr << region.first << " " << region.second << endl;
    //         // }
    //         // cerr << " regions in parallel_regions_gbwt_updates: " << endl;
    //         // for (auto region : parallel_regions_gbwt_updates)
    //         // {
    //         //     cerr << "before: " << endl;
    //         //     for (auto node : region.first)
    //         //     {
    //         //         cerr << gbwt::Node::id(node) << " ";
    //         //     }
    //         //     cerr << endl;
    //         //     cerr << "after: " << endl;
    //         //     for (auto node : region.second)
    //         //     {
    //         //         cerr << gbwt::Node::id(node) << " ";
    //         //     }
    //         //     cerr << endl;
    //         // }
    //         // cerr << " cands in desegregation_candidates: " << endl;
    //         // for (auto cand : desegregation_candidates)
    //         // {
    //         //     cerr << "before: " << endl;
    //         //     cerr << cand.second << endl;
    //         //     cerr << "after: " << endl;
    //         //     cerr << cand.first.first << " " << cand.first.second << endl;
    //         // }
    //         assert(parallel_regions_gbwt_updates.size()==desegregation_candidates.size());
    //         for (int i = 0; i != parallel_regions_gbwt_updates.size(); i++)
    //         {
    //             //the "before" version of the update should be 1 in length:
    //             assert(parallel_regions_gbwt_updates[i].first.size()==1);
    //             //and equivalent to the original id in desegregation_candidates:
    //             assert(gbwt::Node::id(parallel_regions_gbwt_updates[i].first.back())==desegregation_candidates[i].second);
    //             //the "after" version of the update should be 2 in length:
    //             assert(parallel_regions_gbwt_updates[i].second.size()==2);
    //             //and equivalent to the new ids in desegregation_candidates:
    //             assert(gbwt::Node::id(parallel_regions_gbwt_updates[i].second.front())==desegregation_candidates[i].first.first);
    //             assert(gbwt::Node::id(parallel_regions_gbwt_updates[i].second.back())==desegregation_candidates[i].first.second);
    //         }
    //     }
    // }
    // else
    // {
        
    //     std::ifstream file( input_segregate_regions_only_file );
    //     string line_str;
    //     while(getline(file, line_str, '\n') && line_str!="-----")
    //     {
    //         stringstream line_ss(line_str);

    //         string region_first; string region_second;
    //         getline(line_ss, region_first, ' ');
    //         getline(line_ss, region_second, ' ');

    //         pair<vg::id_t, vg::id_t> region = make_pair(parse<int>(region_first), parse<int>(region_second));
    //         parallel_normalize_regions.push_back(region);
    //         // cerr << "working with region: " << region.first << " " << region.second << endl;
    //         // cerr << "in gbwtgraph: " << gbwt_graph->has_node(region.first) << " " << gbwt_graph->has_node(region.second) << endl;
    //         // cerr << "in graph: " << graph->has_node(region.first) << " " << graph->has_node(region.second) << endl;
    //     }
    //     while(getline(file, line_str, '\n'))
    //     {
    //         stringstream line_ss(line_str);
            
    //         string region_first; string region_second; string region_original;
    //         getline(line_ss, region_first, ' ');
    //         getline(line_ss, region_second, ' ');
    //         getline(line_ss, region_original, ' ');

    //         pair<vg::id_t, vg::id_t> desegregation_candidate = make_pair(parse<int>(region_first), parse<int>(region_second));
    //         pair< pair< vg::id_t, vg::id_t >, vg::id_t > candidate_with_original = make_pair(desegregation_candidate, parse<int>(region_original));

    //         desegregation_candidates.push_back(candidate_with_original);

    //     }
    //     file.close();
    // }



    // cerr << "updating the gbwt and gbwt_graph with the isolated regions" << endl;
    // auto isolated_regions_gbwt_update_start = chrono::high_resolution_clock::now();    

    // gbwt::GBWT parallel_regions_gbwt = vg::algorithms::apply_gbwt_changelog(*gbwt_graph, parallel_regions_gbwt_updates, *gbwt, gbwt_threads, false);

    // auto isolated_regions_gbwt_update_end = chrono::high_resolution_clock::now();    
    // chrono::duration<double> gbwt_changelog_time = isolated_regions_gbwt_update_end - isolated_regions_gbwt_update_start;
    // cerr << "Time spent generating the updated gbwt after isolating regions for parallelization: " << gbwt_changelog_time.count() << " s" << endl;

    // cerr << "generating the updated gbwt graph" << endl;
    // // make a new gbwt_graph for the parallel_regions_gbwt.
    // gbwtgraph::GBWTGraph parallel_regions_gbwt_graph;
    // if (parallel_regions_gbwt_updates.size()==0)
    // {
    //     parallel_regions_gbwt_graph = *gbwt_graph;
    // }
    // else
    // {
    //     parallel_regions_gbwt_graph = gbwtgraph::GBWTGraph(parallel_regions_gbwt, *graph);
    // }

    // // cerr << "output_segregate_regions_only_file.size() " << output_segregate_regions_only_file.size() << endl;
    // if (output_segregate_regions_only_file.size()>0)
    // {
    //     // cerr << "saving updated graph to file" << endl;
    //     //save normalized graph
    //     vg::io::save_handle_graph(graph.get(), std::cout);

    //     // cerr << "saving updated gbwt" << endl;
    //     save_gbwt(parallel_regions_gbwt, output_gbwt_file, true);

    //     // cerr << "saving extra segregate regions data to file " << output_segregate_regions_only_file << endl;
    //     std::ofstream output_file(output_segregate_regions_only_file);
    //     for (auto const& region : parallel_normalize_regions)
    //     {
    //         // cerr << "saving the following to output: " << region.first << " " << region.second << endl;
    //         output_file << region.first << " " << region.second << endl;
    //     }
    //     output_file << "-----" << endl; //this marks the end of parallel_normalize_regions and the beginning of desegregation_candidates.
    //     for (auto const& cand : desegregation_candidates)
    //     {
    //         // cerr << "saving the following to desegregation_candidates: " << endl;
    //         output_file << cand.first.first << " " << cand.first.second << " " << cand.second << endl;
    //     }
    //     exit(0);
        
    //     // std::ostream_iterator<std::string> output_iterator(output_file, "\n");
        
    //     // // std::ostream_iterator<std::pair<gbwtgraph::nid_t, gbwtgraph::nid_t>> output_iterator(output_file, "\n");
    //     // std::copy(std::begin(parallel_normalize_regions), std::end(parallel_normalize_regions), output_iterator);

    //     // for (auto const& x : vec1) 
    //     // {
    //     //     fout << x.first << ": "; 
    //     //     for (float f : x.second) fout << f << " ";
    //     //     fout << '\n';
    //     // }
        
        
    //     // exit(0);

    // }

    // if (run_tests)
    // {
    //     //test for all desegregation candidates having made it into the graph, gbwt_graph, and gbwt.
    //     for (auto cand : desegregation_candidates)
    //     {
    //         // cerr << "desegregation_candidates: " << cand.first.first << " and " << cand.first.second << endl;
    //         assert(parallel_regions_gbwt_graph.has_node(cand.first.first));
    //         assert(parallel_regions_gbwt_graph.has_node(cand.first.second));
    //         assert(graph->has_node(cand.first.first));
    //         assert(graph->has_node(cand.first.second));
    //         assert(parallel_regions_gbwt.contains(gbwt::Node::encode(cand.first.first, false)));
    //         assert(parallel_regions_gbwt.contains(gbwt::Node::encode(cand.first.second, false)));
            
    //     }

        
    // }

    // // free some memory, since we don't need the original gbwt anymore. Assuming that the parallel_regions_gbwt is actually new because there were updates made to it.
    // if (!parallel_regions_gbwt_updates.size()==0)
    // {
    //     gbwt.reset();
    // }

    // // =======running normalize=======
    // cerr << "=======running normalize=======" << endl;
    // string alignment_algorithm="sPOA"; bool disable_gbwt_update=false; //todo: remove these options or else ensure they are implemented correctly.

    // if (parallel_normalize_regions.size() == 0)
    // {
    //     cerr << "There appears to be no snarls in the graph. Ending normalize without returning a graph." << endl;
    //     exit(0);
    // }
    

    // //todo: note: debug code.
    // // std::set<vg::id_t> nodes_to_delete;
    // // vg::algorithms::SnarlNormalizer normalizer = vg::algorithms::SnarlNormalizer(
    // //   *graph, *gbwt, *gbwt_graph, nodes_to_delete, max_handle_size, max_region_size, threads, max_strings_per_alignment, "GBWT", alignment_algorithm, disable_gbwt_update, debug_print);
    // // std::vector<vg::RebuildJob::mapping_type> gbwt_normalize_updates = normalizer.parallel_normalization(snarl_roots);
    // vg::algorithms::SnarlNormalizer normalizer = vg::algorithms::SnarlNormalizer(
    //   *graph, parallel_regions_gbwt, parallel_regions_gbwt_graph, max_handle_size, max_region_size, threads, max_strings_per_alignment, "GBWT", alignment_algorithm, disable_gbwt_update, debug_print);

    // std::vector<vg::RebuildJob::mapping_type> gbwt_normalize_updates = normalizer.parallel_normalization(parallel_normalize_regions);

    // cerr << "=======updating gbwt after normalization and before desegregating regions=======" << endl;
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
    // vg::algorithms::NormalizeRegionFinder post_norm_region_finder = vg::algorithms::NormalizeRegionFinder(*graph, max_region_size, max_region_gap);
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


    // return 0;
}


// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize) ;

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

    // v2 distance index
    cerr << "loading distance index" << endl;
    auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_index_file);

    // =======isolating normalize regions for multithreading=======
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


    cerr << "isolating regions" << endl;
    vg::algorithms::NormalizeRegionFinder region_finder = vg::algorithms::NormalizeRegionFinder(*graph, max_snarls_per_region, max_snarl_spacing);
    std::set<vg::id_t> nodes_to_delete;
    std::vector<std::pair<vg::id_t, vg::id_t>> parallel_normalize_regions;
    std::vector<vg::RebuildJob::mapping_type> parallel_regions_gbwt_update_info = region_finder.get_parallel_normalize_regions(snarl_roots, parallel_normalize_regions, nodes_to_delete);

    cerr << "updating the gbwt and gbwt_graph with the isolated regions" << endl;
    auto gbwt_update_start = chrono::high_resolution_clock::now();    

    int gbwt_threads = threads;
    if (threads > 12)
    {
      gbwt_threads = 12;
    }

    gbwt::GBWT parallel_regions_gbwt = vg::algorithms::apply_gbwt_changelog(*gbwt_graph, parallel_regions_gbwt_update_info, *gbwt, gbwt_threads, false);
    
    auto gbwt_update_end = chrono::high_resolution_clock::now();    
    chrono::duration<double> gbwt_changelog_time = gbwt_update_end - gbwt_update_start;
    cerr << "Time spent generating the updated gbwt: " << gbwt_changelog_time.count() << " s" << endl;

    cerr << "generating the updated gbwt graph" << endl;
    // make a new gbwt_graph for the parallel_regions_gbwt.
    gbwtgraph::GBWTGraph parallel_regions_gbwt_graph = gbwtgraph::GBWTGraph(parallel_regions_gbwt, *graph);

    // free some memory, since we don't need the original gbwt anymore.
    gbwt.reset();

    //todo: delete debug code
    cerr << " parallel_regions_gbwt_graph.get_sequence(parallel_regions_gbwt_graph.get_handle(18, false)) " << parallel_regions_gbwt_graph.get_sequence(parallel_regions_gbwt_graph.get_handle(18, false)) << endl;

    // =======running normalize=======
    cerr << "=======running normalize=======" << endl;
    string alignment_algorithm="spoa"; bool disable_gbwt_update=false; bool debug_print=false; //todo: remove these options or else ensure they are implemented correctly.
    vg::algorithms::SnarlNormalizer normalizer = vg::algorithms::SnarlNormalizer(
      *graph, parallel_regions_gbwt, parallel_regions_gbwt_graph, nodes_to_delete, max_handle_size, max_snarls_per_region, max_snarl_spacing, threads, max_strings_per_alignment, "GBWT", alignment_algorithm, disable_gbwt_update, debug_print);

    std::vector<vg::RebuildJob::mapping_type> gbwt_normalize_updates = normalizer.parallel_normalization(parallel_normalize_regions);

    cerr << "======preparing and saving output=======" << endl;
    cerr << "saving updated graph" << endl;
    //save normalized graph
    vg::io::save_handle_graph(graph.get(), std::cout);

    cerr << "updating gbwt" << endl;
    auto _gbwt_update_start = chrono::high_resolution_clock::now();    
    gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(parallel_regions_gbwt_graph, gbwt_normalize_updates, parallel_regions_gbwt, gbwt_threads, false);

    cerr << "saving updated gbwt" << endl;
    save_gbwt(normalized_gbwt, output_gbwt_file, true);
    auto _gbwt_update_end = chrono::high_resolution_clock::now();    

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cerr << "Total elapsed time: " << elapsed.count() << " s" << endl;


    return 0;
}


// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize) ;
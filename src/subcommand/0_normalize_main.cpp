// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce
// more efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"

#include "subcommand.hpp"

// todo: should be able to remove '../../include/...' and replace with e.g.
// <bdsg/hash...>
#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"

#include "../io/save_handle_graph.hpp"

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

//todo: remove disable_gbwt_update option.

// gbwt::GBWT run_norm(vector<const Snarl *> snarl_roots, int optind, int argc, char** argv, string gbwt_file, string gbwt_graph_file, int max_handle_size, int max_alignment_size, int max_region_size, int max_snarl_spacing, int threads,  bool disable_gbwt_update, bool debug_print)
// // tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> run_norm(vector<const Snarl *> snarl_roots, int optind, int argc, char** argv, string gbwt_file, string gbwt_graph_file, int max_handle_size, int max_alignment_size, int max_region_size, int max_snarl_spacing, int threads,  bool disable_gbwt_update, bool debug_print)
// {
//   // getting graph of any type, except non-mutable graphs (e.g., xg)
//   shared_ptr<MutablePathDeletableHandleGraph> graph;
//   get_input_file(optind, argc, argv, [&](istream &in) {
//     graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
//   });


//     /// Build the gbwt:
//   ifstream gbwt_stream;
//   gbwt_stream.open(gbwt_file);

//   // Load the GBWT from its container
//   unique_ptr<gbwt::GBWT> gbwt;
//   gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
  
//   // gbwtgraph::GBWTGraph gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
//   // algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
//   // *graph, *gbwt, gbwt_graph, max_handle_size, max_alignment_size);
  
//   // gbwtgraph::GBWTGraph gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
//   // save_gbwtgraph(gbwt_graph, gbwt_file+".gg");
//   // unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph_p;
//   // gbwt_graph_p = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_file+".gg");
//   // algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
//   // *graph, *gbwt, *gbwt_graph_p, max_handle_size, max_alignment_size);

//   unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph;
//   if (gbwt_graph_file.size() == 0)
//   {
//     string gbwt_graph_output_file = gbwt_file + ".gg";
//     cerr << "gbwt_graph option is empty. Making new GBWTGraph from gbwt and graph. Saving as " << gbwt_graph_output_file << endl;
//     gbwtgraph::GBWTGraph new_gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
//     save_gbwtgraph(new_gbwt_graph, gbwt_graph_output_file);
//     //todo: find way to load gbwtgraph's unique pointer without saving and then reloading file.
//     gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_output_file);
//     gbwt_graph->set_gbwt(*gbwt);
//   }
//   else 
//   {
//     gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_file);
//     gbwt_graph->set_gbwt(*gbwt);
//   }


//   algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
//     *graph, *gbwt, *gbwt_graph, max_handle_size, max_region_size, max_snarl_spacing, threads, max_alignment_size, "GBWT", disable_gbwt_update, debug_print);

//   tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> gbwt_update_items = normalizer.normalize_snarls(snarl_roots);

//   //todo: integrate this into norm. And make the handlegraph into a non-shared pointer.
//   //todo: return gbwt_update ingredients, then call update outside this function.
//   // Save the modified graph
//   vg::io::save_handle_graph(graph.get(), std::cout);
  
//   // free memory:
//   if (!graph.unique()) //todo: verify this is never true
//   {
//     cerr <<" --------------------------------------error: graph ptr is not unique. Fix code." << endl;
//   }
//   graph.reset();
//   snarl_roots.clear();
  
//   if (!disable_gbwt_update)
//   {
//     gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(get<0>(gbwt_update_items), get<1>(gbwt_update_items), get<2>(gbwt_update_items), threads, debug_print);
//     return normalized_gbwt;
//   }
//   else
//   {
//     gbwt::GBWT empty_gbwt;
//     return empty_gbwt;
//   }
//   // gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(get<0>(gbwt_update_items), get<1>(gbwt_update_items), get<2>(gbwt_update_items), threads, debug_print);
  
//   // return gbwt_update_items;

//   // if (disable_gbwt_update) 
//   // {
//   //   gbwt::GBWT empty_gbwt;
//   //   return make_pair(graph, empty_gbwt);
//   // }
//   // else
//   // {
//   //   gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(get<0>(gbwt_update_items), get<1>(gbwt_update_items), get<2>(gbwt_update_items), threads, debug_print);
//   //   return make_pair(graph, normalized_gbwt);
//   // }
//   // gbwt::GBWT normalized_gbwt = normalizer.normalize_snarls(snarl_roots);
// }


// //binary search:
// void binary_search_norm(vector<const Snarl *> chosen_snarls, int snarl_start, int snarl_end, int optind, int argc, char** argv, string gbwt_file, string gbwt_graph, int max_handle_size, int max_alignment_size, int max_region_size, int max_snarl_spacing, int threads,  bool disable_gbwt_update, bool debug_print) // pass 0 and chosen_snarls.size() for first snarl_start/end. 
// {
//   cerr << "\nnormalizing snarls between: " << snarl_start << " and " << snarl_end << endl;
//   if (snarl_end - snarl_start  == 0)
//   {
//     return;
//   }
//   else
//   {
//     try
//     {
//       run_norm(chosen_snarls, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size, max_region_size, max_snarl_spacing, threads, disable_gbwt_update, debug_print);
//     } catch (const std::out_of_range& e) {
//       vector<const Snarl *>::const_iterator first = chosen_snarls.begin();
//       vector<const Snarl *>::const_iterator mid = chosen_snarls.begin() + chosen_snarls.size()/2;
//       vector<const Snarl *>::const_iterator last = chosen_snarls.end();
//       vector<const Snarl *> first_half(first, mid);
//       vector<const Snarl *> second_half(mid, last);
//       cerr << "normalizing snarls between: " << snarl_start << " and " << snarl_end << " failed. Splitting." << endl;
//       binary_search_norm(first_half, snarl_start, snarl_start + chosen_snarls.size()/2, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size, max_region_size, max_snarl_spacing, threads, disable_gbwt_update, debug_print);
//       binary_search_norm(second_half, snarl_start + chosen_snarls.size()/2, snarl_end, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size, max_region_size, max_snarl_spacing, threads, disable_gbwt_update, debug_print);
//     } catch (...) {
//       cerr << "snarls between: " << snarl_start << " and " << snarl_end << " ended successfully." << endl;
//     }
//   }
//   return;
// } 

void help_normalize(char **argv) {
  cerr
      << "usage: " << argv[0] << " normalize [options] <graph.vg> >[normalized.vg]"
      << endl
      << "Modifies snarls, outputs modified on stdout." << endl
      << endl
      << "options:" << endl
      << "    -g, --gbwt       gbwt corresponding to input graph." << endl
      << "    -r, --gbwt_graph       (optional) the gbwt_graph corresponding to the input graph and gbwt. If not supplied, will be temporarily generated." << endl
      << "    -s, --snarls       snarls file corresponding to hashgraph. Must include trivial snarls (i.e. made with vg snarls -T)."
      << endl
      << "    -o, --output_gbwt   name for the gbwt corresponding to the normalized graph. Default is normalized.gbwt."
      << "    -n, --normalize_region        Directly normalizes the graph between nodes -a and -b. Note: code assumes that the regions between a and b is only reachable from nodes within the region."
      << endl
      << "    -a, --leftmost_id        only used when --normalize_type is 'one', or when using -x option. "
         "The leftmost_id of the single node to be normalized. (for orientation of a snarl, check the 'backwards' field of the snarl.)"
      << endl
      << "    -b, --rightmost_id      only used when --normalize_type is 'one' or when using -x option. The "
         "rightmost_id of the single node to be normalized. (for orientation of a snarl, check the 'backwards' field of the snarl.)"
      << endl
      << "    -B, --normalize_bed_overlap      only normalizes snarls that are fully "
         "contained in at least one region in the provided bed file."
      << endl
      << "    -p, --paths_right_to_left       only used when --normalize_type "
         "is 'one'. This must be passed to the normalizer when the gbwt "
         "embedded paths are moving through the snarl from right to left, "
         "rather than left to right. This is frequently the case when the "
         "node_ids are larger on the left than on the right in the graph."
      << endl
      << "    -m, --max_alignment_size       limits the size of the snarl that "
         "will be normalized. If the number of gbwt paths through the snarl "
         "exceeds the max size, that snarl will be skipped for normalization. "
         "(Any additional embedded paths stretching from source to sink that "
         "aren't represented in the gbwt will also contribute to the max size.)"
      << endl
      << "    -z, --snarl_sizes       identifies the size of every top-level "
         "snarl of the inputted graph, outputs in a document specified with "
         "format 'source\tsink\tsize\n'. When passed this argument, there is "
         "no normalization."
      << endl
      << "    -j, --snarl_sizes_skip_source_sink      when passed in "
         "conjunction with --snarl_sizes/-i, the snarl size counting ignores "
         "the sequence in the source and sink. This avoids double-counting the "
         "sequence in a source/sink node between two adjacent snarls."
      << endl
      << "    -l, --max_handle_size       currently, default is INT_MAX. This "
         "is for compatibility with changes in graph size measures (file "
         "0_snarl_analyzer, arg '-i'). Eventually should change to handle size "
         "standard."
      << endl
      << "    -x, --handles_in_snarl      used in conjunction with arguments "
         "source and sink. Will print all the node ids in between source and "
         "sink, inclusive."
      << endl
      << "    -y, --max_search_dist      optional argument used in conjunction with "
         "argument handles_in_snarl. Determines how far in one direction (both left and right) the graph will search for the other node (a or b) before giving up. Default is 500*32 bases, or 500 standard handles."
         "source and sink. Will print all the node ids in between source and "
         "sink, inclusive."
      << endl
      << "    -c, --start_snarl_num      counting starting from 1, determines which snarls from the .snarls.pb will be normalized if normalize_type='all'."
      << endl
      << "    -v, --end_snarl_num      counting starting from 1, determines which snarls from the .snarls.pb will be normalized if normalize_type='all'."
      << endl
      << "    -k, --max_region_size      The max number of snarls to be realigned together at once, assuming they are close enough together to be reasonably realigned together. Default: 1. If specified, strongly consider changing max_snarl_spacing, too." //todo: change to somehow be a metric of the size of a snarl cluster to be aligned. Maybe max number of snarls to be aligned at once? No, something better than that...
      << endl
      << "    -i, --max_snarl_spacing      The max number of nodes between snarls permitted, before starting a new alignment region. Default:1." //todo: make it so that it's max number of *sequence* between snarls permitted, rather than the somewhat arbitrary node size of the graph.
      << endl
      << "    -t, --threads      The number of threads used in the normalization process. Currently only affects the update_gbwt step, which can grow very long without multithreading. Default:1."
      << endl
      << "    -E, --extract_paths      A debugging tool. Given a snarl, prints all the paths through the snarl to cout. Snarl is defined with -a and -b. Also requires -g. (-r saves time but isn't required)." //todo: include an extension of search for the "after normalization" phase which looks in the next snarl for possible sequence-shifting between snarls.
      << endl
      << "    -C, --compare_full_gbwt_paths      A debugging tool. Determines if two gbwts have identical haplotype content. Pass one of the gbwts to -g, and the other to this argument (-C). Optionally, to bypass remaking the gbwt_graphs, pass the gbwts corresponding gbwt_graphs to -r and -D, respectively."
      << endl
      // << "    -C_old, --compare_extracts      A debugging tool. Identifies if two outputs from two different extract_paths are comparable by performing a substring comparison. If all the strings in one file are substrings to the strings in the other file, or vice-versa, it asserts that the two snarls contain the same sequence information. Pass one -E output to this argument, and the other -E output to '-D'." //todo: include an extension of search for the "after normalization" phase which looks in the next snarl for possible sequence-shifting between snarls.
      // << endl
      // << "    -D, --compare_extracts_second_file      Use in conjuction with -C. Pass one of the output files from -E to -C, and the other to this argument." //todo: include an extension of search for the "after normalization" phase which looks in the next snarl for possible sequence-shifting between snarls.
      // << endl
      << "    -q, --disable_gbwt_update      skips the (sometimes lengthy) gbwt editing process after normalization. Note that there is no way to generate a new, haplotype-based gbwt based on a normalized graph without exporting the updated gbwt at this stage."
      << endl
      << "    -d, --debug_print      prints the location of every snarl/region being normalized. It also announces certain special cases - of particular note, when a snarl has increased in size after normalization."
      << endl
      << "    -h, --help      print this help info." << endl;
}

int main_normalize(int argc, char **argv) {

  if (argc == 2) {
    help_normalize(argv);
    return 1;
  }

  int max_alignment_size =
      INT_MAX; // default cutoff used to be 200 threads in a snarl.
  string gbwt_file;
  string gbwt_graph_file;
  string snarls;
  string output_gbwt_file = "normalized.gbwt";
  bool normalize_region = false;
  int leftmost_id = NULL; // todo: do something other than NULL to avoid the compiler
                     // warnings.
  int rightmost_id = NULL;
  string bed_file;
  bool paths_right_to_left = false;
  bool evaluate = false;
  string snarl_sizes;
  bool snarl_sizes_skip_source_sink = false;
  int max_handle_size = 32;
  bool handles_in_snarl = false;
  int max_search_dist = 500*32;
  //todo: note: options mostly for debugging:
  int start_snarl_num = 0;
  int end_snarl_num = 0;
  int max_region_size = 1;
  int max_snarl_spacing = 1;
  int threads = 1;
  bool extract_paths = false;
  string compare_full_gbwt_paths;
  // string extract_paths_a; //associated with -C_old argument
  // string extract_paths_b; //associated with -D argument
  bool disable_gbwt_update = false;
  bool debug_print = false;
  int c;
  string normalize_type = "all"; //note: only used internally; not directly set by an argument.
  optind = 2; // force optind past command positional argument
  while (true) {
    static struct option long_options[] =

        {{"help", no_argument, 0, 'h'},
         {"gbwt", required_argument, 0, 'g'},
         {"gbwt_graph", required_argument, 0, 'r'},
         {"snarls", required_argument, 0, 's'},
         {"output_gbwt", required_argument, 0, 'o'},
         {"normalize_region", required_argument, 0, 'n'},
         {"source", required_argument, 0, 'a'},
         {"sink", required_argument, 0, 'b'},
         {"normalize_bed_overlap", required_argument, 0, 'B'},
         {"paths_right_to_left", no_argument, 0, 'p'},
         {"max_alignment_size", required_argument, 0, 'm'},
         {"snarl_sizes", required_argument, 0, 'z'},
         {"snarl_sizes_skip_source_sink", no_argument, 0, 'j'},
         {"max_handle_size", required_argument, 0, 'l'},
         {"handles_in_snarl", no_argument, 0, 'x'},
         {"max_search_dist", no_argument, 0, 'y'},
         {"start_snarl_num", required_argument, 0, 'c'},
         {"end_snarl_num", required_argument, 0, 'v'},
         {"max_region_size", required_argument, 0, 'k'},
         {"max_snarl_spacing", required_argument, 0, 'i'},
         {"threads", required_argument, 0, 't'},
         {"extract_paths", no_argument, 0, 'E'},
         {"compare_full_gbwt_paths", required_argument, 0, 'C'},
        //  {"compare_extracts", required_argument, 0, 'C_old'},
        //  {"compare_extracts_second_file", required_argument, 0, 'D'},
         {"disable_gbwt_update", no_argument, 0, 'q'},
         {"debug_print", no_argument, 0, 'd'},
         {0, 0, 0, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "hg:r:s:o:na:b:B:pm:i:jl:xy:c:v:k:i:t:EC:qd", long_options,
                    &option_index);

    // Detect the end of the options.
    if (c == -1)
      break;

    switch (c) {

    case 'g':
      gbwt_file = optarg;
      break;

    case 'r':
      gbwt_graph_file = optarg;
      break;

    case 's':
      snarls = optarg;
      break;

    case 'o':
      output_gbwt_file = optarg;
      break;

    case 'n':
      normalize_region = true;
      break;

    case 'a':
      leftmost_id = parse<int>(optarg);
      break;

    case 'b':
      rightmost_id = parse<int>(optarg);
      break;

    case 'B':
      bed_file = optarg;
      break;

    case 'p':
      paths_right_to_left = true;
      break;

    case 'm':
      max_alignment_size = parse<int>(optarg);
      // if max_alignment_size is 0, then that signifies that it should actually
      // be infinite, i.e. that we should not exclude any snarls.
      if (max_alignment_size == 0) {
        max_alignment_size = INT_MAX;
      }
      break;

    case 'z':
      snarl_sizes = optarg;
      // normalize_type = "none";
      break;

    case 'j':
      snarl_sizes_skip_source_sink = true;
      break;

    case 'l':
      max_handle_size = parse<int>(optarg);
      break;

    case 'x':
      handles_in_snarl = true;
      normalize_type = "none";
      break;

    case 'y':
      max_search_dist = parse<int>(optarg);
      break;

    case 'c':
      start_snarl_num = parse<int>(optarg);
      break;

    case 'v':
      end_snarl_num = parse<int>(optarg);
      break;

    case 'k':
      max_region_size = parse<int>(optarg);
      break;

    case 'i':
      max_snarl_spacing = parse<int>(optarg);
      break;
      
    case 't':
      threads = parse<int>(optarg);
      break;
      
    case 'E':
      extract_paths = true;
      normalize_type = "none";
      break;
      
    case 'C':
      compare_full_gbwt_paths = optarg;
      normalize_type = "none";
      break;

    // case 'C_old':
    //   extract_paths_a = optarg;
    //   normalize_type = "none";
    //   break;
      
    // case 'D':
    //   extract_paths_b = optarg;
    //   normalize_type = "none";
    //   break;

    case 'q':
      disable_gbwt_update = true;
      break;

    case 'd':
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


  
  if (normalize_type != "all" && normalize_type != "none") 
  {
    cerr << "please enter a valid normalize_type: all or none." << endl;
  }

  if (normalize_type == "all") 
  {
      
    // cerr << "snarl_roots.size()" << snarl_roots.size() << endl;

    //todo: use for debugging:
    /// for select snarl:
    // vector<const Snarl *> chosen_snarls(snarl_roots.begin() + 1370135, snarl_roots.begin() + 1370136); 
    // run_norm(chosen_snarls, optind, argc, argv, gbwt, max_handle_size, max_alignment_size);
    /// for binary search:
    // binary_search_norm(snarl_roots, 0, snarl_roots.size(), optind, argc, argv, gbwt, max_handle_size, max_alignment_size, max_region_size, max_snarl_spacing, threads, disable_gbwt_update, debug_print);

    
    cerr << "running normalize!" << endl;

    // Record start time
    auto start = chrono::high_resolution_clock::now();
    pair<shared_ptr<MutablePathDeletableHandleGraph>, gbwt::GBWT> output;

    // unique_ptr<MutablePathDeletableHandleGraph> graph;
    // gbwt::GBWT normalized_gbwt;
    // pair<MutablePathDeletableHandleGraph, gbwt::GBWT> output = make_pair(*graph, normalized_gbwt);
    

    ///////// Input objects: ///////////
    // graph
    cerr << "messing with graph pointers here!" << endl;
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    // gbwt
    ifstream gbwt_stream;
    gbwt_stream.open(gbwt_file);    
    unique_ptr<gbwt::GBWT> gbwt;
    gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
 
    // gbwt graph 
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

    // snarl_stream:
    std::ifstream snarl_stream;
    string snarl_file = snarls;
    snarl_stream.open(snarl_file);
    if (!snarl_stream) {
      cerr << "error:[vg normalize] Cannot open Snarls file " << snarl_file
            << endl;
      exit(1);
    }
    cerr << "before snarl manager is created" << endl;
    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);

    // snarl_roots 
    // Depending on user specifications, we may want to normalize some specified snarls:
    vector<const Snarl *> snarl_roots;

    if (start_snarl_num != 0 || end_snarl_num != 0)
    {
      cerr << "normalizing user-specified snarls" << endl;
      //prepare to normalize user-specified snarls:
      //check that start and end are valid:
      if (start_snarl_num >= end_snarl_num)
      {
        cerr << "error:[vg normalize] start_snarl_num >= end_snarl_num."  << endl;
        exit(1);
      }
      else if (end_snarl_num > snarl_roots.size())
      {
        cerr << "WARNING:[vg normalize] end_snarl_num greater than snarl_roots.size(). Will normalize starting at start_snarl_num and ending at last snarl in snarl_roots." << endl;
      }
      
      snarl_roots = snarl_manager->top_level_snarls();
      vector<const Snarl *> chosen_snarls(snarl_roots.begin() + start_snarl_num, snarl_roots.begin() + end_snarl_num); 
      //todo: check that the following line properly overwrites snarl_roots with a copy of chosen_snarls. 
      snarl_roots = chosen_snarls;
      cerr << "of snarl selection that is " << chosen_snarls.size() << " long, first snarl selected has source: " << (*chosen_snarls.front()).start().node_id() << " and sink: " << (*chosen_snarls.front()).end().node_id() << endl; 
    }
    // or only normalize snarls that are fully contained by at least one region of a bed file:
    else if (bed_file.size() > 0)
    {
      cerr << "normalizing bed specified snarls" << endl;
      //todo: add bedfile functionality.
      bdsg::PositionOverlay position_graph = bdsg::PositionOverlay(graph.get());
      vector<Region> bed_regions;
      //todo: non-hardcode bed region
      parse_bed_regions(bed_file, bed_regions);
      // cerr << "Number of regions in the bedfile we will use to filter out undesired snarls: " << bed_regions.size() << " " << bed_regions.front().start << " " << bed_regions.front().end << " " << bed_regions.front().seq << " " << endl;
      cerr << "Number of regions in the bedfile we will use to filter out undesired snarls: " << bed_regions.size() << endl;
      //todo: check that contain_endpoints bool does what I think it does. I think I want contain_endpoints=true.
      visit_contained_snarls(&position_graph, bed_regions, *snarl_manager, true, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step, int64_t start_int, int64_t end_int, bool, const Region* region)
      {
        snarl_roots.push_back(snarl);
      });
    }
    // otherwise, normalize all snarls: (except if we have the region specified by -a and -b, in which case leave snarl_roots empty.)
    else if (!normalize_region)
    {
      cerr << "normalizing all snarls" << endl;
      snarl_roots = snarl_manager->top_level_snarls();
    }
    cerr << "number of snarls in targeted normalize regions: " << snarl_roots.size() << endl;
    

    /////////// normalize graph /////////////////
    // normalize
    cerr << "running normalize" << endl;
    algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
      *graph, *gbwt, *gbwt_graph, max_handle_size, max_region_size, max_snarl_spacing, threads, max_alignment_size, "GBWT", disable_gbwt_update, debug_print);
    if (normalize_region)
    {
      // if we've specified the normalize_region with -a and -b, then we just normalize that region. //todo: add gbwt update?
      normalizer.normalize_snarl(leftmost_id, rightmost_id, false, 0);
    }
    else
    {
      // perform the normalization.
      tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> gbwt_update_items = normalizer.normalize_snarls(snarl_roots);

      //save normalized graph
      //todo: integrate this into norm
      vg::io::save_handle_graph(graph.get(), std::cout);

      /////////// perform apply_gbwt_changelog /////////////////
      // free memory:
      if (!graph.unique()) //todo: verify this is never true
      {
        cerr <<" --------------------------------------error: graph ptr is not unique. Fix code." << endl;
        exit(1);
      }
      graph.reset();
      snarl_roots.clear();

      // apply_gbwt_changelog
      if (!disable_gbwt_update)
      {
        cerr << "applying changes from normalization to the gbwt." << endl;
        gbwt::GBWT normalized_gbwt = vg::algorithms::apply_gbwt_changelog(get<0>(gbwt_update_items), get<1>(gbwt_update_items), get<2>(gbwt_update_items), threads, debug_print);
        save_gbwt(normalized_gbwt, output_gbwt_file, true);
      }
    }

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cerr << "Elapsed time: " << elapsed.count() << " s" << endl;



    // // Save the modified graph
    // vg::io::save_handle_graph(output.first.get(), std::cout);
  }

  // snarl_analyzer identifies the size of every top-level snarl, outputs in a
  // document specified with format "source\tsink\tsize\n"
  if (snarl_sizes.size() != 0) {
    cerr << "messing with graph pointers here!" << endl;
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });
    
    std::ifstream snarl_stream;
    snarl_stream.open(snarls);
    if (!snarl_stream) {
      cerr << "error:[vg normalize] Cannot open Snarls file " << snarls << endl;
      exit(1);
    }

    algorithms::SnarlAnalyzer sizes = algorithms::SnarlAnalyzer(
        *graph, snarl_stream, snarl_sizes_skip_source_sink);

    sizes.output_snarl_sizes(snarl_sizes);
  }

  if (handles_in_snarl) 
  {
    if (leftmost_id == NULL && rightmost_id == NULL) 
    {
      cerr << "error:[vg normalize] please enter a values for leftmost_id and rightmost_id "
              "to define the snarl."
           << endl;
    } 
    else 
    {
      cerr << "messing with graph pointers here!" << endl;
      shared_ptr<MutablePathDeletableHandleGraph> graph;
      get_input_file(optind, argc, argv, [&](istream &in) {
        graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
      });
      algorithms::print_handles_in_snarl(*graph, leftmost_id, rightmost_id, max_search_dist);
    }
  }

  // debugging tool to find all haplotypes in a snarl:
  if (extract_paths)
  {
    ///////// Input objects: ///////////
    // graph
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    // gbwt
    ifstream gbwt_stream;
    gbwt_stream.open(gbwt_file);    
    unique_ptr<gbwt::GBWT> gbwt;
    gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
 
    // gbwt graph 
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

    // snarl_roots 
    // Depending on user specifications, we may want to normalize some specified snarls:
    vector<const Snarl *> snarl_roots;


    SubHandleGraph snarl = algorithms::SnarlNormalizer::extract_subgraph(*graph, leftmost_id, rightmost_id);

    // gbwt_haps is of format:
    // 0: a set of all the haplotypes which stretch from source to sink, in string format.
    //   - it's a set, so doesn't contain duplicates
    // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
    // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
    algorithms::SnarlSequenceFinder sequence_finder = algorithms::SnarlSequenceFinder(*graph, snarl, *gbwt_graph, leftmost_id, rightmost_id, false);
    tuple<vector<vector<gbwtgraph::handle_t>>, vector<vector<gbwtgraph::handle_t>>, unordered_set<handlegraph::handle_t>> gbwt_haps = sequence_finder.find_gbwt_haps();
    // cerr << get<1>(gbwt_haps).size() << " " << get<2>(gbwt_haps).size() << " " << get<3>(gbwt_haps).size() << " " << endl;
    cerr << " " << get<0>(gbwt_haps).size() << " " <<  get<1>(gbwt_haps).size() << " " << get<2>(gbwt_haps).size()  << endl;

    cerr << "about to print strings:" << endl;
    unordered_set<string> hap_strings = algorithms::SnarlNormalizer::format_handle_haplotypes_to_strings(*graph, *gbwt_graph, get<0>(gbwt_haps));
    for (auto string : hap_strings)
    {
      cerr << string << endl;
    }
  }
  if (compare_full_gbwt_paths.size() != 0) 
  {
    ///////// Input objects: ///////////
    // graph
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    cerr << "hi" << endl;
    // gbwt
    ifstream gbwt_stream;
    gbwt_stream.open(gbwt_file);    
    unique_ptr<gbwt::GBWT> gbwt;
    gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);

    // gbwt graph 
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

    // use min_node_id as the starting point. Get every path on that node, and go to their path_begin spot. Extract each of their paths.
    // Every time we traverse a new node, check to make sure that all paths on that node are not new. If there is a new path there, put its path.begin() in the list of paths to extend.
    //todo: ask lab if I can skip the "check each node if its new" thing, so long as I guarantee that the first node I check for paths is a top level snarl root. 
    //todo:     I.e.: "are gbwt_haplotypes all guaranteed to stretch from source to sink? Or can they start in the middle of a graph?"
    // gbwt_graph->for_each_path_handle(gbwt_graph->get_handle(gbwt_graph->min_node_id()), [&](path_handle_t path_handle) 
    vector<string> path_strings;
    cerr << "reality check: does gbwt_graph have nodes? Min node: " << gbwt_graph->min_node_id() << " max: " << gbwt_graph->max_node_id() << endl;
    handle_t test = gbwt_graph->get_handle(gbwt_graph->max_node_id());
    
    SubHandleGraph empty_snarl(&(*graph));
    // gbwt_haps is of format:
    // 0: a set of all the haplotypes which stretch from source to sink, in string format.
    //   - it's a set, so doesn't contain duplicates
    // 1: a vector of all the other haps in the snarl (in vector<handle_t> format)
    // 2: a vector of all the handles ever touched by the SnarlSequenceFinder.
    algorithms::SnarlSequenceFinder sequence_finder = algorithms::SnarlSequenceFinder(*graph, empty_snarl, *gbwt_graph, gbwt_graph->min_node_id(), gbwt_graph->max_node_id(), false);
    tuple<vector<vector<gbwtgraph::handle_t>>, vector<vector<gbwtgraph::handle_t>>, unordered_set<handlegraph::handle_t>> gbwt_haps = sequence_finder.find_gbwt_haps();

    cerr << "size of haps 0, 1, and 2: " << get<0>(gbwt_haps).size() << " " << get<1>(gbwt_haps).size() << " " << get<2>(gbwt_haps).size() << endl;
    cerr << "size of full length haps:" << endl;
    for (auto hap : get<0>(gbwt_haps))
    {
      cerr << hap.size() << endl;
    }
    cerr << "size of partial length haps:" << endl;
    for (auto hap : get<1>(gbwt_haps))
    {
      cerr << hap.size() << endl;
    }
    cerr << "to the right: " << endl;
    gbwt_graph->follow_edges(test, false, [&](handle_t next_handle){
      cerr << "handle to the right of id: " << gbwt_graph->get_id(next_handle) << endl; 
    });
    cerr << "to the left: " << endl;
    gbwt_graph->follow_edges(test, true, [&](handle_t next_handle){
      cerr << "handle to the left of id: " << gbwt_graph->get_id(next_handle) << endl; 
    });

    vector<step_handle_t> steps = gbwt_graph->steps_of_handle(gbwt_graph->get_handle(998743));
    cerr << "example handle's sequence" << gbwt_graph->get_sequence(gbwt_graph->get_handle(998743)) << endl;
    cerr << "all steps' paths: " << endl;
    for (auto step : steps)
    {
      cerr << "b" << endl;
      cerr << gbwt_graph->get_path_name(gbwt_graph->get_path_handle_of_step(step)) << endl;
    }

    cerr << "about to use for_each_path_handle:" << endl;
    gbwt_graph->for_each_path_handle([&](path_handle_t path_handle) 
    {
      cerr << "looking at path_handle " << gbwt_graph->get_path_name(path_handle) << endl;
      step_handle_t cur_step = gbwt_graph->path_begin(path_handle);
      string path_string;
      while (cur_step != gbwt_graph->path_end(path_handle))
      {
        path_string += gbwt_graph->get_sequence(gbwt_graph->get_handle_of_step(cur_step));
        cur_step = gbwt_graph->get_next_step(cur_step);
      }
      path_strings.insert(upper_bound(path_strings.begin(), path_strings.end(), path_string, [&] (const string& first, const string& second)
      {
        return (first.size() < second.size());
      }), path_string);
    });
    
    cerr << "string lengths: " << endl;
    for (auto path : path_strings)
    {
      cerr << path.size() << endl;
    }
    

  }

    
              // vector<pair<step_handle_t, step_handle_t>> embedded_paths = sequence_finder.find_gbwt_haps();
              // for (auto path : embedded_paths)
              // {
              //   step_handle_t cur_step = path.first;
              //   step_handle_t last_step = path.second;
              //   // cerr << "getting id of cur_step: " << (*graph).get_id((*graph).get_handle_of_step(cur_step)) << endl;
              //   // cerr << "getting handle of final...: " << endl;
              //   // (*graph).get_handle_of_step(path.second);
              //   // cerr << "getting id of final step: " << (*graph).get_id((*graph).get_handle_of_step((*graph).get_previous_step(last_step))) << endl;
              //   while (cur_step != last_step)
              //   {
              //     cout << (*graph).get_sequence((*graph).get_handle_of_step(cur_step)); //todo: check that if the handle's .reverse() is true, then does it include the proper reversed comp of the sequence?
              //     cur_step = (*graph).get_next_step(cur_step);
              //     // cerr << "(In while loop: " << endl;
              //     // cerr << "getting id of cur_step: " << (*graph).get_id((*graph).get_handle_of_step(cur_step)) << endl;
              //     // cerr << "getting id of final step: " << (*graph).get_id((*graph).get_handle_of_step(last_step)) << endl;
              //     // cerr << ")" << endl;
              //   }
              //   cout << endl;
              // }


  // if (extract_paths_a.size() != 0 && extract_paths_b.size() != 0)
  // {
  //   vector<string> paths_a;
  //   ifstream a_file(extract_paths_a);
  //   // check stream status
  //   if (!a_file)
  //   {
  //     cerr << "Can't open input file " << extract_paths_a << "." << endl;
  //     exit(1);
  //   }

  //   string line;
  //   while (getline(a_file, line)) {

  //     // store each line in the vector
  //     paths_a.push_back(line);
  //   }
  //   a_file.close();

  //   vector<string> paths_b;
  //   ifstream extract_paths_a_file(extract_paths_a);
  //   //if all the strings in a have a corresponding string in b, then we consider the two files to have functionally identical sequence content.
  //   //we determine "corresponding string" as a unique string in a that matches as a substring to a unique string in b, or vice-versa.
  //   //To prevent multiple small paths in one file matching to the same large string in the other file, I will preferentially match the string in a 

  //   //I just realized this still won't work, because there's no requirement that strings in a all be substrings to strings in b, or vice-versa. Normalization might increase the source but shrink the sink at the same time - leading to merely partially-overlapping paths extracted from either snarl.
  //   //BUT, the substring relation is restored if I drop the source and sink of the paths in the extract_paths algorithm.

  //   //so:
  //   //todo: drop source and sink in extract_paths.
  //   //todo: determine that all strings in a are substrings to string in b. If that fails, do the reverse (strings in b are substrings to a).
  //     // Note that every string-substring match I make should then be removed from the list of possible matches I make.
  //     // This is because each path in the original snarl should have its own designated path in the normalized snarl.
  //     // To avoid popping the wrong string, I should preferentially pop the smallest string that has the desired substring.

    
  // }

  return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize) ;
#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "../vg.hpp"
#include <string>
#include <gbwtgraph/gbwtgraph.h>
#include <spoa/spoa/spoa.hpp>
#include "bdsg/hash_graph.hpp"

namespace vg {
namespace algorithms {

class SnarlNormalizer {
  public:
    virtual ~SnarlNormalizer() = default;

    SnarlNormalizer(MutablePathDeletableHandleGraph &graph, const gbwt::GBWT &gbwt, const gbwtgraph::GBWTGraph & gbwt_graph,
                    const int max_handle_size, 
                    const int max_region_size,
                    const int max_snarl_spacing,
                    const int threads,
                    const int max_alignment_size = INT_MAX, //TODO: add a _max_handle_length default length
                    const string path_finder = "GBWT", /*alternative is "exhaustive"*/
                    const string alignment_algorithm = "sPOA",
                    const bool disable_gbwt_update = false,
                    const bool debug_print = false);


    virtual tuple<gbwtgraph::GBWTGraph, std::vector<vg::RebuildJob::mapping_type>, gbwt::GBWT> normalize_snarls(const vector<const Snarl *>& snarl_roots);

    virtual vector<int> normalize_snarl(const id_t source_id, const id_t sink_id, const bool backwards, const int snarl_num);

    static SubHandleGraph extract_subgraph(const HandleGraph &graph, const id_t leftmost_id, const id_t rightmost_id);

    virtual void output_msa(const id_t leftmost_id, const id_t rightmost_id);

    //////////////////////////////////////////////////////////////////////////////////////
    // format-type switching:
    //////////////////////////////////////////////////////////////////////////////////////
    static unordered_set<string> format_handle_haplotypes_to_strings(const HandleGraph& graph, const gbwtgraph::GBWTGraph & gbwt_graph,
        const vector<vector<handle_t>> &haplotype_handle_vectors);

  protected:
    // member variables:
    // the handle graph with snarls to normalize
    MutablePathDeletableHandleGraph &_graph;
    // GBWT graph with snarls to normalize, includes the embedded threads needed for the
    // GBWTPathFinder approach.
    const gbwt::GBWT &_gbwt;
    const gbwtgraph::GBWTGraph &_gbwt_graph;
    // const gbwtgraph::GBWTGraph _gbwt_graph = gbwtgraph::GBWTGraph(_gbwt, _graph);
    // const gbwtgraph::GBWTGraph &_gbwt_graph = gbwtgraph::GBWTGraph(_gbwt, _graph);
    // vector<pair<vector<std::uint32_t>, vector<std::uint32_t>>> _gbwt_changelog; //todo: delete this
    // vector<pair<vector<id_t>, vector<id_t>>> _gbwt_changelog;
    vector<pair<gbwt::vector_type, gbwt::vector_type>> _gbwt_changelog;

    // _touched_border_nodes tracks which nodes have been already considered when 
    // calculating the before-normalization total snarl size.
    unordered_set<id_t> _touched_border_nodes;
    int _pre_norm_net_snarl_size = 0;
    //for tracking handles that exist in the graph, but have no corresponding path in the gbwt.
    int _handles_not_touched_by_gbwt = 0;
    int _sequence_not_touched_by_gbwt = 0;
    // int _snarls_skipped_because_gbwt_misses_handles = 0; //currently handled by messy system of error_record vectors. See: error_record[7] in normalize_snarls(). //todo: change?
    int _skipped_snarl_sizes = 0;
    int _unskipped_snarl_sizes = 0;
    int _unskipped_snarl_num = 0;
    unordered_set<pair<id_t, id_t>> _skipped_snarls;
    unordered_set<pair<id_t, id_t>> _unskipped_snarls;
//     vector<int> skipped_snarl_sizes.push_back(snarl_size); //todo: remove this later for efficiency? Or at least turn into a rolling calculation of averages. (just a rolling sum + a tracker of total number skipped).
//     
    // vector<int> unskipped_snarl_sizes.push_back(snarl_size); //todo: remove for increased efficiency? Or at least turn into a rolling calculation of averages. (just a rolling sum + a tracker of total number skipped). 

    // a separate tracker for measuring independent snarl increase/decrease in size. Note 
    // that because snarls overlap, the sum of all the independent snarl 
    // increase/decreases will not equal the total increase/decrease in graph size.
    // Also note that the "after normalization" size is measured immediately after normalization.
    // Normalization of neighboring snarls may result in further changes of snarl size as 
    // a result of changes in the border handles where snarls overlap.
    // format: key: pair of leftmost,rightmost ids. value: size of snarl before,after normalization. 
    map<pair<id_t, id_t>, pair<int, int>> _snarl_size_changes;
    vector<int> _alignments_calling_for_abpoa; //todo: delete this after I've implemented abpoa?

    // the maximum number of haplotypes allowed to align in a given snarl. If the number of
    // haplotypes exceeds this threshold, the snarl is skipped.
    const int _max_alignment_size;
    const int _max_handle_size;
    const int _max_region_size;
    const int _max_snarl_spacing;
    const string _path_finder;
    const bool _disable_gbwt_update;
    const bool _debug_print; // for printing info that isn't necessarily something gone wrong.
    const int _threads;
    const string _alignment_algorithm;
    //////////////////////////////////////////////////////////////////////////////////////
    // finding information on original graph:
    //////////////////////////////////////////////////////////////////////////////////////

    // SubHandleGraph extract_subgraph(const HandleGraph &graph, const id_t leftmost_id, const id_t rightmost_id);
                                    
    vector<int> check_handle_as_start_of_path_seq(const string &handle_seq,
                                                  const string &path_seq);


    // clustering snarls based on _max_region_size and //todo add variable!
    vector<pair<id_t, id_t>> get_normalize_regions(const vector<const Snarl *>& snarl_roots);

    vector<pair<id_t, id_t>> get_single_snarl_normalize_regions(const vector<const Snarl *> &snarl_roots);

    bool is_trivial(const Snarl& snarl);

    bool snarls_adjacent(const Snarl& snarl_1, const Snarl& snarl_2);

    vector<vector<const Snarl *> > cluster_snarls(const vector<const Snarl *> &snarl_roots);

    vector<pair<id_t, id_t>> convert_snarl_clusters_to_regions(const vector<vector<const Snarl *> >& clusters);


    //////////////////////////////////////////////////////////////////////////////////////
    // creation of new graph:
    //////////////////////////////////////////////////////////////////////////////////////
    VG poa_source_to_sink_haplotypes(const unordered_set<string>& source_to_sink_haplotypes, const int snarl_num, const bool output_msa=false);

    VG kalign_source_to_sink_haplotypes(const unordered_set<string>& source_to_sink_haplotypes);

    VG align_source_to_sink_haplotypes(const unordered_set<string>& source_to_sink_haplotypes);

    pair<handle_t, handle_t> integrate_snarl(SubHandleGraph &old_snarl, const HandleGraph &new_snarl,
                         vector<pair<step_handle_t, step_handle_t>>& embedded_paths,
                         const id_t source_id, const id_t sink_id, const bool backwards);

    handle_t overwrite_node_id(const id_t old_node_id, const id_t new_node_id);

    void log_gbwt_changes(const vector<pair<vector<gbwt::vector_type::value_type>, string>>& source_to_sink_gbwt_paths, const HandleGraph &new_snarl);

    bool source_and_sink_handles_map_properly(
        const HandleGraph &graph, const id_t new_source_id, const id_t new_sink_id,
        const bool touching_source, const bool touching_sink,
        const handle_t potential_source, const handle_t potential_sink);

    void force_maximum_handle_size(MutableHandleGraph &graph);

    // moving paths to new graph (new draft functions)
    vector<pair<vector<handle_t>, int> > find_possible_path_starts (const handle_t leftmost_handle, const handle_t rightmost_handle, const pair<bool, bool> path_spans_left_right);

    vector<handle_t> extend_possible_paths(vector<pair<vector<handle_t>, int>> &possible_path_starts, const string &path_str, const handle_t leftmost_handle, const handle_t rightmost_handle, const pair<bool, bool> path_spans_left_right, const pair<id_t, id_t> main_graph_source_and_sink);

    pair<step_handle_t, step_handle_t> move_path_to_new_snarl(const pair<step_handle_t, step_handle_t> old_path, const id_t source, const id_t sink, const pair<bool, bool> path_spans_left_right, const bool path_directed_left_to_right, const pair<id_t, id_t> main_graph_source_and_sink);

    // updating the gbwt:
    gbwt::GBWT apply_gbwt_changelog();

    std::unordered_map<nid_t, size_t> get_node_to_job(const vector<unordered_set<nid_t>>& weakly_connected_components);

    std::vector<RebuildJob> divide_changelog_into_jobs(const std::unordered_map<nid_t, size_t>& node_to_job, const vector<unordered_set<nid_t>>& weakly_connected_components);

    RebuildParameters set_parameters();
    





    // -------------------------------- DEBUG CODE BELOW:
    // ------------------------------------

    pair<vector<handle_t>, vector<handle_t>>
    debug_get_sources_and_sinks(const HandleGraph &graph);

    void make_one_edit(id_t leftmost_id, id_t rightmost_id);
    void get_all_gbwt_sequences(id_t source, id_t sink_id, bool backwards);

    
};
}
} // namespace vg

#include "../handle.hpp"
#include "../gbwt_helper.hpp"
#include <gbwtgraph/gbwtgraph.h>
#include "../subgraph.hpp"
#include "../snarl_distance_index.hpp"

#include "../vg.hpp"
namespace vg {
namespace algorithms {

class NormalizeRegionFinder {
  public:
    virtual ~NormalizeRegionFinder() = default;

    NormalizeRegionFinder(MutablePathDeletableHandleGraph &graph,
                                 const int max_region_size,
                                 const int max_region_gap);

    /// Runs get_normalize_regions, and then split_sources_and_sinks on the normalize regions.
    // to ensure that each region is non-overlapping with each other region.
    // Also provides a list of nodes that need to be recombined into a single node after 
    // normalization, since some nodes that are shared between two snarls are one base 
    // long, but still need to be divided into two in order to guarantee that the two 
    // snarls aren't overlapping. 
    // arguments parallel_normalize_regions and nodes_to_remove are empty vectors 
    // that are filled by the function.
    // Returns the changes that need to be made ot the gbwt to account for split nodes.
    std::vector<vg::RebuildJob::mapping_type> get_parallel_normalize_regions(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots, const SnarlDistanceIndex& distance_index, vector<pair<id_t, id_t>>& parallel_normalize_regions, vector< pair< pair< id_t, id_t >, id_t > >& desegregation_candidates);

    // /// function called by get_parallel_normalize_regions:
    // vector<pair<id_t, id_t>> get_normalize_regions(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots);

    ///Function used to remove the splits introduced in split_sources_and_sinks. (NOTE: this will overwrite the current _gbwt_changelog.
    std::vector<vg::RebuildJob::mapping_type> desegregate_nodes(std::vector<std::pair<std::pair<id_t, id_t>, id_t>> desegregation_candidates);

  protected:

    //member variables from construction:
    MutablePathDeletableHandleGraph &_graph;
    const int _max_region_size;
    const int _max_region_gap;

    //other member variables:
    vector<vg::RebuildJob::mapping_type> _gbwt_changelog; //this is equivalent to vector<pair<gbwt::vector_type, gbwt::vector_type>> _gbwt_changelog;

    /// other functionS called by get_parallel_normalize_regions:
    vector<pair<id_t, id_t>> split_sources_and_sinks(vector<pair<id_t, id_t>> normalize_regions, vector< pair< pair< id_t, id_t >, id_t > >& desegregation_candidates);

    vector<pair<id_t, id_t>> cluster_snarls(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots, const SnarlDistanceIndex& distance_index);

    //used by cluster_roots:
    SubHandleGraph extract_subgraph(const HandleGraph &graph,
                                                 const id_t leftmost_id,
                                                 const id_t rightmost_id);


    bool test_snarl_for_clustering(const HandleGraph &graph, const int snarl_size);

    // /// functions called by get_normalize_regions.
    // vector< vector<pair<vg::id_t, vg::id_t>> > cluster_snarls(const vector<pair<vg::id_t, vg::id_t>> &snarl_roots);

    // vector<pair<id_t, id_t>> convert_snarl_clusters_to_regions(const vector<vector<pair<vg::id_t, vg::id_t>> >& clusters);

    // /// Functions used in cluster_snarls:
    // // function used in cluster_snarls, and copied from SnarlNormalizer.
    // SubHandleGraph extract_subgraph(const HandleGraph &graph,
    //                                              const id_t leftmost_id,
    //                                              const id_t rightmost_id);

    // // Another function used in cluster_snarls:
    // bool snarls_adjacent(const pair<id_t, id_t>& snarl_1, const pair<id_t, id_t>& snarl_2);

    // bool is_trivial(const pair<id_t, id_t>& snarl);
    
    ///Functions used in desegregate_nodes:
    handle_t overwrite_node_id(const id_t old_node_id, const id_t new_node_id);
    
    ///Debug function:
    // std::vector<vg::RebuildJob::mapping_type> debug_replace_node_six()

    
};
}
}

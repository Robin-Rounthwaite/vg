#include "../gbwt_helper.hpp"
// #include "../handle.hpp"
#include <gbwtgraph/gbwtgraph.h>
// #include "../handle.hpp"
// #include "../vg.hpp"
// #include "../types.hpp"


////////////////////////////////////////////////////////////
/// Use the gbwt_changelog to perform paralellized rebuild_gbwt():
////////////////////////////////////////////////////////////
namespace vg {
namespace algorithms {

gbwt::GBWT apply_gbwt_changelog(const gbwtgraph::GBWTGraph &gbwt_graph, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog, const gbwt::GBWT& gbwt, const int& threads, const bool& debug_print);

std::unordered_map<nid_t, size_t> get_node_to_job(const vector<unordered_set<nid_t>>& weakly_connected_components);

std::vector<RebuildJob> divide_changelog_into_jobs(const std::unordered_map<nid_t, size_t>& node_to_job, const vector<unordered_set<nid_t>>& weakly_connected_components, const std::vector<vg::RebuildJob::mapping_type>& gbwt_changelog);

RebuildParameters set_update_gbwt_parameters(const int& threads, const bool& debug_print);

}
}
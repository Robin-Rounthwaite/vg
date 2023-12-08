#ifndef VG_MSA_CONVERTER_HPP_INCLUDED
#define VG_MSA_CONVERTER_HPP_INCLUDED

/** \file
 * msa_converter.hpp: contains a class that can construct VGs from clustal MSAs
 */

#include <vector>
#include <set>
#include <map>
#include <cstdlib>
#include <functional>
#include <regex>

#include <vg/vg.pb.h>

#include "vg.hpp"

namespace vg {

using namespace std;

    class MSAConverter : public Progressive  {
    public:
        
        MSAConverter();
        ~MSAConverter();
        
        void load_alignments(istream& in, string format = "fasta");

        void load_alignments_from_vector(vector<string> input);
        
        unique_ptr<MutablePathDeletableHandleGraph> make_graph(bool keep_paths = true, size_t max_node_length = numeric_limits<size_t>::max());


        class MSANode {
            public:
                virtual ~MSANode() = default;
                MSANode(int id, string seq);
                MSANode(const MSANode&) = delete; //this makes sure that whenever I copy, compiler will protest that.
                MSANode& operator=(const MSANode&) = delete;

                // note: use right_neighbors.emplace(), which lets me put things in the set even if I cant move or copy it.
                set<MSANode*> right_neighbors; //todo-done: bad because msanodes will move but the pointer wont track new location. 
                set<MSANode*> left_neighbors;


                int get_id() const;

                string get_sequence();
                string* get_mutable_sequence();
                
            // private:
                int _id = -1;
                string _seq = "";

                bool operator==(const MSANode& otherNode) const
                {
                    if (this->_id == otherNode.get_id()) return true;
                    else return false;
                }

                struct HashFunction
                {
                    size_t operator()(const MSANode& node) const
                    {
                        return std::hash<int>()(node.get_id());
                        // return node.get_id();
                    // size_t xHash = std::hash<int>()(node.x);
                    // size_t yHash = std::hash<int>()(node.y) << 1;
                    // return xHash ^ yHash;
                    }
                };
        };

        class MSAGraph {
        public:
            virtual ~MSAGraph() = default;
            MSAGraph();

            MSANode* create_node(string str);
            // MSANode* create_node(string str);
            // MSANode& create_node(string str);
            // MSANode* create_node(int id, string str);
            void create_edge(MSANode* left, MSANode* right);

            void destroy_node(MSANode& node);

            unique_ptr<MutablePathDeletableHandleGraph> make_handlegraph();
            // int get_node_count();
            // int get_edge_count();

        // protected:

            int _highest_node_id = 0;


            unordered_map<int, shared_ptr<MSANode>> _ids_to_nodes; //todo-done: make obj guarantee that you can grow wihtout deallocating (e.g. deck) (or unique pointers put into the unordered_map. otherwise deconstruct wouldn't work; shared pointer would work too.)
            // set<pair<MSANode*, MSANode*>> _edge_set;

        };



    private:
        
        vector<unordered_map<string, string>> alignments;
        



    };

}

#endif

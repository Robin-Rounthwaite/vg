//
//  suffix_tree.cpp
//  
//
//  Created by Jordan Eizenga on 1/25/17.
//
//

#include <stdio.h>
#include <iostream>
#include <unordered_map>

#include "msa_converter.hpp"
#include "catch.hpp"

using namespace std;

namespace vg {
    namespace unittest {
        
        TEST_CASE("MSAConverter can build from each kind of input", "[msa]") {
            
            
            SECTION("MSAConverter can build from FASTA input") {
                
                string input = ">seq1 description 1\nAAGTGATAGAGATAAA\nGTGATAGAGATA\n>seq2 description2\n-AAGTATAGACATAAA\nGTTAGAGATA--\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
            }
            
            SECTION("MSAConverter can build from Clustal input") {
                
                string input = "CLUSTAL O(1.2.4) multiple sequence alignment\n\nGI262359905      TGTACCTTGATTTCGTATTCTGAGAGGCTGCTGCTTAGCGGTAGCCCCT-TGGTTTCCGT\nGI528476558      ---------------------------GTGGAAGTGTTTGCTACCAAGTTTATTTGCAGT\nref              ---------------------------GTGGAAGTGTTTGCTACCAAGTTTATTTGCAGT\n                **    *    * ** *   * *  ** * **\n\nGI262359905      GGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCT\nGI528476558      GTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCT----\nref              GTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCT----\n                *  *** * * * *      *  **  **    ** **    * ** * * **                \n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "clustal");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
            }
        }
        
        TEST_CASE("MSAConverter produces correct graphs from MSAs", "[msa]") {
            
            SECTION("MSAConverter produces one node from a completely matching alignment") {
                string input = ">seq1\nAAA\n>seq2\nAAA\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
                
                // REQUIRE(graph->get_node_count() == 1);
                // REQUIRE(graph->get_edge_count() == 0);
                // REQUIRE(g.node(0).sequence() == "AAA");

                REQUIRE(graph->get_node_count() == 1);
                REQUIRE(graph->get_edge_count() == 0);
                REQUIRE(graph->get_sequence(graph->get_handle(1)) == "AAA");
            }
            
            SECTION("MSAConverter respects the max node length") {
                
                string input = ">seq1\nAAA\n>seq2\nAAA\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                // unique_ptr<HandleGraph> graph = msa_converter.make_graph(true, 1);
                unique_ptr<HandleGraph> graph = msa_converter.make_graph(1);
                
                REQUIRE(graph->get_node_count() == 3);
                REQUIRE(graph->get_edge_count() == 2);
            }
            
            SECTION("MSAConverter splits columns that mismatch into multiple nodes") {
                
                string input = ">seq1\nATG\n>seq2\nACG\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
                
                unordered_map<string, handle_t> handles;
                
                graph->for_each_handle([&](handle_t handle)
                {
                    handles[graph->get_sequence(handle)] = handle;
                });
                
                REQUIRE(graph->get_node_count() == 4);
                REQUIRE(graph->get_edge_count() == 4);
                
                REQUIRE(handles.count("A"));
                REQUIRE(handles.count("C"));
                REQUIRE(handles.count("G"));
                REQUIRE(handles.count("T"));
                
                REQUIRE(graph->has_edge(handles["A"], handles["C"]));
                REQUIRE(graph->has_edge(handles["A"], handles["T"]));
                REQUIRE(graph->has_edge(handles["T"], handles["G"]));
                REQUIRE(graph->has_edge(handles["C"], handles["G"]));
            }
            
            SECTION("MSAConverter adds edges over gap") {
                
                string input = ">seq1\nA-G\n>seq2\nACG\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
                
                unordered_map<string, handle_t> handles;
                graph->for_each_handle([&](handle_t handle)
                {
                    handles[graph->get_sequence(handle)] = handle;
                });
                

                REQUIRE(graph->get_node_count() == 3);
                REQUIRE(graph->get_edge_count() == 3);
                
                REQUIRE(handles.count("A"));
                REQUIRE(handles.count("C"));
                REQUIRE(handles.count("G"));
                
                REQUIRE(graph->has_edge(handles["A"], handles["C"]));
                REQUIRE(graph->has_edge(handles["A"], handles["G"]));
                REQUIRE(graph->has_edge(handles["C"], handles["G"]));
            }
            
            SECTION("MSAConverter handles overlapping gaps") {
                
                string input = ">seq1\nAA--GTT\n>seq2\nAAACGTT\n>seq3\nAAA--TT\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
                
                unordered_map<string, handle_t> handles;
                graph->for_each_handle([&](handle_t handle)
                {
                    handles[graph->get_sequence(handle)] = handle;
                });
                
                REQUIRE(graph->get_node_count() == 5);
                REQUIRE(graph->get_edge_count() == 6);
                
                REQUIRE(handles.count("AA"));
                REQUIRE(handles.count("A"));
                REQUIRE(handles.count("C"));
                REQUIRE(handles.count("G"));
                REQUIRE(handles.count("TT"));
                
                REQUIRE(graph->has_edge(handles["AA"], handles["A"]));
                REQUIRE(graph->has_edge(handles["AA"], handles["G"]));
                REQUIRE(graph->has_edge(handles["A"], handles["C"]));
                REQUIRE(graph->has_edge(handles["A"], handles["TT"]));
                REQUIRE(graph->has_edge(handles["C"], handles["G"]));
                REQUIRE(graph->has_edge(handles["G"], handles["TT"]));
            }
            
            SECTION("MSAConverter handles nested gaps") {
                
                string input = ">seq1\nAAACGTT\n>seq2\nAA---TT\n>seq3\nAAA-GTT\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                unique_ptr<HandleGraph> graph = msa_converter.make_graph();
                
                unordered_map<string, handle_t> handles;
                graph->for_each_handle([&](handle_t handle)
                {
                    handles[graph->get_sequence(handle)] = handle;
                });
                
                REQUIRE(graph->get_node_count() == 5);
                REQUIRE(graph->get_edge_count() == 6);
                
                REQUIRE(handles.count("AA"));
                REQUIRE(handles.count("A"));
                REQUIRE(handles.count("C"));
                REQUIRE(handles.count("G"));
                REQUIRE(handles.count("TT"));
                
                REQUIRE(graph->has_edge(handles["AA"], handles["A"]));
                REQUIRE(graph->has_edge(handles["AA"], handles["TT"]));
                REQUIRE(graph->has_edge(handles["A"], handles["C"]));
                REQUIRE(graph->has_edge(handles["A"], handles["G"]));
                REQUIRE(graph->has_edge(handles["C"], handles["G"]));
                REQUIRE(graph->has_edge(handles["G"], handles["TT"]));
            }
        }
    }
}

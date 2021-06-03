// #include "../../../git_repositories/dozeu/dozeu.h"
// #include "0_oo_normalize_snarls.hpp"

// #include <seqan/align.h>
// #include <seqan/graph_align.h>
// #include <seqan/graph_msa.h>

// #include "../msa_converter.hpp"

// namespace vg {
// namespace algorithms{

// VG SnarlNormalizer::align_haplotypes(const unordered_set<string>& source_to_sink_haplotypes, const unordered_set<string>& source_only_haplotypes, const unordered_set<string>& sink_only_haplotypes, const unordered_set<string>& other_haplotypes) 
// {
//     //todo: add the proper haplotypes to the proper sections of align_haplotypes. Currently passing all source_to_sink_haps to each argument.
//     //first, align the source-to-sink haplotypes using seqan to create a graph.
//     VG new_snarl = align_source_to_sink_haplotypes(source_to_sink_haplotypes);


//     // // if there are haplotypes besides source-to-sink haps, convert new_snarl to a dozeu-compatible format
//     // // and align all other haplotypes to the dozeu graph, updating the graph as we go.
//     // if (source_only_haplotypes.size() != 0 || sink_only_haplotypes.size() != 0 || other_haplotypes.size() != 0) 
//     // {
//     //     dozeu_align_to_graph()
//     // }



//     // convert the dozeu graph back to a handlegraph.

//     return new_snarl;
// }
//         // dozeu_align_to_graph()


// // TODO: eventually change to deal with haplotypes that start/end in middle of snarl.
// // Aligns haplotypes to create a new _graph using MSAConverter's seqan converter.
// //      Assumes that each haplotype stretches from source to sink.
// // Arguments:
// //      source_to_sink_haplotypes: a vector of haplotypes in string format (concat of
// //      handle sequences).
// // Returns:
// //      VG object representing the newly realigned snarl.
// VG SnarlNormalizer::align_source_to_sink_haplotypes(
//     const unordered_set<string>& source_to_sink_haplotypes) 
// {
//     // cerr << "align_source_to_sink_haplotypes" << endl;
//     // cerr << " haplotypes in source_to_sink_haplotypes: " << endl;
//     // for (string hap : source_to_sink_haplotypes) {
//     //     cerr << hap << endl;
//     // }
//     // cerr << "number of strings to align: " << source_to_sink_haplotypes.size() << endl;
//     // TODO: make the following comment true, so that I can normalize haplotypes that
//     // TODO:    aren't source_to_sink by adding a similar special character to strings in
//     // TODO:    the middle of the snarl.
//     // modify source_to_sink_haplotypes to replace the leading and
//     // trailing character with a special character. This ensures that the leading char of
//     // the haplotype becomes the first character in the newly aligned snarl's source - it
//     // maintains the context of the snarl.

//     // store the source/sink chars for later reattachment to source and sink.
//     string random_element;
//     for (auto hap : source_to_sink_haplotypes){
//         random_element = hap;
//         break;
//     }
//     string source_char(1, random_element.front());
//     string sink_char(1, random_element.back());

//     // cerr << "strings in path_seq before replacing final character: " << endl;
//     // for (auto path : get<0>(haplotypes))
//     // {
//     //     cerr << path << endl;f
//     // }

//     // replace the source and sink chars with X, to force match at source and sink.
//     unordered_set<string> edited_source_to_sink_haplotypes;
//     // for (auto it = source_to_sink_haplotypes.begin(); it != source_to_sink_haplotypes.end(); it++)
//     for (auto hap : source_to_sink_haplotypes)
//     {
//         // cerr << "hap before replace: " << hap << endl;
//         hap.replace(0, 1, "X");
//         hap.replace(hap.size() - 1, 1, "X");
//         // cerr << "hap after replace: " << hap << endl;
//         edited_source_to_sink_haplotypes.emplace(hap);
//     }
//     // cerr << "source_char: " << source_char << endl;
//     // cerr << "sink_char: " << sink_char << endl;

//     // //todo: debug_statement
//     // source_to_sink_haplotypes.emplace_back("XX");

//     // /// make a new scoring matrix with _match=5, _mismatch = -3, _gap_extend = -1, and
//     // _gap_open = -3, EXCEPT that Q has to be matched with Q (so match score between Q
//     // and Q =len(seq)+1)
//     // // 1. Define type and constants.
//     // //
//     // // Define types for the score value and the scoring scheme.
//     // typedef int TValue;
//     // typedef seqan::Score<TValue, seqan::ScoreMatrix<seqan::Dna5, seqan::Default> >
//     // TScoringScheme;
//     // // Define our gap scores in some constants.
//     // int const gapOpenScore = -1;
//     // int const gapExtendScore = -1;

//     // static int const _data[TAB_SIZE] =
//     //     {
//     //         1, 0, 0, 0, 0,
//     //         0, 1, 0, 0, 0,
//     //         0, 0, 1, 0, 0,
//     //         0, 0, 0, 1, 0,
//     //         0, 0, 0, 0, 0
//     //     };

//     // create seqan multiple_sequence_alignment object
//     //// seqan::Align<seqan::DnaString>   align;
//     seqan::Align<seqan::CharString> align;

//     seqan::resize(rows(align), edited_source_to_sink_haplotypes.size());
//     int i = 0;
//     for (auto hap : edited_source_to_sink_haplotypes) {
//         assignSource(row(align, i), hap.c_str());
//         i++;
//     }

//     globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));

//     vector<string> row_strings;
//     for (auto &row : rows(align)) {
//         string row_string;
//         auto it = begin(row);
//         auto itEnd = end(row);
//         for (; it != itEnd; it++) {
//             row_string += *it;
//         }
//         // todo: debug_statement
//         // cerr << "ROW_STRING: " << row_string << endl;
//         // edit the row so that the proper source and sink chars are added to the
//         // haplotype instead of the special characters added to ensure correct alignment
//         // of source and sink.
//         // cerr << "row_string before: " << row_string << endl;
//         row_string.replace(0, 1, source_char);
//         row_string.replace(row_string.size() - 1, 1, sink_char);
//         row_strings.push_back(row_string);
//         // cerr << "row_string after: " << row_string << endl;
//     }

//     stringstream ss;
//     for (string seq : row_strings) {
//         // todo: debug_statement
//         // cerr << "seq in alignment:" << seq << endl;
//         ss << endl << seq;
//     }
//     // ss << align;
//     MSAConverter myMSAConverter = MSAConverter();
//     myMSAConverter.load_alignments(ss, "seqan");
//     VG snarl = myMSAConverter.make_graph();
//     snarl.clear_paths();

//     pair<vector<handle_t>, vector<handle_t>> source_and_sink =
//         debug_get_sources_and_sinks(snarl);

//     // TODO: throw exception(?) instead of cerr, or remove these messages if I'm confident
//     // TODO:    code works.
//     if (source_and_sink.first.size() != 1) {
//         cerr << "WARNING! Snarl realignment has generated "
//              << source_and_sink.first.size() << " source nodes." << endl;
//     }

//     if (source_and_sink.second.size() != 1) {
//         cerr << "WARNING! Snarl realignment has generated "
//              << source_and_sink.second.size() << " sink nodes." << endl;
//     }

//     return snarl;
// }


// }
// }
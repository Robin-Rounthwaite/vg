#include "0_oo_normalize_snarls.hpp"
#include <spoa/spoa/spoa.hpp>
#include <abpoa/abpoa.h>
#include "../msa_converter.hpp"
// extern "C" {
// #include <libkalign.h>
// }
namespace vg {
namespace algorithms{


//Note: you probably don't want to set output_msa to true while sending other things to cout. Because the msa is outputted to cout as well.
VG SnarlNormalizer::poa_source_to_sink_haplotypes(const unordered_set<string>& source_to_sink_haplotypes, const int snarl_num, const bool output_msa) {
    // uint8_t ***msa_seq = NULL;

    
    // replace the source and sink chars with X, to force match at source and sink.
    //todo: if copying out the set takes significant time, try using unordered_set::extract() here. (https://stackoverflow.com/questions/42519867/efficiently-moving-contents-of-stdunordered-set-to-stdvector)
    vector<string> edited_source_to_sink_haplotypes(source_to_sink_haplotypes.begin(), source_to_sink_haplotypes.end());
    sort(edited_source_to_sink_haplotypes.begin(), edited_source_to_sink_haplotypes.end(), []
        (const std::string& first, const std::string& second){
            return first.size() > second.size();
        });
        
    char first_char = edited_source_to_sink_haplotypes.back().front();
    char last_char = edited_source_to_sink_haplotypes.back().back();
    bool two_char_hap_exists = 0;
    // cerr << "haps, in sorted order, to be aligned: " << endl;
    for (string hap : edited_source_to_sink_haplotypes)
    {
        // cerr << hap << endl;
        if (hap.size() == 2) 
        {
            // then we need to remove this sequence from the alignment process. Because 
            // part of the preprocessing involves trimming off the first and last 
            // characters of the string. And if we do that to this string, there will be 
            // nothing left!
            if (source_to_sink_haplotypes.size() == 2)
            {
                // then we're eliminating one of the two sequences to be aligned, and so 
                // there's no point in doing the alignment at all. Just manually construct
                // the alignment instead.
                vector<string> msa_output;
                for (string hap_to_align: edited_source_to_sink_haplotypes)
                {
                    if (hap_to_align.size() > 2)
                    {
                        // char dash_char = "-";
                        string two_char_hap((hap_to_align.size() - 2), '-');
                        // cerr << "two_char_hap pre-character additions: " << two_char_hap << endl;
                        two_char_hap.insert(0, 1, first_char);
                        // cerr << "first char and last char: " << first_char << " " << last_char << endl;
                        two_char_hap.insert(two_char_hap.size(), 1, last_char);
                        // cerr << "two_char_hap generated for manual MSA: " << two_char_hap << endl;

                        msa_output.push_back(two_char_hap);
                        msa_output.push_back(hap_to_align);
                    }
                }
                

                MSAConverter myMSAConverter = MSAConverter();
                myMSAConverter.load_alignments_from_vector(msa_output);
                VG snarl = myMSAConverter.make_graph();
                
                snarl.clear_paths();

                return snarl;
            }
            else
            {
                // remove the sequence from the alignment, but remember to add it back in 
                // before making the new graph. 
                two_char_hap_exists = 1;
                // it's the smallest hap, so it'll be at the back.
                edited_source_to_sink_haplotypes.pop_back(); 
            }
        }
    }
    // cerr << "deleting first and last characters:" << endl;
    for (auto it = edited_source_to_sink_haplotypes.begin(); it != edited_source_to_sink_haplotypes.end(); it++)
    {
        // cerr << "before delete: " << *it << endl;
        it->erase(0, 1);
        it->erase(it->size() - 1, 1);
        // cerr << "after delete: " << *it << endl;
    }

    //todo: make sPOA and abPOA sections separate functions?
    // if the alignment is for <=750 bases, use sPOA. 
    vector<string> msa_output;
    if (edited_source_to_sink_haplotypes.front().size() <= 750)
    { 
        // auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -3, -3, -1); // original inputs. m n g e, probably means "match, mismatch, gap, extend"?
        // auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -5, -5, -1); // modified inputs, with increased gap-open cost.
        auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -10, -10, -1); // modified inputs, with increased gap-open cost.
        spoa::Graph spoa_graph{};

        for (string sequence : edited_source_to_sink_haplotypes)
        {
            auto alignment = alignment_engine->Align(sequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, sequence);
        }
        //todo: figure out good threshold for when to add a consensus sequence. Greater than 3 or more? Or just 2?
        if (edited_source_to_sink_haplotypes.size() > 2)
        {
            // if there are more than two(?) sequences being aligned, use the consensus sequence as a scaffold.
            // use initial POA to get a consensus sequence
            auto consensus = spoa_graph.GenerateConsensus();

            // Seed the alignment with the consensus
            spoa::Graph seeded_spoa_graph{};
            auto seeded_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 6, -2, -4, -1);
            auto cons_alignment = seeded_alignment_engine->Align(consensus, seeded_spoa_graph);
            seeded_spoa_graph.AddAlignment(cons_alignment, consensus);

            // Iterate a second time on alignment, this time with consensus as the seed
            for (string sequence : edited_source_to_sink_haplotypes)
            {
                auto alignment = alignment_engine->Align(sequence, seeded_spoa_graph);
                seeded_spoa_graph.AddAlignment(alignment, sequence);
            }
            msa_output = seeded_spoa_graph.GenerateMultipleSequenceAlignment();
        }
        else
        {
            //otherwise, just directly generate the MSA output.
            msa_output = spoa_graph.GenerateMultipleSequenceAlignment();
        }
    }
                                            //for alignments of sequences > 750 bases, use abPOA
                                            // else
                                            // {
                                            //     // this is annoyingly not accessible by any header, so i'll just copy it
                                            //     static const uint8_t ab_nt4_table[256] = {
                                            //         0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                                            //         4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
                                            //     };

                                            //     abpoa_t* abpoa = abpoa_init();
                                                
                                            //     vector<uint8_t*> encoded_seqs;
                                            //     vector<char*> seq_names;
                                            //     vector<int> seq_lens;
                                            //     int seq_count = 0;
                                            //     for (string& sequence : edited_source_to_sink_haplotypes)
                                            //     {
                                            //         // convert sequence to 0-4 vals and C array
                                            //         uint8_t* enc_seq = (uint8_t*) malloc(sequence.size() * sizeof(uint8_t));
                                            //         for (size_t j = 0; j < sequence.size(); ++j) {
                                            //             enc_seq[j] = ab_nt4_table[(uint8_t)sequence[j]];
                                            //         }
                                            //         // make a C-style string for the name
                                            //         vector<char> name_vec;
                                            //         string name_str = to_string(seq_count);
                                            //         for (char name_char : name_str)
                                            //         {
                                            //             name_vec.push_back(name_char);
                                            //         }
                                            //         name_vec.push_back('\0');
                                                    
                                            //         // remember the results
                                            //         encoded_seqs.push_back(enc_seq);
                                            //         seq_names.push_back(&name_vec[0]);
                                            //         seq_lens.push_back(sequence.size());
                                            //     }
                                                
                                            //     //todo: figure out the abpoa params. Here are where the params are set in Bluntifier.cpp:
                                            //     // set up abPOA
                                            //     abpoa_para_t* abpoa_params = abpoa_init_para();
                                                
                                            //     abpoa_params->align_mode = 0; // global alignment
                                            //     abpoa_params->match = 6;
                                            //     abpoa_params->mismatch = 2;
                                            //     abpoa_params->gap_mode = ABPOA_AFFINE_GAP;
                                            //     abpoa_params->gap_open1 = 4;
                                            //     abpoa_params->gap_ext1 = 2;
                                            //     abpoa_params->wb = 10; // band width for adaptive banding
                                            //     abpoa_params->end_bonus = 0;
                                                
                                            //     // we don't care about the MSA
                                            //     // TODO: maybe we should do a consensus realignment like with SPOA though?
                                            //     abpoa_params->out_msa = 1;
                                            //     // abpoa_params->out_msa_header = ??
                                            //     abpoa_params->out_gfa = 0;
                                            //     abpoa_params->progressive_poa = 1;
                                            //     abpoa_params->use_read_ids = 1;
                                                
                                            //     // // even though we only care about the POA graph, we have to tell abpoa
                                            //     // // that we want some output, or else it won't do the alignment
                                            //     // abpoa_params->out_cons = 1;
                                                
                                            //     // finalize the params
                                            //     abpoa_post_set_para(abpoa_params);

                                            //     FILE* abpoa_output = tmpfile();
                                            //     uint8_t ***cons_seq;
                                            //     int ***cons_cov;
                                            //     int **cons_l;
                                            //     int *cons_n;
                                            //     uint8_t ***msa_seq;
                                            //     int *msa_l;
                                            //     *msa_l = 1;
                                                
                                            //     //todo: anything else I need to add before calling abpoa_msa?
                                            //     // int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seqs, char **seq_names, int *seq_lens, uint8_t **seqs, FILE *out_fp, uint8_t ***cons_seq, int ***cons_cov, int **cons_l, int *cons_n, 
                                            //     //uint8_t ***msa_seq, int *msa_l);
                                            
                                            //     // and now actually perform the MSA now that we've added everything
                                            //     abpoa_msa(abpoa,
                                            //             abpoa_params,
                                            //             encoded_seqs.size(),
                                            //             seq_names.data(),
                                            //             seq_lens.data(),
                                            //             encoded_seqs.data(),
                                            //             abpoa_output,
                                            //             cons_seq,
                                            //             cons_cov,
                                            //             cons_l,
                                            //             cons_n,
                                            //             msa_seq,
                                            //             msa_l);

                                            //     // abpoa_output(abpoa,
                                            //     //         abpoa_params,
                                            //     //         encoded_seqs.size(),
                                            //     //         seq_names.data(),
                                            //     //         seq_lens.data(),
                                            //     //         encoded_seqs.data(),
                                            //     //         abpoa_output,
                                            //     //         cons_seq,
                                            //     //         cons_cov,
                                            //     //         cons_l,
                                            //     //         cons_n,
                                            //     //         msa_seq,
                                            //     //         msa_l);

                                            //     // abpoa_generate_rc_msa(abpoa, abpoa_params, abpoa_output, msa_seq, msa_l);

                                            //     cerr << "here (hopefully) is the abpoa output: " << endl;
                                            //     char buffer [100];
                                            //     while ( ! feof (abpoa_output) )
                                            //     {
                                            //     if ( fgets (buffer , 100 , abpoa_output) == NULL ) break;
                                            //     fputs (buffer , stderr);
                                            //     }
                                            //     fclose (abpoa_output);
                                                
                                            //     // annoyingly not available in header, so we copy it here
                                            //     static const char ab_nt256_table[256] = {
                                            //         'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
                                            //         'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
                                            //     };
                                                






                                            //     // std::string line;
                                            //     // while (getline(abpoa_output, line))
                                            //     // {
                                            //     //     std::istringstream iss(line);
                                            //     //     int a, b;
                                            //     //     if (!(iss >> a >> b)) { break; } // error

                                            //     //     // process pair (a,b)
                                            //     // }

                                            //     cerr << "here (hopefully) is the contents of msa_seq: " << endl; 
                                            //     cerr << "msa_seq 0 0 0" << *msa_seq << endl;
                                            //     cerr << "msa_seq 0 0 0" << **msa_seq << endl;
                                            //     cerr << "msa_seq 0 0 0" << ***msa_seq << endl;
                                            //     // for (uint8_t*** i = msa_seq; i != encoded_seqs.size(); i ++)
                                            //     for(size_t i=0;i<encoded_seqs.size();++i)
                                            //     {
                                            //         for (int j = 0; j != *msa_l; j ++)
                                            //         {
                                            //             cerr << *msa_seq[i][j]; 
                                            //         }
                                            //         cerr << endl;
                                            //     }

                                            //     abpoa_free_para(abpoa_params);
                                                
                                            //     exit(1);
                                            // }

                                            // seeded_spoa_graph.ExtractSubgraph

                                            // cerr << "aligned haps in sorted order before X's removed: " << endl;
                                            // for (auto hap : msa_output)
                                            // {
                                            //     cerr << hap << endl;
                                            // }

    // cerr << "here is the msa_output from seeded_spoa_graph." << endl;
    // for (string& output_line : msa_output)
    for (auto it = msa_output.begin(); it != msa_output.end(); it++)
    {
        // it->replace(0, 1, string(1, first_char));
        // it->replace(it->size() - 1, 1, string(1, last_char)); 
        // cerr << "before insert: " << *it << endl;
        it->insert(0, 1, first_char);
        it->insert(it->size(), 1, last_char);
        // cerr << "after insert: " << *it << endl;
    }

    // cerr << "aligned haps to be given to MSAConverter: " << endl;
    // for (auto hap : msa_output)
    // {
    //     cerr << hap << endl;
    // }

    if (two_char_hap_exists)
    {
        string two_char_hap((msa_output.back().size() - 2), '-');
        two_char_hap.insert(0, 1, first_char);
        two_char_hap.insert(two_char_hap.size(), 1, last_char);
        // cerr << "two_char_hap generated after the MSA: " << two_char_hap << endl;
        msa_output.push_back(two_char_hap);
    }

    if (output_msa)
    {
        //then print the msa straight to cout.
        for (auto hap : msa_output)
        {
            // cout << hap << endl;
            cerr << hap << endl;
        }
    }

    MSAConverter myMSAConverter = MSAConverter();
    myMSAConverter.load_alignments_from_vector(msa_output);
    VG snarl = myMSAConverter.make_graph();
    
    snarl.clear_paths();
    // cerr << "Hi! Testing." << endl;

    return snarl;

}

// VG SnarlNormalizer::kalign_source_to_sink_haplotypes(const unordered_set<string>& source_to_sink_haplotypes) {
//     //c_str_array will contain reformatted source_to_sink_haplotypes as a char* array:
//     vector<char *> c_str_array;
//     c_str_array.reserve(source_to_sink_haplotypes.size());
//     //c_str_sizes will contain the size of each string in c_str_array.
//     vector<int> c_str_sizes;
    
//     for (string hap : source_to_sink_haplotypes)
//     {
//         vector<char> str_copy;
//         str_copy.reserve(hap.size());
//         for (char base : hap)
//         {
//             str_copy.push_back(base);
//         }
//         c_str_sizes.push_back(hap.size());
//         c_str_array.push_back(&str_copy[0]);
//     }

//     //todo: will the objects placed in this vector be destructed when the vector is destructed? Or do I need to do that manually?
//     vector<char **> align_output;
//     int aln_len;
//     kalign(&c_str_array[0], &c_str_sizes[0], source_to_sink_haplotypes.size(), &align_output[0], &aln_len);

//     vector<string> align_output_strings;
//     for (auto align_c_str_array : align_output)
//     {
//         string align_output_string;
//         align_output_string.reserve(aln_len);
//         for (int i = 0; i != aln_len; i++)
//         {
//             align_output_string += align_c_str_array[i];
//         }
//         cerr << "align_output_string is: " << align_output_string << endl;
        
//     }

//     //todo: make this the proper VG output from myMSAConverter.
//     VG new_vg;
//     return new_vg;
// }

}
}

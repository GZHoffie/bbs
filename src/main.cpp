#include <sharg/all.hpp>
#include "./bbs.hpp"
#include "./utils.hpp"


std::vector<seqan3::dna4_vector> _string_vector_to_dna4(const std::vector<std::string>& sequences) {
    std::vector<seqan3::dna4_vector> res;
    for (auto& s : sequences) {
        seqan3::dna4_vector dna_sequence;
        for (auto &nt : s) {
            seqan3::dna4 dna = seqan3::assign_char_to(nt, seqan3::dna4{});
            dna_sequence.push_back(dna);
        }
        res.push_back(dna_sequence);
        //seqan3::debug_stream << res.size() << " " << sequences.size() << "\n";
    }
    return res;
}

std::string _dna4_to_string(const seqan3::dna4_vector& sequence) {
    std::string res;
    for (auto & nt : sequence) {
        res += nt.to_char();
    }
    return res;
}

 
void run_program(std::filesystem::path const & input_path, const std::string& separator, greedy_ensembler<unsigned int>* g_ensembler, 
                 unsigned int read_length, unsigned int beam_width, bool single_sided, std::filesystem::path const & output_path, uint8_t fixed_k)
{
    Timer timer;

    timer.tick();

    // Create the ensembler
    std::ifstream clusters_stream(input_path);

    std::vector<std::string> cluster;
    unsigned int num_clusters = 0;
    unsigned int num_sequences = 0;

    std::ofstream output_stream;

    if (!output_path.empty()) {
        output_stream.open(output_path);
        output_stream << "read_id,reconstruction_result,k,path_weight,confidence\n";
    }

    while (clusters_stream) {
        //seqan3::debug_stream << center;
        cluster.clear();
        
        // read the cluster of sequences
        std::string temp;
        while (clusters_stream >> temp) {
            // if the separator is found as substring of temp, break
            if (temp.find(separator) != std::string::npos) break;
            cluster.push_back(temp);
        }

        if (cluster.empty() && num_clusters == 0) {
            continue;
        }

        // turn the cluster to seqan3::dna4_vector's.
        auto dna4_cluster = _string_vector_to_dna4(cluster);

        // ensembling
        if (dna4_cluster.size() >= 0) {
            auto [res, k, weight, confidence] = g_ensembler->ensemble(dna4_cluster, read_length, beam_width, !single_sided, fixed_k);
            
            // output the results to the output file
            for (auto nt : res) {
                std::cout << nt.to_char();
            }
            num_clusters++;
            num_sequences += dna4_cluster.size();
            if (!output_path.empty()) {
                output_stream << num_clusters << "," << _dna4_to_string(res) << "," << k << "," << weight << "," << confidence << "\n";
            }
            if (clusters_stream) {
                std::cout << "\n";
            }
        }
    }
    timer.tock();

    seqan3::debug_stream << "Time used: " << timer.elapsed_seconds() << " s.\n";
    seqan3::debug_stream << "Total number of clusters: " << num_clusters << ".\n";
    seqan3::debug_stream << "Average number of sequences per cluster: " << (float)num_sequences / num_clusters << ".\n";
    seqan3::debug_stream << "Average time per cluster: " << (float)timer.elapsed_seconds() / num_clusters << " s.\n";
    if (!output_path.empty()) {
        output_stream.close();
        seqan3::debug_stream << "Output written to " << output_path << ".\n";
    }
}
 
struct cmd_arguments
{
    unsigned int read_length{100};
    unsigned int beam_width{20};
    std::string input_path{};
    std::string separator{"CLUSTER"};
    unsigned int k_lower_bound{5};
    unsigned int k_upper_bound{14};
    float alpha{1};
    bool single_sided{false};
    std::filesystem::path output_path{};
    uint8_t fixed_k{0};
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Gu Zhenhao";
    parser.info.short_description = "Bidirectional beam search for efficient trace reconstruction.";
    parser.info.version = "0.1.0";
    parser.add_option(args.input_path,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input",
                                    .description = "The path to the input cluster file.",
                                    .required = true});
    
    parser.add_option(args.read_length,
                      sharg::config{.short_id = 'l',
                                    .long_id = "read_length",
                                    .description = "The length of the reads.",
                                    .required = true});
    
    parser.add_option(args.separator,
                      sharg::config{.short_id = 's',
                                    .long_id = "separator",
                                    .description = "String that indicate a new cluster.",
                                    .required = true});
    
    parser.add_option(args.beam_width,
                      sharg::config{.short_id = 'b',
                                    .long_id = "beam_width",
                                    .description = "The width of the beam search."});

    parser.add_option(args.k_lower_bound,
                        sharg::config{.short_id = 'k',
                                        .long_id = "k_lower_bound",
                                        .description = "The lower bound of the k-mer length."});

    parser.add_option(args.k_upper_bound,
                        sharg::config{.short_id = 'K',
                                        .long_id = "k_upper_bound",
                                        .description = "The upper bound of the k-mer length."});
    
    parser.add_option(args.fixed_k,
                        sharg::config{.short_id = 'f',
                                        .long_id = "fixed_k",
                                        .description = "If this is set to a positive value (>= 1), the BBS algorithm ignores values of -k and -K option and runs the reconstruction program with the fixed k value."});

    parser.add_option(args.alpha,
                        sharg::config{.short_id = 'a',
                                        .long_id = "alpha",
                                        .description = "The alpha value for Laplace smoothing in estimation of conditional probabilities."});

    parser.add_option(args.single_sided,
                        sharg::config{.short_id = 'S',
                                        .long_id = "single_sided",
                                        .description = "Whether to use single sided beam search."});
    
    parser.add_option(args.output_path,
                        sharg::config{.short_id = 'o',
                                        .long_id = "output",
                                        .description = "The path to the detailed output file.",
                                        .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"csv"}}});
}
 
int main(int argc, char const ** argv)
{
    sharg::parser parser("BBS", argc, argv, sharg::update_notifications::off);
    cmd_arguments args{};
 
    initialise_parser(parser, args);
 
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    // Check the arguments
    if (args.k_lower_bound > args.k_upper_bound) {
        std::cerr << "The lower bound of k-mer length should be smaller than the upper bound.\n";
        return -1;
    }
    if (args.k_upper_bound >= 15 || args.fixed_k >= 15) {
        std::cerr << "The current version of BBS implementation only allows k-mer length smaller than 15.\n";
        return -1;
    }

    // Create the ensembler
    greedy_ensembler<unsigned int>* g_ensembler;
    g_ensembler = new greedy_ensembler<unsigned int>(args.k_lower_bound, args.k_upper_bound, args.alpha);


    run_program(args.input_path, args.separator, g_ensembler, args.read_length, args.beam_width, args.single_sided, args.output_path, args.fixed_k);

    delete g_ensembler;
 
    return 0;
}
#ifndef GREEDY_ALIGN_UTILS_HPP
#define GREEDY_ALIGN_UTILS_HPP


#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

#include <chrono>
#include <thread>
#include <assert.h>
#include <filesystem>


using seqan3::operator""_cigar_operation;

/**
 * @brief Class for benchmark of runtime.
 *        Adopted from https://coliru.stacked-crooked.com/a/508ec779dea5a28d.
 */
template <class DT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer
{
    using timep_t = decltype(ClockT::now());
    
    // variable storing start and end time.
    timep_t _start = ClockT::now();
    timep_t _end = {};

    // total duration between each (_start, _end) pair
    DT _elapsed = DT::zero();

public:
    auto duration() const { 
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time"); 
        return _elapsed; 
    }

    void tick() { 
        _end = timep_t{};
        _start = ClockT::now(); 
    }
    
    void tock() {
        _end = ClockT::now();
        _elapsed += std::chrono::duration_cast<DT>(_end - _start);
    }
    
    float elapsed_seconds() const {
        return ((float) _elapsed.count()) / 1000;
    }
};

/**
 * @brief Check if the specified directory exists. If it does, check if files with certain
 *        file extension name exist. If not, create the directory.
 * @param index_directory a path to a directory to be checked.
 * @param ext a string starting with a ".", indicating file extension name.
 * @returns true if the directory doesn't exist or no such file extension found. False otherwise.
 */
bool check_extension_in(std::filesystem::path const & index_directory,
                            std::string ext) {
    if (!std::filesystem::create_directories(index_directory)) {
        for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
            if (entry.path().extension() == ext) {
                seqan3::debug_stream << "[ERROR]\t\t" << "The file with extension " << ext <<" already exists in directory: " 
                                     << entry.path() << "." << '\n';
                return false;
            }
        }    
    }
    return true;    
}

/**
 * @brief Check if the specified directory exists. If it does, check if files with certain
 *        name exist. If not, create the directory.
 * @param index_directory a path to a directory to be checked.
 * @param filename a string, including the extension name.
 * @returns true if the directory doesn't exist or no such file found. False otherwise.
 */
bool check_filename_in(std::filesystem::path const & index_directory,
                       const std::string& filename) {
    if (!std::filesystem::create_directories(index_directory)) {
        for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
            if (entry.path() == index_directory / filename) {
                seqan3::debug_stream << "[ERROR]\t\t" << "The specified file already exists in directory: " 
                                     << index_directory / filename << "." << '\n';
                return false;
            }
        }    
    }
    return true;    
}

struct _dna4_traits : seqan3::sequence_file_input_default_traits_dna {
    /**
     * @brief Syntax for reading the query file.
     */
    using sequence_alphabet = seqan3::dna4; // instead of dna5
 
    template <typename alph>
    using sequence_container = std::vector<alph>; // must be defined as a template!
};

class CIGAR {
    /**
     * @brief Util class for the cigar string, with support of inserting cigar operation.
     */
private:
    std::vector<seqan3::cigar> CIGAR_vector;

public:
    CIGAR() = default;

    explicit CIGAR(std::vector<seqan3::cigar> CIGAR_init) {
        /**
         * @brief Initialize the class with a vector of seqan3::cigar objects.
         * 
         */
        CIGAR_vector = CIGAR_init;
    }

    /**
     * @brief Initialize the class with `size` number of `op`s.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    CIGAR(unsigned int size, seqan3::cigar::operation op) {
        /**
         * @brief Initialize the class with `size` number of `op`s.
         * 
         */
        for (int i = 0; i < size; i++) {
            seqan3::cigar operation{1, op};
            CIGAR_vector.push_back(operation);
        }
    }

    /**
     * @brief Insert a cigar operation at a certain index.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    void insert(unsigned int index, seqan3::cigar::operation op) {
        seqan3::cigar operation{1, op};
        CIGAR_vector.insert(CIGAR_vector.begin() + index, operation);
    }

    /**
     * @brief Insert one or several cigar operation at the back of all operations.
     * 
     * @param size The number of operations in the CIGAR string.
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc.
     */
    void push_back(unsigned int size, seqan3::cigar::operation op) {
        seqan3::cigar operation{size, op};
        CIGAR_vector.push_back(operation);
    }

    /**
     * @brief Replace the operation at `index` to `op`.
     * 
     * @param index the index at which the operation is replaced
     * @param op the operation, which should be 'M', '=', 'X', 'I', 'D', etc. 
     */
    void replace(unsigned int index, seqan3::cigar::operation op) {
        seqan3::cigar operation{1, op};
        CIGAR_vector[index] = operation;
    }

    /**
     * @brief Convert the vector to a CIGAR string for output.
     * 
     * @return std::string the CIGAR string that captures the whole vector.
     */
    std::string to_string() {
        std::string res;
        unsigned int size = 0;
        seqan3::cigar::operation last_op, curr_op;
        if (!CIGAR_vector.empty()) {
            last_op = get<1>(CIGAR_vector[0]);
            for (auto cigar : CIGAR_vector) {
                auto &[curr_size, curr_op] = cigar;
                if (last_op == curr_op) {
                    size += curr_size;
                } else if (size != 0) {
                    res += std::to_string(size) + last_op.to_char();
                    size = curr_size;
                    last_op = curr_op;
                }
            }
            if (size != 0) {
                res += std::to_string(size) + last_op.to_char();
            }
        }
        return res;
    }
    
    /**
     * @brief Find the k_matching of the indicated alignment based on the CIGAR string.
     * 
     * @return unsigned int the number of k-mers that are matched in the alignment.
     */
    unsigned int k_matching(unsigned int k) {
        unsigned int res = 0;
        int size = 0;
        seqan3::cigar::operation last_op, curr_op;
        if (!CIGAR_vector.empty()) {
            last_op = get<1>(CIGAR_vector[0]);
            for (auto cigar : CIGAR_vector) {
                auto &[curr_size, curr_op] = cigar;
                if (last_op == curr_op) {
                    size += curr_size;
                } else if (size != 0) {
                    if (last_op == '='_cigar_operation && size >= k) {
                        res += size - k + 1;
                    }
                    size = curr_size;
                    last_op = curr_op;
                }
            }
            if (last_op == '='_cigar_operation && size >= k) {
                res += size - k + 1;
            }
        }
        return res;
    }

    std::vector<seqan3::cigar> get_cigar() {
        return CIGAR_vector;
    }
};

/**
 * @brief find the corresponding k-mer given its hash value
 * @param hash the hash value of the k-mer
 * @param k the length of k-mer.
 */
seqan3::dna4_vector hash_to_kmer(unsigned int hash, unsigned int k) {
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> res;
    for (unsigned int i = 0; i < k; i++) {
        auto rank = (hash >> ((k - 1 - i) * 2)) & 3;
        auto nt = seqan3::assign_rank_to(rank, seqan3::dna4{});
        res.push_back(nt);
    }
    return res;
}

#endif
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <zlib.h>
#include <cstring>
#include "sequitur.hpp"
#include "huffman.hpp"
#include <fstream>
#include <string>

using namespace jw;

// Used for bitpacked sequitur serialization.
// (Not in use right now)
class BitWriter {
public:
    void writeBits(uint32_t value, uint8_t numBits) {
        buffer |= (value << bitCount);
        bitCount += numBits;

        while (bitCount >= 8) {
            output.push_back(static_cast<uint8_t>(buffer & 0xFF));
            buffer >>= 8;
            bitCount -= 8;
        }
    }

    void flush() {
        if (bitCount > 0) {
            output.push_back(static_cast<uint8_t>(buffer & 0xFF));
            buffer = 0;
            bitCount = 0;
        }
    }

    const std::vector<uint8_t>& getOutput() const {
        return output;
    }

private:
    uint64_t buffer = 0;
    uint8_t bitCount = 0;
    std::vector<uint8_t> output;
};

// Save any vector of doubles to a CSV file
void write_to_file(const std::vector<double>& data, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }
    for (const auto& val : data) {
        out << val << "\n";
    }
    out.close();
}

// Functions made for re-pair, to read in their files and file sizes.
std::vector<uint8_t> read_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);  // open in binary mode
    if (!file) {
        std::cerr << "Error opening file: " << filename << "\n";
        return {};
    }
    std::vector<uint8_t> data((std::istreambuf_iterator<char>(file)),
                               std::istreambuf_iterator<char>());
    return data;
}

// For re-pair.
void print_combined_compression_ratio(const std::string& label,
                                      size_t original_size,
                                      const std::vector<uint8_t>& prel,
                                      const std::vector<uint8_t>& seq) {
    size_t combined_size = prel.size() + seq.size();
    double ratio = static_cast<double>(combined_size) / original_size;

    std::cout << label << ":\n";
    std::cout << "Prel size: " << prel.size() << " bytes\n";
    std::cout << "Seq size:  " << seq.size() << " bytes\n";
    std::cout << "Combined:  " << combined_size << " bytes\n";
    std::cout << "Compression ratio (combined / original): " << ratio << "\n\n";
}

// Parameters
const double frequency = 1.0;
const int sampling_rate = 500;
const double duration = 2000.0;
const double randomizer_range = 0.15; // Range for semi-random sine
const double max_step = 1.0; // Step size of random walk
const int signal_length = static_cast<int>(sampling_rate * duration);

// Generate time vector
std::vector<double> generate_time_series() {
    std::vector<double> t(signal_length);
    for (int i = 0; i < signal_length; ++i) {
        t[i] = i * (duration / signal_length);
    }
    return t;
}

// Generate sine wave
std::vector<double> generate_sine(const std::vector<double>& t) {
    std::vector<double> amplitude(signal_length);
    for (int i = 0; i < signal_length; ++i) {
        amplitude[i] = std::sin(2.0 * M_PI * frequency * t[i]);
    }
    write_to_file(amplitude, "sine.csv");
    return amplitude;
}

// Add small perturbation
std::vector<double> perturb_sine(const std::vector<double>& sine) {
    std::vector<double> perturbed = sine;
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(-randomizer_range, randomizer_range);

    for (double &val : perturbed) {
        val += dist(gen);
    }
    write_to_file(perturbed, "sine_perturbed.csv");
    return perturbed;
}

// Generate random walk
std::vector<double> generate_random_walk() {
    std::vector<double> walk(signal_length);
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(-max_step, max_step);

    double point = 0.0;
    for (int i = 0; i < signal_length; ++i) {
        point += dist(gen);
        walk[i] = point;
    }
    write_to_file(walk, "random_walk.csv");
    return walk;
}

// Quantize to 0–25
std::vector<uint8_t> quantize_signal(const std::vector<double>& signal, const double symbol_number) {
    std::vector<uint8_t> result(signal_length);
    for (int i = 0; i < signal_length; ++i) {
        // Scale to [0, 25], then shift to [1, 26]
        double scaled = (signal[i] + 1.0) / 2.0 * symbol_number;
        int quantized = static_cast<int>(std::round(scaled)) + 1; // range 1–26
        result[i] = static_cast<uint8_t>(std::clamp(quantized, 1, static_cast<int>(symbol_number)));
    }
    return result;
}

// Convert symbols to alphabet for grammar compression
std::string convert_to_symbols(const std::vector<uint8_t>& quantized) {
    std::string symbols;
    symbols.reserve(quantized.size());

    for (uint8_t val : quantized) {
        // Map 1 → 'A', 2 → 'B', ..., 26 → 'Z'
        if (val >= 1 && val <= 26) {
            symbols += static_cast<char>('A' + val - 1);
        } else {
            symbols += '?'; // fallback for out-of-bounds (shouldn't happen)
        }
    }
    return symbols;
}

//
// Compression schemes (sequitur, re-pair, and zlib)
//

// Compress using zlib
std::vector<uint8_t> compress_zlib(const std::vector<uint8_t>& input) {
    uLongf compressed_size = compressBound(input.size());
    std::vector<uint8_t> compressed(compressed_size);

    if (compress(compressed.data(), &compressed_size, input.data(), input.size()) != Z_OK) {
        throw std::runtime_error("Zlib compression failed");
    }

    compressed.resize(compressed_size);
    return compressed;
}

// Compress using sequitur
Sequitur<char> compress_sequitur(std::string symbols) {

    // Create sequitur instance
    Sequitur<char> s; 
    for (int i = 0; i < symbols.length(); i++) {
        // Create ruleset as symbols are pushed
        s.push_back(symbols[i]);
    }

    return s;
}

// I need a way of converting the objects generated by sequitur into actual bytes
// So I can measure the compression ratio.
std::vector<uint8_t> serialize_sequitur_bitpacked(Sequitur<char>& s) {
    std::vector<uint8_t> output;
    BitWriter writer;

    // Map each rule head to a rule ID
    std::unordered_map<Symbol*, uint8_t> rule_ids;
    uint8_t next_id = 0;

    auto rules = s.getRules();
    for (const auto& entry : rules) {
        assert(next_id < 64 && "Too many rules for 6-bit encoding");
        rule_ids[entry.second] = next_id++;
    }

    // Serialize number of rules (full byte)
    output.push_back(static_cast<uint8_t>(rule_ids.size()));

    // Write each rule: [rule_id][symbol_count][symbols...]
    for (const auto& entry : rules) {
        uint rule_id = entry.first;
        Symbol* head = entry.second;
        Symbol* tail = static_cast<RuleHead*>(head)->getTail();

        // Temporary buffer of 6-bit encoded symbols
        std::vector<uint8_t> symbols;
        for (Symbol* sym = head->next(); sym != tail; sym = sym->next()) {
            if (auto* rule_sym = dynamic_cast<RuleSymbol*>(sym)) {
                uint8_t id = rule_sym->getID();
                symbols.push_back(static_cast<uint8_t>(id + 26)); // rules = 26+
            } else if (auto* val_sym = dynamic_cast<ValueSymbol<char>*>(sym)) {
                char ch = val_sym->getValue();
                assert(ch >= 'A' && ch <= 'Z');
                symbols.push_back(static_cast<uint8_t>(ch - 'A')); // terminals = 0–25
            } else {
                std::cerr << "Unknown symbol type\n";
            }
        }

        // Write rule ID (6 bits)
        writer.writeBits(rule_id + 26, 6);  // stored as non-terminal
        // Write rule length (8 bits)
        output.push_back(static_cast<uint8_t>(symbols.size()));
        // Write symbols (6 bits each)
        for (uint8_t s : symbols) {
            writer.writeBits(s, 6);
        }
    }

    // Write main expansion sequence (same 6-bit encoding)
    Symbol* main_rule_head = s.getRules().at(0);
    Symbol* tail = static_cast<RuleHead*>(main_rule_head)->getTail();

    for (Symbol* sym = main_rule_head->next(); sym != tail; sym = sym->next()) {
        if (auto* rule_sym = dynamic_cast<RuleSymbol*>(sym)) {
            writer.writeBits(rule_sym->getID() + 26, 6);
        } else if (auto* val_sym = dynamic_cast<ValueSymbol<char>*>(sym)) {
            char ch = val_sym->getValue();
            writer.writeBits(ch - 'A', 6);
        } else {
            std::cerr << "Unknown symbol in start sequence\n";
        }
    }

    writer.flush();
    const auto& bitpacked = writer.getOutput();
    output.insert(output.end(), bitpacked.begin(), bitpacked.end());

    return output;
}

std::vector<uint8_t> serialize_sequitur(Sequitur<char>& s) {
    std::vector<uint8_t> output;

    // Map rules to numeric IDs starting at 128 (terminal ASCII chars < 128)
    std::unordered_map<Symbol*, uint8_t> rule_ids;
    uint8_t next_rule_id = 128;

    auto rules = s.getRules();
    for (const auto& entry : rules) {
        rule_ids[entry.second] = next_rule_id++;
    }

    // Serialize number of rules
    output.push_back(static_cast<uint8_t>(rule_ids.size()));

    // Serialize each rule
    for (const auto& entry : rules) {
        uint rule_id = entry.first;
        Symbol* head = entry.second;

        std::vector<uint8_t> body;

        for (Symbol* sym = head->next(); sym != static_cast<jw::RuleHead*>(head)->getTail(); sym = sym->next()) {
            if (auto* rule_sym = dynamic_cast<jw::RuleSymbol*>(sym)) {
                body.push_back(static_cast<uint8_t>(rule_sym->getID() + 128));  // Rule ID
            } else if (auto* val_sym = dynamic_cast<jw::ValueSymbol<char>*>(sym)) {
                body.push_back(static_cast<uint8_t>(val_sym->getValue()));  // Terminal as ASCII
            } else {
                std::cerr << "Unknown symbol type encountered\n";
            }
        }

        output.push_back(static_cast<uint8_t>(rule_id + 128));      // Rule ID
        output.push_back(static_cast<uint8_t>(body.size()));        // Rule length
        output.insert(output.end(), body.begin(), body.end());      // Rule body
    }

    // Serialize the main sequence (start expansion)
    Symbol* main_rule_head = s.getRules().at(0);  // get the head of rule 0
    Symbol* tail = static_cast<RuleHead*>(main_rule_head)->getTail();

    for (Symbol* sym = main_rule_head->next(); sym != tail; sym = sym->next()) {
        if (auto* rule_sym = dynamic_cast<RuleSymbol*>(sym)) {
            output.push_back(rule_ids[rule_sym]);
        } else if (auto* val_sym = dynamic_cast<ValueSymbol<char>*>(sym)) {
            output.push_back(static_cast<uint8_t>(val_sym->getValue()));
        } else {
            std::cerr << "Unknown symbol type in main sequence.\n";
        }
    }

    return output;
}


// Get compression size of ruleset of sequitur
double get_sequitur_compression_size(Sequitur<char> s) {
    // TODO: Implement
    return 0.0;
}

// Debug helper
void print_compression_ratio(const std::string& name, size_t original, size_t compressed) {
    std::cout << name << ":\n"
              << "Original size:   " << original << " bytes\n"
              << "Compressed size: " << compressed << " bytes\n"
              << "Compression ratio: " << static_cast<double>(compressed) / original << "\n\n";
}

int main() {
    auto t = generate_time_series();

    auto sine = generate_sine(t);
    auto semi_rand = perturb_sine(sine);
    auto walk = generate_random_walk();

    // Quantize and compress each signal
    double symbol_number = 100.0;
    auto q_sine = quantize_signal(sine, symbol_number);
    auto q_semi_rand = quantize_signal(semi_rand, symbol_number);
    auto q_walk = quantize_signal(walk, symbol_number);
    
    // Convert quantizations to symbols for sequitur
    auto symbol_sine = convert_to_symbols(q_sine);
    auto symbol_semi_rand = convert_to_symbols(q_semi_rand);
    auto symbol_walk = convert_to_symbols(q_walk);

    // Compress signals using zlib
    auto z_sine = compress_zlib(q_sine);
    auto z_semi_rand = compress_zlib(q_semi_rand);
    auto z_walk = compress_zlib(q_walk);
    // Compress using sequitur
    auto seq_sine = compress_sequitur(symbol_sine);
    auto seq_semi_rand = compress_sequitur(symbol_semi_rand);
    auto seq_walk = compress_sequitur(symbol_walk);
    // Serialize to compare compression ratio
    auto seq_sine_serialized = serialize_sequitur(seq_sine);
    auto seq_semi_rand_serialized = serialize_sequitur(seq_semi_rand);
    auto seq_walk_serialized = serialize_sequitur(seq_walk);
    // Use Huffman coding to further compress data
    // This is working awfully and 10x bigger than the original data. Oops
    auto sine_seq_huffman = huffmanCompress(seq_sine_serialized);
    auto semi_rand_seq_huffman = huffmanCompress(seq_semi_rand_serialized);
    auto walk_seq_huffman = huffmanCompress(seq_walk_serialized);
    // Just use LZ to compress. Let's see what happens
    auto seq_sine_z = compress_zlib(seq_sine_serialized);
    auto seq_semi_rand_z = compress_zlib(seq_semi_rand_serialized);
    auto seq_walk_z = compress_zlib(seq_walk_serialized);
    // Report
    std::cout << "Traditional:\n" << std::endl;
    print_compression_ratio("Sine Wave", q_sine.size(), z_sine.size());
    print_compression_ratio("Semi-Random Sine", q_semi_rand.size(), z_semi_rand.size());
    print_compression_ratio("Random Walk", q_walk.size(), z_walk.size());
    std::cout << "\nGrammar-Based RAW:\n" << std::endl;
    print_compression_ratio("Sequitur Sine", q_sine.size(), seq_sine_serialized.size());
    print_compression_ratio("Sequitur Semi-Random", q_semi_rand.size(), seq_semi_rand_serialized.size());
    print_compression_ratio("Sequitur Random Walk", q_walk.size(), seq_walk_serialized.size());
    std::cout << "\nGrammar-Based ZLIB:\n" << std::endl;
    print_compression_ratio("Sequitur Sine", q_sine.size(), seq_sine_z.size());
    print_compression_ratio("Sequitur Semi-Random", q_semi_rand.size(), seq_semi_rand_z.size());
    print_compression_ratio("Sequitur Random Walk", q_walk.size(), seq_walk_z.size());
    std::cout << "\nGrammar-Based HUFFMAN:\n" << std::endl;
    print_compression_ratio("Sequitur Sine", q_sine.size(), sine_seq_huffman.size());
    print_compression_ratio("Sequitur Semi-Random", q_semi_rand.size(), semi_rand_seq_huffman.size());
    print_compression_ratio("Sequitur Random Walk", q_walk.size(), walk_seq_huffman.size());

    // Now check re-pair (This did not work at all. Seems )
    std::cout << "\nGrammar-Based Re-Pair:\n" << std::endl;
    auto prel_sine = read_file("./repair_outputs/sine.csv.prel");
    auto seq_sine_f  = read_file("./repair_outputs/sine.csv.seq");
    auto prel_semi_rand = read_file("./repair_outputs/sine_perturbed.csv.prel");
    auto seq_semi_rand_f  = read_file("./repair_outputs/sine_perturbed.csv.seq");
    auto prel_walk = read_file("./repair_outputs/random_walk.csv.prel");
    auto seq_walk_f  = read_file("./repair_outputs/random_walk.csv.seq");

    print_combined_compression_ratio("Sine Combined", q_sine.size(), prel_sine, seq_sine_f);
    print_combined_compression_ratio("Sine Perturbed Combined", q_sine.size(), prel_sine, seq_semi_rand_f);
    print_combined_compression_ratio("Walk Combined", q_sine.size(), prel_sine, seq_walk_f);




    // Debug: Print sequitur rulesets
    // std::cout << "Sine rules:" << std::endl;
    // seq_sine.printRules();
    // std::cout << "Random sine rules:" << std::endl;
    // seq_semi_rand.printRules();
    // std::cout << "Random walk rules:" << std::endl;
    // seq_walk.printRules();

    return 0;
}

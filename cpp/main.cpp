#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <zlib.h> // Link with -lz
#include <cstring>

// Parameters
const double frequency = 1.0;
const int sampling_rate = 50;
const double duration = 5.0;
const double randomizer_range = 0.1;
const double max_step = 0.25;
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
    return walk;
}

// Quantize to 0–25
std::vector<uint8_t> quantize_signal(const std::vector<double>& signal) {
    std::vector<uint8_t> result(signal_length);
    for (int i = 0; i < signal_length; ++i) {
        // Scale to [0, 25], then shift to [1, 26]
        double scaled = (signal[i] + 1.0) / 2.0 * 25.0;
        int quantized = static_cast<int>(std::round(scaled)) + 1; // range 1–26
        result[i] = static_cast<uint8_t>(std::clamp(quantized, 1, 26));
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
    auto q_sine = quantize_signal(sine);
    auto q_semi_rand = quantize_signal(semi_rand);
    auto q_walk = quantize_signal(walk);

    auto z_sine = compress_zlib(q_sine);
    auto z_semi_rand = compress_zlib(q_semi_rand);
    auto z_walk = compress_zlib(q_walk);

    // Report
    print_compression_ratio("Sine Wave", q_sine.size(), z_sine.size());
    print_compression_ratio("Semi-Random Sine", q_semi_rand.size(), z_semi_rand.size());
    print_compression_ratio("Random Walk", q_walk.size(), z_walk.size());

    return 0;
}

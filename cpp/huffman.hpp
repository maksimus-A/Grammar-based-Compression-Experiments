#ifndef HUFFMAN_HPP
#define HUFFMAN_HPP

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>
#include <string>

// Node in Huffman tree
struct HuffmanNode {
    uint8_t value;
    size_t freq;
    std::shared_ptr<HuffmanNode> left = nullptr;
    std::shared_ptr<HuffmanNode> right = nullptr;

    bool isLeaf() const;
};

// Comparison operator for priority queue
struct CompareNode {
    bool operator()(const std::shared_ptr<HuffmanNode>& a,
                    const std::shared_ptr<HuffmanNode>& b);
};

// Build Huffman code map from a tree
void buildCodeMap(const std::shared_ptr<HuffmanNode>& node,
                  std::unordered_map<uint8_t, std::string>& codeMap,
                  std::string prefix = "");

// Compress data using Huffman coding
std::vector<uint8_t> huffmanCompress(const std::vector<uint8_t>& input);

#endif // HUFFMAN_HPP

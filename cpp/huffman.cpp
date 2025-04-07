#include "huffman.hpp"
#include <queue>
#include <bitset>
#include <iostream>

bool HuffmanNode::isLeaf() const {
    return !left && !right;
}

bool CompareNode::operator()(const std::shared_ptr<HuffmanNode>& a,
                             const std::shared_ptr<HuffmanNode>& b) {
    return a->freq > b->freq;
}

void buildCodeMap(const std::shared_ptr<HuffmanNode>& node,
                  std::unordered_map<uint8_t, std::string>& codeMap,
                  std::string prefix) {
    if (!node) return;
    if (node->isLeaf()) {
        codeMap[node->value] = prefix.empty() ? "0" : prefix;
        return;
    }
    buildCodeMap(node->left, codeMap, prefix + "0");
    buildCodeMap(node->right, codeMap, prefix + "1");
}

std::vector<uint8_t> huffmanCompress(const std::vector<uint8_t>& input) {
    std::unordered_map<uint8_t, size_t> freq;
    for (uint8_t byte : input) freq[byte]++;

    std::priority_queue<std::shared_ptr<HuffmanNode>,
                        std::vector<std::shared_ptr<HuffmanNode>>,
                        CompareNode> queue;

    for (const auto& [byte, f] : freq) {
        queue.push(std::make_shared<HuffmanNode>(HuffmanNode{byte, f}));
    }

    while (queue.size() > 1) {
        auto left = queue.top(); queue.pop();
        auto right = queue.top(); queue.pop();
        auto parent = std::make_shared<HuffmanNode>(HuffmanNode{0, left->freq + right->freq, left, right});
        queue.push(parent);
    }

    auto root = queue.top();
    std::unordered_map<uint8_t, std::string> codeMap;
    buildCodeMap(root, codeMap);

    std::string bitstream;
    for (uint8_t byte : input) {
        bitstream += codeMap[byte];
    }

    std::vector<uint8_t> encoded;
    uint8_t current = 0;
    int bitsFilled = 0;

    for (char bit : bitstream) {
        current <<= 1;
        if (bit == '1') current |= 1;
        bitsFilled++;
        if (bitsFilled == 8) {
            encoded.push_back(current);
            current = 0;
            bitsFilled = 0;
        }
    }

    if (bitsFilled > 0) {
        current <<= (8 - bitsFilled);
        encoded.push_back(current);
    }

    // Serialize header (simple): number of codes and each [byte][length][bits...]
    std::vector<uint8_t> output;
    uint16_t table_size = static_cast<uint16_t>(codeMap.size());
    output.push_back(table_size >> 8);
    output.push_back(table_size & 0xFF);

    for (const auto& [byte, code] : codeMap) {
        output.push_back(byte);
        output.push_back(static_cast<uint8_t>(code.size()));
        uint8_t buf = 0;
        int filled = 0;
        for (char c : code) {
            buf <<= 1;
            if (c == '1') buf |= 1;
            filled++;
            if (filled == 8) {
                output.push_back(buf);
                buf = 0;
                filled = 0;
            }
        }
        if (filled > 0) {
            buf <<= (8 - filled);
            output.push_back(buf);
        }
    }

    // Append compressed data
    output.insert(output.end(), encoded.begin(), encoded.end());
    return output;
}

import numpy as np
import matplotlib.pyplot as plt
import random
# Used for SAX (floating point to symbols)
from pyts.approximation import SymbolicAggregateApproximation as SAX
from sksequitur import parse # Used for Sequitur
import sys
import zlib # Use to compress data

def utf8len(s):
    # Returns the length of a string 's' in bytes.
    return len(s.encode('utf-8'))

# Just messing around for now to get a handle of things.
# Ideally, I will rewrite all of this in C++ for efficiency.

debug = False # Just to check debugging
show_plot = False
# Parase command line arguments
for i in range(len(sys.argv)):
    if sys.argv[i][0] == '-':
        if sys.argv[i][1:] == "debug" or sys.argv[i][1] == "d":
            debug = True
        if sys.argv[i][1:] == "plot" or sys.argv[i][1] == "p":
            show_plot = True

# Parameters
frequency = 1        # Hz
sampling_rate = 50  # Hz
duration = 5         # seconds
# SAX Parameters (unused rn)
n_segments = 50      # SAX segments

# Generate time values
# TODO: Add real data to compress and test on.
t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

"""
Generating 3 different datasets: Sine wave, semi-random sine wave, and random walk.
"""
# Generate sine wave
amplitude = np.sin(2 * np.pi * frequency * t)

# Generate semi-random sine wave
randomizer_range = 0.1 # modify each sine val by +- randomizer_range
amplitude_semi_rand = amplitude.copy()
for i in range(len(amplitude)):
    perturbation = random.uniform(-randomizer_range, randomizer_range)
    amplitude_semi_rand[i] += perturbation
    
# Generate random walk
random_walk = np.empty_like(amplitude)
point = 0
for i in range(len(random_walk)):
    step = random.uniform(-0.25, 0.25)
    point += step
    random_walk[i] = point

"""
 SAX
"""
# Use on generic sine wave
transformer = SAX()
amplitude_2d = amplitude.reshape(1, -1)
symbolic_sine = transformer.fit_transform(amplitude_2d)
if debug:
    np.set_printoptions(threshold=np.inf)
    print(symbolic_sine)
    
# Use on perturbed sine wave
transformer = SAX()
amplitude_2d = amplitude_semi_rand.reshape(1, -1)
symbolic_perturbed_sine = transformer.fit_transform(amplitude_2d)
if debug:
    print(symbolic_perturbed_sine)
    
# Use SAX on random walk
transformer = SAX()
amplitude_2d = amplitude_semi_rand.reshape(1, -1)
symbolic_random_walk= transformer.fit_transform(amplitude_2d)
if debug:
    print(symbolic_random_walk)

"""
 Compress using Sequitur
"""
# Compress sine wave
symbols = symbolic_sine[0].tolist()  # list of symbols
grammar_sine = parse(symbols)
print(grammar_sine)

# Compress perturbed sine wave
symbols = symbolic_perturbed_sine[0].tolist()  # list of symbols
grammar_perturbed = parse(symbols)
print("perturbed: ", grammar_perturbed)

# Compress random walk
symbols = symbolic_random_walk[0].tolist()  # list of symbols
grammar_random_walk = parse(symbols)
print("random walk: ", grammar_random_walk)

"""
In this section, we use traditional compression methods. We will test the 'zlib'
library's capabilities which uses LZ77 and Huffman to compress.
"""

# Step 1: Quantize the float data
# Convert to 1-byte integers (e.g., 0–255)
quantized_sine = ((amplitude + 1) / 2 * 255).astype(np.uint8)

# Step 2: Convert to bytes
byte_data_sine = quantized_sine.tobytes()

# Step 3: Compress using zlib (uses LZ77 + Huffman)
compressed_sine = zlib.compress(byte_data_sine)

# Compare sizes
print(f"Original size: {len(byte_data_sine)} bytes")
print(f"Compressed size: {len(compressed_sine)} bytes")
print(f"Compression ratio: {len(compressed_sine) / len(byte_data_sine):.2f}")

# Plot sine wave
if show_plot:
    plt.plot(t, amplitude)
    plt.plot(t, amplitude_semi_rand)
    plt.plot(t, random_walk)
    plt.title(f"{frequency} Hz Sine Wave")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.show()


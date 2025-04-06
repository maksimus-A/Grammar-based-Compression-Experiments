import numpy as np
import matplotlib.pyplot as plt
import random
# Used for SAX (floating point to symbols)
from pyts.approximation import SymbolicAggregateApproximation as SAX
from sksequitur import parse # Used for Sequitur
import sys

# Just messing around for now to get a handle of things.
# Ideally, I will rewrite all of this in C++ for efficiency.

debug = False # Just to check debugging
show_plot = False
for i in range(len(sys.argv)):
    if sys.argv[i][0] == '-':
        if sys.argv[i][1:] == "debug":
            debug = True
        if sys.argv[i][1:] == "plot":
            show_plot = True

# Parameters
frequency = 1        # Hz
sampling_rate = 50  # Hz
duration = 1         # seconds
# SAX Parameters (unused rn)
n_segments = 50      # SAX segments

# Generate time values
t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

# Generating 3 different datasets: Sine wave, semi-random sine wave, and random walk.

# Generate sine wave
amplitude = np.sin(2 * np.pi * frequency * t)

# Generate semi-random sine wave
randomizer_range = 0.1 # modify each sine val by +- randomizer_range
amplitude_semi_rand = amplitude.copy()
for i in range(len(amplitude)):
    perturbation = random.uniform(-randomizer_range, randomizer_range)
    amplitude_semi_rand[i] += perturbation
    

# SAX

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

# Compress using Sequitur

# Compress sine wave
symbols = symbolic_sine[0].tolist()  # list of symbols
grammar = parse(symbols)
print(grammar)

# Compress perturbed sine wave
symbols = symbolic_perturbed_sine[0].tolist()  # list of symbols
grammar_perturbed = parse(symbols)
print("perturbed: ", grammar_perturbed)


# Plot sine wave
if show_plot:
    plt.plot(t, amplitude)
    plt.plot(t, amplitude_semi_rand)
    plt.title(f"{frequency} Hz Sine Wave")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.show()


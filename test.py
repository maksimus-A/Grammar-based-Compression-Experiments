from sksequitur import parse
import numpy as np
import matplotlib.pyplot as plt
# Used for SAX (floating point to symbols)
from pyts.approximation import SymbolicAggregateApproximation as SAX
from sksequitur import parse # Used for Sequitur
import sys

# Just messing around for now to get a handle of things.
# Ideally, I will rewrite all of this in C++ for efficiency.

debug = False # Just to check debugging
for i in range(len(sys.argv)):
    if sys.argv[i][0] == '-':
        if sys.argv[i][1:] == "debug":
            debug = True

# Parameters
frequency = 1        # Hz
sampling_rate = 500  # Hz
duration = 5         # seconds
n_segments = 50      # SAX segments

# Generate time values
t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

# Generate sine wave
amplitude = np.sin(2 * np.pi * frequency * t)

# SAX
transformer = SAX()
amplitude_2d = amplitude.reshape(1, -1)
symbolic_sine = transformer.fit_transform(amplitude_2d)
if debug:
    np.set_printoptions(threshold=np.inf)
    print(symbolic_sine)

# Compress using Sequitur
symbols = symbolic_sine[0].tolist()  # list of symbols
grammar = parse(symbols)
print(grammar)


# Plot sine wave
if debug:
    plt.plot(t, amplitude)
    plt.title(f"{frequency} Hz Sine Wave")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.show()


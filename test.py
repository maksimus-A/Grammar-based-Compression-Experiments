from sksequitur import parse
import numpy as np
import matplotlib.pyplot as plt
# Used for SAX (floating point to symbols)
from pyts.approximation import SymbolicAggregateApproximation as SAX

# Just messing around for now to get a handle of things.
# Ideally, I will rewrite all of this in C++ for efficiency.


# Parameters
frequency = 2        # Hz
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
np.set_printoptions(threshold=np.inf)
print(symbolic_sine)

# Plot sine wave
# plt.plot(t, amplitude)
# plt.title(f"{frequency} Hz Sine Wave")
# plt.xlabel("Time [s]")
# plt.ylabel("Amplitude")
# plt.grid(True)
# plt.show()


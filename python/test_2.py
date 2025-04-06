import numpy as np
import matplotlib.pyplot as plt
from pyts.approximation import SymbolicAggregateApproximation

# Parameters
frequency = 2        # Hz
sampling_rate = 500  # Hz
duration = 5         # seconds
n_segments = 20      # SAX segments
alphabet_size = 7    # Number of symbols

# Time and signal
t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)
amplitude = np.sin(2 * np.pi * frequency * t)

# SAX transform
# Define transformer
window_size = 50   # number of time points per segment
transformer = SymbolicAggregateApproximation(n_bins=n_segments)

# Transform the signal
amplitude_2d = amplitude.reshape(1, -1)
symbolic_sine = transformer.fit_transform(amplitude_2d)[0]
symbolic_sine = transformer.fit_transform(amplitude.reshape(1, -1))[0]

# Plot
segment_length = len(t) // n_segments
symbol_times = np.arange(n_segments) * segment_length

plt.plot(t, amplitude, label="Original Sine Wave", alpha=0.5)
for i, symbol in enumerate(symbolic_sine):
    idx = min(i * segment_length, len(t) - 1)
    x_pos = t[idx]
    plt.text(x_pos, 0.9, symbol, fontsize=10, color='red')

plt.title("SAX Symbol Overlay on Sine Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend()
plt.show()

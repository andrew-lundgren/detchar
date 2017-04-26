# Code Development

# Filtering

FIR filtering is more flexible and plenty efficient in offline applications. We should focus on that but also try to allow conversion of Foton filters to Python, and support ZPK and SOS forms of IIR filtering.

All of these should have serious tests on fake data (or real when applicable). They should be written in a way that they could be contributed to scipy if at all practicable.

### FIR filtering

* Efficient application of FIR filters - break into largest blocks that can be efficiently FFTed, walk through accounting for filter corruption at edges of blocks
* Conversion of spectrum to whitening filter (should deal nicely with badness near DC and Nyquist)
* Conversion of zero-phase to causal minimum phase filter (uses cepstrum - Andy has code)
* Creation of derivative filters with bandpass (remember factors of i)
* Do we need a good Hilbert transform filter?
* Notch design (use height and width of peak to design filter)

We might want a good way to indicate when a filter is zero-phase (either centered or wrapped around), or causal. This probably affects what is corrupted and how to best align the output. If we use overlap-add or overlap-save the overlapping of the blocks isn't such an issue.

### Applications of filtering

* Reliable methods for measuring coherence and transfer functions even with high dynamic range (use ASD2 method - Andy has code)
* Subtract noise from signals using transfer function
* Line tracking and removal

### IIR filters

* Foton to scipy converter, in pure Python (if possible)
* IIR to FIR filter
* Well-tested Bode plot that can be called in a single line (for FIR or IIR filters); must make frequency convention (angular or normal, how to set Nyquist) very clear

### Code snippets

```
	# Take the cepstrum - IFFT of log spectrum
    cepstrum=irfft(log(invasd))
    # Fold the anti-causal part over the causal part
    halfidx=len(cepstrum)/2
    folded=concatenate([[cepstrum[0]],
    					cepstrum[1:halfidx]+cepstrum[-1:-halfidx:-1],
    					zeros(halfidx)])
    # Make sure the response truncates smoothly in TD
    folded[:halfidx]*=signal.hann(2*halfidx)[-halfidx:]
    # Get the minimum-phase filter with the same amplitude
    # as the original inverse ASD
    newinvasd=exp(rfft(folded))
    # The time-domain version is easiest to work with
    # because you can pad with zeros to filter a longer segment of data
    return irfft(newinvasd)
```


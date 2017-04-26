# Detchar Desiderata

These are the things we need for a more advanced statistical / machine learning approach to detchar.

## What we want to predict

1. Times (precise to ~ 0.1 seconds) when a particular kind of glitch is likely; also predict how likely, what amplitude
2. Segments of time (precise to ~ 1 second) where the noise is nonstationary or non-Gaussian (and how nonstationary, and in which bands)
3. Segments of time when the background of a search has a particular distribution; determine accurately the shape of the distribution
4. Causes of locklosses

Included in the first three will be some clue to _why_ these kinds of noise happened. We may have to reconstruct that from which auxiliary channels are used, and it will be an interesting problem on its own for complicated models.

### Glitches

The important point is that it's unlikely that we can predict most glitches with a high degree of confidence. All our measurements are noisy and incomplete, so we are sometimes happy with even a 90% false alarm probability of predicting the worst kinds of glitches. Very bad glitches are hopefully the tails of a distribution, and not just black swans (happen once with no way of predicting).

### Non-idealized noise

The noise may be stationary but non-Gaussian, which shouldn't be too hard to deal with. The noise could have heavier tails, or the whitened residuals may be uncorrelated but not independent. There's some interesting work on detecting causality/the arrow of time from non-Gaussian residuals.

But the simplest may be to use ARIMA, and compare with causal rather than zero-phase whitening, and look at the residuals for some of the channels, parameterizing deviations from Gaussianity and independence.

Much more difficult is non-stationary noise. There's a literature on localized stationarity. Just doing some kind of robust PCA on the spectrogram and finding lumps that move together is a pretty good start. The robustness is needed because of two kinds of outliers: loud glitches distorting a single spectrum, and lines polluting a particular frequency band.


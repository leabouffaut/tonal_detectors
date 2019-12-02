# tonal_detectors
Tonal detectors for low frequency vocalization of whales
[![DOI](https://zenodo.org/badge/182302880.svg)](https://zenodo.org/badge/latestdoi/182302880)

Methods are described and compared for the detection of low-frequency whale signals (Antarctic Blue Whale) in: 
L. Bouffaut, S. Madhusudhana, V. Labat, A. Boudraa and H. Klinck (In prep.), “Comparative study of tonal detectors for
low frequency vocalizations of blue whales,” prepared for J. Acoust. Soc. Am.

# Methods implemented
- Instantaneous frequency estimator, based on Boashash, "Estimating and interpreting the instantaneous frequency of a signal. II. algorithms and applications," Proc. of the IEEE 80(4), 540-568 (1992) doi: 10.1109/5.135378.

- YIN estimator, based on A. De Cheveigné and H. Kawahara, "YIN, a fundamental frequency estimator for speech and music," J. Acoust. Soc. Am. 111(4), 1917-1930 (2002) doi: 10.1121/1.1458024.

- Harmonic product spectrum, A. M. Noll, "Pitch determination of human speech by the harmonic product spectrum, the harmonic sum spectrum, and a maximum likelihood estimate," in Symposium on Computer Processing in Communication, ed., University of Brooklyn Press, New York, Vol. 19, pp. 779-797 (1969).

- Cost-function-based detector, based on M. F. Baumgartner and S. E. Mussoline, "A generalized baleen whale call detection and classification system," J. Acoust. Soc. Am. 129(5), 2889-2902 (2011) doi: 10.1121/1.3562166.

- Ridge detector, based on S. K. Madhusudhana, "Automatic detectors for underwater soundscape measurements," Ph.D. thesis, Curtin University, 2015.

* An additional program is added to allow SNR-control on simulations.

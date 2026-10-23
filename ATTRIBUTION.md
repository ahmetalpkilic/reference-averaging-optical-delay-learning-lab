# Attribution and source review

Adapted lab experiment and documentation, and reused model.py, optics.py, whitening.py, 26 tests, runner and requirements from ahmetalpkilic/shared-reference-optical-delay-learning-lab at a4f6ea63fb3a72365261ddcc423923984f1d9bb3 (MIT). shared_reference.py preserves that commit's lab.py unchanged. LICENSE retains the author's notice; PREDECESSOR_ATTRIBUTION.md preserves the earlier chain. New work: independent reference-bank averaging, paired L experiment, covariance derivation/check and five tests. All data synthetic.

Reviewed live public, unarchived repositories on October 23, 2026. Exact commits, activity and license URLs are in sources.json; full license texts are in provenance/.

- brandondube/prysm (MIT): physical-optics reference; latest reviewed work uses signed-distance segment apertures. No copied implementation or runtime dependency; this lab uses an analytic Gaussian pupil.
- scipy/scipy (BSD-3-Clause): ndtr and least_squares runtime dependencies. Recent ShortTimeFFT short-signal padding fix inspected in commit metadata.
- numpy/numpy (BSD-3-Clause plus bundled notices): RNG, arrays and linear algebra. Recent GitHub Actions maintenance.
- matplotlib/matplotlib (PSF-based Matplotlib license): plotting. Recent pandas-date test maintenance.
- ahmetalpkilic/shared-reference-optical-delay-learning-lab (MIT): direct reused predecessor, one shared reference per batch.

Upstream activity records do not imply that these changes are used by the pinned runtime versions. No numerical-library implementation was copied. Dependency licenses and bundled notices remain applicable; this project's MIT license does not relicense them.

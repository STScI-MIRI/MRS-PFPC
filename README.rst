MIRI MRS Point Fixed Pattern Correction (MRS-PFPC)
==================================================

Goal is to derive static 1D fixed pattern correction vectors
for point sources observed with the standard 4 point dither pattern.
Use calibration observations for flux standards and asteroids to 
derive the average fixed pattern correction for each channel and grating.

In Development!
---------------

Active development.
Data still in changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Details
-------

Data reduced using `utils/proc_mrs.py`.  Details of settings are in 
`utils/mrs_helpers.py`.  Leak correction on individual dither spectra
using `utils/leakcor_files.py`.

Models downloaded from CALSPEC and stored in `models/` subdir.

Average residual fringe corrections for each channel and grating
setting created using `plotting/plot_spectra_dithers.py`.

Once the static residual fringe correction is done, then data can 
corrected using `mrs_srfcor` and plotted using `mrs_plot`.

import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.modeling import models, fitting


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # S/N regions
    snreg = {"1short": [5.28, 5.35], "1medium": [6.0, 6.1], "1long": [7.1, 7.25],
             "2short": [8.3, 8.4], "2medium": [9.1, 9.3], "2long": [10.9, 11.1],
             "3short": [12.8, 13.0], "3medium": [14.3, 14.5], "3long": [15.6, 16.0],
             "4short": [19.2, 19.7], "4medium": [21.5, 22.0], "4long": [25.0, 26.0],
             }

    names = ["HD163466_c1",
             "HD163466_c2_e1",
             "HD163466_c2_e2",
             "HD163466_c2_e3",
             # "HD163466_c2_e4",
             "HD163466_c2_e5",
             "HD163466_c2_e6",
             # "HD163466_c2_e7",
             "HD163466_c2_e8",
             "HD163466_c2_e9",
             "HD163466_c2_e10",
             "HD163466_c2_e11",
             ]
    nobs = len(names)

    fontsize = 16

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    for i in range(3):  # channels
        otab = QTable()
        # for j in range(3):  # grating settings
        for gr in ["short", "medium", "long"]:

            allwave = None
            dithave = np.zeros(nobs)
            for k, sname in enumerate(names):
                cfile = f"{sname}/{sname}_pfpc_ch{i+1}-{gr}_x1d.fits"

                pipefile = cfile.replace("pfpc", "level3")
                # get details of segment so the right color can be used
                h = fits.getheader(pipefile, hdu=1)
                chn = int(h["CHANNEL"])
                band = h["BAND"].lower()
                if band == "short":
                    bnum = 0
                elif band == "medium":
                    bnum = 1
                else:
                    bnum = 2
                # pcol = pcolors[(chn - 1) * 3 + bnum]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=u.UnitsWarning)
                    itab = QTable.read(cfile, hdu=1)
                    pipetab = QTable.read(pipefile, hdu=1)

                if allwave is None:
                    allwave = itab["WAVELENGTH"].value
                    nwaves = len(allwave)
                    allspec = np.full((nwaves, nobs), np.nan)
                    allspec_rf = np.full((nwaves, nobs), np.nan)

                allspec[:, k] = itab["FLUX"].value
                allspec_rf[:, k] = itab["RF_FLUX"].value

                dithave[k] = np.nanmedian(allspec[:, k])

            dave = np.average(dithave)

            for k in range(nobs):
                allspec[:, k] *= dave / dithave[k]
                allspec_rf[:, k] *= dave / dithave[k]

                pflux = allwave * allwave * allspec[:, k]
                ax.plot(allwave, pflux)

                pflux = allwave * allwave * allspec_rf[:, k]
                ax.plot(allwave, pflux + 3.0)

            sigfac = 4.0
            stdfunc = "mad_std"
            grow = None

            specclipped = sigma_clipped_stats(
                allspec, axis=1, sigma=sigfac, stdfunc=stdfunc, grow=grow
            )
            specclipped_rf = sigma_clipped_stats(
                allspec_rf, axis=1, sigma=sigfac, stdfunc=stdfunc, grow=grow
            )

            avespec = specclipped[0]
            pflux = allwave * allwave * avespec
            ax.plot(allwave, pflux + 1.0, "k-")
 
            avespec_rf = specclipped_rf[0]
            pflux_rf = allwave * allwave * avespec_rf
            ax.plot(allwave, pflux_rf + 4.0, "k-")

            ckey = f"{chn}{band}"
            if ckey in snreg.keys():
                fit = fitting.LinearLSQFitter()
                line_init = models.Linear1D()
                gvals = (allwave >= snreg[ckey][0]) & (allwave <= snreg[ckey][1])

                # final
                fitted_line = fit(line_init, allwave[gvals], pflux_rf[gvals])
                tratio = pflux_rf[gvals] / fitted_line(allwave[gvals])
                sstats_rf = sigma_clipped_stats(tratio)

                # before residual fringe
                fitted_line = fit(line_init, allwave[gvals], pflux[gvals])
                tratio = pflux[gvals] / fitted_line(allwave[gvals])
                sstats = sigma_clipped_stats(tratio)

                print(ckey, sstats_rf[0]/ sstats_rf[2], sstats[0]/ sstats[2])

    fig.tight_layout()

    save_str = "tmp"
    # save_str = f"{args.objname}/{args.objname}_dither_divide{extstr}_chn{channame}"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
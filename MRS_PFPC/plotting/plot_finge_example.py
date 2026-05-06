import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits

from MRS_PFPC.utils.helpers import pcolors, rebin_constres, get_h_waves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--showchan4", help="show channel 4 with other channels", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    sname = "HD163466_c1"

    files = []
    for ch in range(4):
        for gr in ["short", "medium", "long"]:
            files.append(f"HD163466_c1_nofringe/{sname}_level3_ch{ch+1}-{gr}_x1d.fits")

    # fmt: on
    for k, pipefile in enumerate(files):
        pipefile_wfringe = pipefile.replace("_nofringe", "")
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
        pcol = pcolors[(chn - 1) * 3 + bnum]

        if (not args.showchan4) & (chn == 4):
            continue

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=u.UnitsWarning)
            pipetab = QTable.read(pipefile, hdu=1)
            pipetab_wfringe = QTable.read(pipefile_wfringe, hdu=1)
        pipeflux = pipetab["FLUX"].value
        pipeflux_wfringe = pipetab_wfringe["FLUX"].value
        pipewave = pipetab["WAVELENGTH"].value

        tpipeflux = pipeflux * np.power(pipewave, 2.0)
        ax.plot(pipewave, tpipeflux, linestyle="-", color=pcol, alpha=0.8)

        offval = 1.0
        tpipeflux_wfringe = pipeflux_wfringe * np.power(pipewave, 2.0)
        ax.plot(pipewave, tpipeflux_wfringe + offval, linestyle="-", color=pcol, alpha=0.8)

        if chn < 4:
            yrange = ax.get_ylim()

    yrange = np.array(yrange)
    yrange[1] = yrange[1] + offval
    ax.set_ylim(yrange)
    xrange = np.array(ax.get_xlim())
    xrange[0] = 3.5
    ax.set_xlim(xrange)

    ttitle = sname
    ax.set_title(ttitle)
    ax.set_xlabel(r"$\lambda$ [$\mu$m]")
    ax.set_ylabel(r"$\lambda^2 F(\nu)$ [$\mu$m$^2$ Jy = RJ units]")

    ax.text(
        3.8, 11.0, "w/ Fringe Flat", fontsize=0.8 * fontsize, rotation=45.0, alpha=0.6
    )

    plt.tight_layout()

    fname = f"figs/full_fringe_example"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()

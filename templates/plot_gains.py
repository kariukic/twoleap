import click
import h5py
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ps_eor import psutil
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm, Normalize

mpl.style.use("default")
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["image.interpolation"] = "none"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.aspect"] = "auto"


def selectPol(sol, pol):
    try:
        if pol == "XX":
            sol_pol = sol[..., 0]
        elif pol == "XY":
            sol_pol = sol[..., 1]
        elif pol == "YX":
            sol_pol = sol[..., 2]
        elif pol == "YY":
            sol_pol = sol[..., 3]

    except IndexError:
        assert pol in ["XX", "YY"]
        if pol == "XX":
            sol_pol = sol[..., 0]
        elif pol == "YY":
            sol_pol = sol[..., 1]

    return sol_pol


def readAnts(sol):
    ant = []
    for a in range(sol["antenna"].shape[0]):
        ant.append(re.split("b|'", str(sol["antenna"][a][0]))[-2])

    return ant


def readSource(sol):
    source = []
    for s in range(sol["source"].shape[0]):
        source.append(str(sol["source"][s][0])[3:-2])

    return source


def ref_phase(sol_ph_noref):
    sol_ph = np.zeros(sol_ph_noref.shape, dtype=np.float64)
    for i in range(sol_ph.shape[2]):
        sol_ph[..., i] = sol_ph_noref[..., i] - sol_ph_noref[..., 0]

    sol_ph[sol_ph > np.pi] -= 2 * np.pi
    sol_ph[sol_ph < -np.pi] += 2 * np.pi

    return sol_ph


def get_ps(d, dx, w="blackmanharris"):
    from scipy.signal import get_window

    w = get_window(w, d.shape[2])[None, None, :]

    ps = abs(np.fft.fftshift(np.fft.fft(d * w, axis=2), axes=2)) ** 2
    M = ps.shape[2]
    if psutil.is_odd(M):
        ps = 0.5 * (ps[:, :, M // 2 + 1 :] + ps[:, :, : M // 2][:, :, ::-1])
    else:
        ps = 0.5 * (ps[:, :, M // 2 + 1 :] + ps[:, :, 1 : M // 2][:, :, ::-1])

    delay = (-np.fft.fftfreq(d.shape[2], dx)[d.shape[2] // 2 + 1 :])[::-1]

    return delay, ps


def readSols(h5list, pols: list = ["XX", "XY", "YX", "YY"]):
    directions = readSource(h5py.File(h5list[0], "r")["sol000"])
    print(directions)

    gains = {}
    for d, sel_dir in enumerate(directions):
        gains[sel_dir] = {}
        gains[sel_dir]["amplitude"] = {}
        gains[sel_dir]["phase"] = {}

        for pol in pols:
            count = 0
            for h5file in h5list:

                sol = h5py.File(h5file, "r")["sol000"]

                sol_amp = selectPol(sol["amplitude000"]["val"], pol)[..., d]

                sol_ph_noref = selectPol(sol["phase000"]["val"], pol)[..., d]
                sol_ph = ref_phase(sol_ph_noref)

                # shape = (time, freq, ant)
                if count == 0:
                    all_sol_amp = np.copy(sol_amp)
                    all_sol_ph = np.copy(sol_ph)

                else:
                    all_sol_amp = np.concatenate([all_sol_amp, sol_amp], axis=0)
                    all_sol_ph = np.concatenate([all_sol_ph, sol_ph], axis=0)

                count += 1

            gains[sel_dir]["amplitude"][pol], gains[sel_dir]["phase"][pol] = (
                all_sol_amp,
                all_sol_ph,
            )

    return gains


def plotAntSols(
    sols,
    time,
    freq,
    ants,
    direction=None,
    component="amplitude",
    vmin=None,
    vmax=None,
    norm="log",
):
    if norm == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(24, 16))

    # Corrected: Use full antenna list for grid dimensions
    nrows, ncols = int(np.ceil(np.sqrt(len(ants)))), int(np.ceil(np.sqrt(len(ants))))
    gs = GridSpec(nrows, ncols, figure=fig)

    a = 0
    for y in range(nrows):
        for x in range(ncols):
            if a < len(ants):  # Corrected: Check against full antenna list
                ax = fig.add_subplot(gs[y, x])
                if len(sols.shape) > 2:
                    im = ax.imshow(
                        sols[:, :, a].T,
                        origin="lower",
                        aspect="auto",
                        interpolation="none",
                        cmap="viridis",
                        extent=[np.min(time), np.max(time), np.min(freq), np.max(freq)],
                        norm=norm,
                    )

                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)

                    if x == 0:
                        ax.set_ylabel("Frequency (MHz)")
                        ax.yaxis.set_tick_params(labelleft=True)
                else:
                    im = ax.plot(
                        sols[:, a].T,
                    )

                ax.set_title(ants[a])

                if y == nrows - 1:
                    ax.set_xlabel("Time (hr)")
                    ax.xaxis.set_tick_params(labelbottom=True)

            a += 1

    if len(sols.shape) > 2:
        fig.subplots_adjust(right=0.90, top=0.95)
        cbar_ax = fig.add_axes([1.01, 0.05, 0.02, 0.92])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label(f"Gain {component}", rotation=90, fontsize="x-large")

    fig.suptitle(f"{direction}")
    fig.tight_layout()
    return fig


def readGaincalSols(h5list, sel_dir=None, pols=["XX", "YY"]):
    gains = {}
    gains[sel_dir] = {}
    gains[sel_dir]["amplitude"] = {}
    gains[sel_dir]["phase"] = {}

    for pol in pols:
        count = 0
        for h5file in h5list:
            sol = h5py.File(h5file, "r")["sol000"]
            # print(sol["amplitude000"]["val"].shape)
            sol_amp = selectPol(sol["amplitude000"]["val"], pol)
            sol_ph_noref = selectPol(sol["phase000"]["val"], pol)
            sol_ph = ref_phase(sol_ph_noref)
            # shape = (time, freq, ant)
            if count == 0:
                all_sol_amp = np.copy(sol_amp)
                all_sol_ph = np.copy(sol_ph)

            else:
                all_sol_amp = np.concatenate([all_sol_amp, sol_amp], axis=1)
                all_sol_ph = np.concatenate([all_sol_ph, sol_ph], axis=1)

            count += 1

        gains[sel_dir]["amplitude"][pol], gains[sel_dir]["phase"][pol] = (
            all_sol_amp,
            all_sol_ph,
        )

    return gains


@click.group()
def main():
    """plot DP3 solutions ..."""


@main.command("plot")
@click.argument(
    "h5list",
    nargs=-1,  # Accepts unlimited positional arguments
    type=click.Path(exists=True),
    required=False,
)
@click.option("--sols_list", help="imagelist", type=str, default="")
@click.option("--filename", help="output filename", type=str, default="gain_sols")
@click.option(
    "--vmin",
    help="vmin",
    type=float,
)
@click.option(
    "--vmax",
    help="vmax",
    type=float,
)
def plot(h5list=[], sols_list="", filename="gains_sols", vmin=None, vmax=None):
    if sols_list:
        with open(sols_list) as t:
            h5list = [line.strip() for line in t]
        if not h5list:
            print("No solutions found in the provided gains solution files list.")
            return
    else:
        assert h5list

    print(f"Processing {len(h5list)} files")

    directions = readSource(h5py.File(h5list[0], "r")["sol000"])
    ants = readAnts(h5py.File(h5list[0], "r")["sol000"])
    ncs = sum(map(lambda x: 1 if x.startswith("CS") else 0, ants))
    nrs = sum(map(lambda x: 1 if x.startswith("RS") else 0, ants))
    print(nrs)

    time_all = [
        list(h5py.File(h, "r")["sol000"]["amplitude000"]["time"]) for h in h5list
    ]
    time_all = [t for tchunk in time_all for t in tchunk]
    time_all = np.asarray([t - time_all[0] for t in time_all]) / 3600

    freq = [
        h5py.File(h, "r")["sol000"]["amplitude000"]["freq"][:] / 1.0e6 for h in h5list
    ]
    freq = np.squeeze(freq)
    print(freq.shape, time_all.shape)

    try:
        gains = readSols(h5list, pols=["XX", "YY"])
    except:
        gains = readGaincalSols(h5list, sel_dir="OINTIN")
        print(gains["OINTIN"]["amplitude"]["XX"].shape, "mainshape")

    # gains = readSols(sols, pols=["XX", "YY"])
    print(gains.keys(), "gains keys")

    comp = "amplitude"
    for d, sel_dir in enumerate(list(directions)):
        fig, axs = plt.subplots(
            nrows=2, ncols=2, dpi=200, sharex=True, sharey=True, figsize=(10, 6)
        )

        for p, pol in enumerate(["XX", "YY"]):
            gns_full = np.squeeze(gains[sel_dir][comp][pol])
            gns = gns_full.copy()
            shape = gns.shape
            print(shape, "shape")
            if len(shape) > 2:
                gns = np.nanmean(gns, axis=0)

            # for idx in idx_to_flag:
            #     gns[:, idx] =  np.nan

            gns_core = gns[..., :ncs]
            gns_remote = gns[..., ncs:]

            axs[0, p].plot(freq, gns_core, lw=1.3, c="gray")
            axs[0, p].set_title(f"{pol} CS")

            axs[1, p].plot(freq, gns_remote, lw=1.3, c="gray")
            axs[1, p].set_title(f"{pol} RS")

            axs[1, p].set_xlabel("Frequency [MHz]")
            # axs[1, p].set_yscale("log")

        dirname = f"{sel_dir.strip('_FIELD')}"
        if dirname != "OINTIN":
            fig.suptitle(dirname)
        fig.tight_layout()
        plt.savefig(
            f"{filename}_{sel_dir.strip('_FIELD')}_{comp}_1d.png",
            bbox_inches="tight",
            dpi=200,
        )

        fig, axs = plt.subplots(
            nrows=2, ncols=1, figsize=(14, 10), dpi=200, sharey=True
        )
        for p, (pol, ax) in enumerate(zip(["XX", "YY"], axs)):
            gns_full = np.squeeze(gains[sel_dir][comp][pol])
            gns = gns_full.copy()
            if len(gns.shape) > 2:
                gns = np.nanmean(gns, axis=0)

            im = ax.imshow(
                gns,
                extent=(
                    0.0,
                    float(len(ants)),
                    float(np.min(freq)),
                    float(np.max(freq)),
                ),
            )
            ax.set_xticks(range(len(ants)))
            ax.set_xticklabels(ants, rotation=90)
            ax.set_title(pol)
            plt.colorbar(im, ax=ax)
            ax.set_ylabel("Frequency (MHz)")
        fig.tight_layout()
        plt.savefig(
            f"{filename}_{sel_dir.strip('_FIELD')}_{comp}_ant_vs_freq.png",
            bbox_inches="tight",
            dpi=200,
        )

        # vmin, vmax = (0.5, 1.5)
        fig = plotAntSols(
            gns,
            time_all,
            freq,
            ants,
            direction=sel_dir,
            # vmin=vmin,
            # vmax=vmax,
            component=comp,
        )
        fig.savefig(
            f"{filename}_{sel_dir.strip('_FIELD')}_{comp}_dynspec.png",
            bbox_inches="tight",
            dpi=200,
        )

        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), dpi=200, sharey=True)
        for p, pol in enumerate(["XX", "YY"]):
            gns = np.squeeze(gains[sel_dir][comp][pol])
            if len(gns.shape) > 2:
                gns = np.swapaxes(gns, 1, 2)
            else:
                gns = np.swapaxes(gns, 0, 1)
                gns = np.expand_dims(gns, axis=0)

            delay, dde_cal_di_ps_core = get_ps(
                gns[:, :48, :], 2452.2
            )  # 36.6)  # 0.196)
            dmax_core = np.nanmax(dde_cal_di_ps_core)

            delay, dde_cal_di_ps_rem = get_ps(gns[:, 48:, :], 0.196)  # 2452.2)  # 36.6)
            dmax_rem = np.nanmax(dde_cal_di_ps_rem)

            axs[p].plot(
                delay,
                np.nanmean(dde_cal_di_ps_core, axis=0).T / dmax_core,
                c="g",
                alpha=0.1,
            )
            axs[p].plot(
                delay,
                np.nanmean(dde_cal_di_ps_core, axis=(0, 1)) / dmax_core,
                c="g",
                alpha=1,
                label="Mean (Core)",
            )

            axs[p].plot(
                delay,
                np.nanmean(dde_cal_di_ps_rem, axis=0).T / dmax_rem,
                c="b",
                alpha=0.1,
            )
            axs[p].plot(
                delay,
                np.nanmean(dde_cal_di_ps_rem, axis=(0, 1)) / dmax_rem,
                c="b",
                alpha=1,
                label="Mean (Remote)",
            )
            axs[p].set_xlabel("Delay (µs)")
            axs[p].grid()
            axs[p].set_yscale("log")
            axs[p].legend()
            axs[p].set_ylim(1e-15, 2)

            axs[p].axhline(1.0, color="k", linestyle="--", lw=0.5)

        axs[0].set_ylabel("Normalised gain spectra")
        fig.tight_layout()
        plt.savefig(
            f"{filename}_{sel_dir.strip('_FIELD')}_{comp}_delay.png",
            bbox_inches="tight",
            dpi=200,
        )


if __name__ == "__main__":
    main()

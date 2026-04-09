import click
import pickle
import numpy as np
import aoquality as ao
import matplotlib.pyplot as plt


class QualityStats(object):

    def __init__(self, qs_file, name, plot_crosses=False):
        self.qs_file = qs_file
        self.name = name
        self.plot_crosses = plot_crosses
        self.pols = ["XX", "XY", "YX", "YY"]

    def pol_idx(self):
        return [0, 1, 2, 3] if self.plot_crosses else [0, 3]

    def get_available_stats(self):
        return {
            "Mean": {"label": "Vis. Mean", "xfactor": 1, "units": "[Jy]", "poly": 3},
            "SNR": {"label": "Vis. SNR", "xfactor": 1, "units": "", "poly": 3},
            "Std": {"label": "Vis Std", "xfactor": 1, "units": "[Jy]", "poly": 3},
            "DStd": {"label": "Vis. Dstd", "xfactor": 1, "units": "[Jy]", "poly": 3},
            "RFIPercentage": {
                "label": "Flagged vis.",
                "xfactor": 1e2,
                "units": r"[%]",
                "poly": 1,
            },
            "Count": {"label": "Vis. count", "xfactor": 1, "units": "", "poly": 3},
        }

    def plot_freq_qstats(self):
        aof = ao.AOQualityFrequencyStat(self.qs_file)
        stats = self.get_available_stats()

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 8), sharex=True)
        axes = axes.ravel()
        for k, stat in enumerate(stats.keys()):
            for p in self.pol_idx():
                axes[k].plot(
                    aof.freqs / 1e6,
                    aof.get_stat(stat)[:, p] * stats[stat]["xfactor"],
                    label=self.pols[p],
                )
            axes[k].set_ylabel(f"{stats[stat]['label']} {stats[stat]['units']}")
            axes[k].set_title(stat)
            # axes[k].set_yscale("log")
            axes[k].legend()
        for ax in axes[-3:]:
            ax.set_xlabel("Frequency (MHz)")
        axes[-1].remove()
        fig.tight_layout()
        plt.savefig("aoq_frequency_stats.pdf", bbox_inches="tight")
        plt.close(fig)
        return

    def plot_freq_qstats_with_outliers(
        self,
        n_mads=4,
    ):
        aof = ao.AOQualityFrequencyStat(self.qs_file)
        stats = self.get_available_stats()

        # Dictionary to store bad frequencies and fit info
        bad_frequencies = {}

        # Create figure for original data + fit
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), sharex=True)
        axes = axes.ravel()

        # Create separate figure for residuals
        res_fig, res_axes = plt.subplots(
            nrows=3, ncols=2, figsize=(12, 12), sharex=True
        )
        res_axes = res_axes.ravel()

        for k, stat in enumerate(stats.keys()):
            bad_frequencies[stat] = {}
            freqs = aof.freqs / 1e6  # Frequency in MHz

            for p in self.pol_idx():
                # Get the data
                ydata = aof.get_stat(stat)[:, p] * stats[stat]["xfactor"]

                # Fit polynomial
                poly_order = stats[stat]["poly"]
                coeffs = np.polyfit(freqs, ydata, poly_order)
                poly = np.poly1d(coeffs)
                yfit = poly(freqs)
                residuals = ydata - yfit

                # Find outliers in residuals

                median = np.nanmedian(residuals)
                mad = 1.4826 * np.nanmedian(np.abs(residuals - median))  # Scaled MAD
                threshold = n_mads * mad

                res_mean = np.nanmean(residuals)
                res_std = np.nanstd(residuals)
                # threshold = 3 * res_std
                outlier_mask = np.abs(residuals) > threshold
                # outlier_mask = np.abs(residuals - median) > threshold
                outliers = np.where(outlier_mask)[0]
                # outliers = np.where(outlier_mask)[0]

                # Store results
                bad_frequencies[stat][self.pols[p]] = {
                    "frequencies": freqs[outliers],
                    "values": ydata[outliers],
                    "residuals": residuals[outliers],
                    "fit_coeffs": coeffs,
                    "res_mean": res_mean,
                    "res_std": res_std,
                }

                # --- Plot 1: Original data with fit ---
                # Plot original data as solid line
                axes[k].plot(freqs, ydata, label=f"{self.pols[p]} data", alpha=0.7)

                # Plot fit as dashed line
                axes[k].plot(freqs, yfit, "--", label=f"{self.pols[p]} fit", alpha=0.9)

                # Mark outliers with red x's
                axes[k].scatter(
                    freqs[outliers],
                    ydata[outliers],
                    color="red",
                    marker="x",
                    s=40,
                    label=f"{self.pols[p]} outliers",
                    zorder=10,
                )

                # --- Plot 2: Residuals ---
                # Plot residuals as solid line
                res_axes[k].plot(
                    freqs, residuals, label=f"{self.pols[p]} residuals", alpha=0.7
                )

                # Mark outliers in residuals
                res_axes[k].scatter(
                    freqs[outliers],
                    residuals[outliers],
                    color="red",
                    marker="x",
                    s=40,
                    label=f"{self.pols[p]} outliers",
                    zorder=10,
                )

                # Add 3sigma thresholds
                res_axes[k].axhline(threshold, color="gray", linestyle=":", alpha=0.5)
                res_axes[k].axhline(-threshold, color="gray", linestyle=":", alpha=0.5)
                res_axes[k].axhline(0, color="black", linestyle="--", alpha=0.5)

            # Format main plot
            axes[k].set_ylabel(f"{stats[stat]['label']} {stats[stat]['units']}")
            axes[k].set_title(f"{stat} (data + fit)")
            axes[k].legend(fontsize="small")

            # Format residual plot
            res_axes[k].set_ylabel(f"Residuals {stats[stat]['units']}")
            res_axes[k].set_title(f"{stat} residuals")
            res_axes[k].legend(fontsize="small")

        # Clean up and save original+fit plot
        for ax in axes[-3:]:
            ax.set_xlabel("Frequency (MHz)")
        if len(axes) % 2 != 0:  # Remove last axis if odd number
            axes[-1].remove()
        fig.tight_layout()
        fig.savefig("aoq_frequency_stats_with_fit.pdf", bbox_inches="tight")
        plt.close(fig)

        # Clean up and save residuals plot
        for ax in res_axes[-3:]:
            ax.set_xlabel("Frequency (MHz)")
        if len(res_axes) % 2 != 0:  # Remove last axis if odd number
            res_axes[-1].remove()
        res_fig.tight_layout()
        res_fig.savefig("aoq_frequency_residuals.pdf", bbox_inches="tight")
        plt.close(res_fig)

        #  Write bad frequencies to JSON file
        output_file = "aoq_frequency_stats_outliers.pkl"

        with open(output_file, "wb") as f:
            pickle.dump(bad_frequencies, f)

        print(f"Bad frequencies written to {output_file} (JSON format)")

        return bad_frequencies

    def plot_time_qstats(self):
        aot = ao.AOQualityTimeStat(self.qs_file)
        stats = self.get_available_stats()

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(16, 8), sharex=True)
        axes = axes.ravel()

        time = np.unique(aot.time)
        ntimesteps = time.shape[0]
        print(ntimesteps)
        duration_sec = time[-1] - time[0]
        divisor = 3600 if duration_sec > 3600 else 60
        time = np.linspace(0, duration_sec, ntimesteps) / divisor

        for k, stat in enumerate(stats.keys()):
            for p in self.pol_idx():
                data = np.nanmean(
                    aot.get_stat(stat)[:, p].reshape(-1, ntimesteps), axis=0
                )
                data *= stats[stat]["xfactor"]
                axes[k].plot(time, data, label=self.pols[p])
            axes[k].set_ylabel(f"{stats[stat]['label']} ({stats[stat]['units']})")
            axes[k].legend()
        for ax in axes[-3:]:
            ax.set_xlabel("Time (hr)" if divisor == 3600 else "Time (min)")

        axes[-1].remove()
        fig.tight_layout()
        plt.savefig("aoq_time_stats.pdf", bbox_inches="tight")
        plt.close(fig)
        return

    def plot_baseline_qstats(self):
        aob = ao.AOQualityBaselineStat(self.qs_file)
        stats = self.get_available_stats()

        for k, stat in enumerate(stats.keys()):
            for p in self.pol_idx():
                fb1 = aob.plot_baseline_stats(stat, log=True, name="", pol=p)
                fb2 = aob.plot_antennae_stats(stat, log=True, name="", pol=p)
                fb3 = aob.plot_baseline_length_stats(stat, log=True, name="", pol=p)

                fb1.savefig(
                    f"aoq_{self.name}_baseline_{stat}_{self.pols[p]}.pdf",
                    bbox_inches="tight",
                )
                fb2.savefig(
                    f"aoq_{self.name}_antenna_{stat}_{self.pols[p]}.pdf",
                    bbox_inches="tight",
                )
                fb3.savefig(
                    f"aoq_{self.name}_baseline_length_{stat}_{self.pols[p]}.pdf",
                    bbox_inches="tight",
                )
                for ff in [fb1, fb2, fb3]:
                    plt.close(ff)
        return


@click.group()
def main():
    """plot data quality statistics..."""


@main.command("plot_aoq")
@click.argument("qstats_file")
@click.option("--name", help="output name", type=str, default="")
def plot_aoqstats(qstats_file, name=""):
    """Plot aoquality statistics"""

    qq = QualityStats(qstats_file, name)

    qq.plot_time_qstats()
    qq.plot_freq_qstats()
    bf = qq.plot_freq_qstats_with_outliers()
    print(bf)
    qq.plot_baseline_qstats()


if __name__ == "__main__":
    main()

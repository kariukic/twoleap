import os
import h5py
import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import comb

from argparse import ArgumentParser
parser = ArgumentParser(description="Fit a bernstein polynomial to solutions")

parser.add_argument(
    "-s",
    "--sols_file",
    help="DP3 h5 solutions file",
    required=True,
)

parser.add_argument(
    "-d",
    "--degree",
    type=int,
    help="Bpol degree",
    required=False
)


# Define the generalized Bernstein polynomial
def bernstein_polynomial(x, *coeffs):
    """
    Generalized Bernstein polynomial.
    Parameters:
        x (array-like): Input values (normalized to [0, 1]).
        coeffs (array-like): Coefficients for the Bernstein basis.
    Returns:
        y (array-like): Evaluated polynomial values.
    """
    n = len(coeffs) - 1  # Degree of the polynomial
    y = np.zeros_like(x)
    for k, c in enumerate(coeffs):
        y += c * comb(n, k) * (x**k) * ((1 - x)**(n - k))
    return y


def fit(data, degree: int = 3):

    x_data = np.linspace(0, 1, len(data))  # Normalized x-values

    # Fit the generalized Bernstein polynomial
    initial_guess = np.ones(degree + 1)  # Initial guess for coefficients
    popt, pcov = curve_fit(bernstein_polynomial, x_data, data, p0=initial_guess)

    y_fit = bernstein_polynomial(x_data, *popt)

    return y_fit


def write_fitted_gains(h5file, degree=3):
    output_path = h5file.replace('.h5', f'_degree{degree}_bpol.h5')
    os.system(f"cp {h5file} {output_path}")

    with h5py.File(output_path, 'r+') as h5_out:

        cc = h5_out['sol000/amplitude000']["val"][:]

        for d in range(cc.shape[3]):  # per direction
            for t in range(cc.shape[0]):  # per tstep
                for a in range(cc.shape[2]):  # per antenna
                    for p in range(cc.shape[4]):
                        # Overwrite the data with fitted solutions
                        raw_gains = h5_out['sol000/amplitude000']["val"][t, :, a, d, p]
                        if np.isnan(raw_gains).any():
                            raw_gains = np.nan_to_num(raw_gains)
                        h5_out['sol000/amplitude000']["val"][t, :, a, d, p] = fit(raw_gains, degree=degree)
    return


if __name__ == "__main__":
    args = parser.parse_args()
    write_fitted_gains(args.sols_file, degree=args.degree)

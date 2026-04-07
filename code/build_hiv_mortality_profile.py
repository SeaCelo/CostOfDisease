"""
Rebuild the frozen age-specific HIV mortality profile from the raw GBD
South Africa HIV/AIDS results CSV.

The checked-in profile lives at the path given by
``get_pop_data.HIV_MORTALITY_PROFILE_PATH`` and is the runtime source of
truth for the paper.  This script only needs to be rerun if the raw GBD
source file changes or the interpolation method is revised.
"""

import csv
import os

import numpy as np
from scipy.interpolate import PchipInterpolator

import get_pop_data

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
GBD_HIV_RATE_DATA_PATH = os.path.abspath(
    os.path.join(
        CUR_DIR,
        "..",
        "source",
        "JDE",
        "hiv-data",
        "IHME-GBD_2023_DATA-ddf37f70-1",
        "IHME-GBD_2023_DATA-ddf37f70-1.csv",
    )
)
GBD_HIV_RATE_GROUPS = [
    ("<1 year", [0], 0.0),
    ("12-23 months", [1], 1.0),
    ("2-4 years", [2, 3, 4], 3.0),
    ("5-9 years", [5, 6, 7, 8, 9], 7.0),
    ("10-14 years", [10, 11, 12, 13, 14], 12.0),
    ("15-19 years", [15, 16, 17, 18, 19], 17.0),
    ("20-24 years", [20, 21, 22, 23, 24], 22.0),
    ("25-29 years", [25, 26, 27, 28, 29], 27.0),
    ("30-34 years", [30, 31, 32, 33, 34], 32.0),
    ("35-39 years", [35, 36, 37, 38, 39], 37.0),
    ("40-44 years", [40, 41, 42, 43, 44], 42.0),
    ("45-49 years", [45, 46, 47, 48, 49], 47.0),
    ("50-54 years", [50, 51, 52, 53, 54], 52.0),
    ("55-59 years", [55, 56, 57, 58, 59], 57.0),
    ("60-64 years", [60, 61, 62, 63, 64], 62.0),
    ("65-69 years", [65, 66, 67, 68, 69], 67.0),
    ("70-74 years", [70, 71, 72, 73, 74], 72.0),
    ("75-79 years", [75, 76, 77, 78, 79], 77.0),
    ("80-84 years", [80, 81, 82, 83, 84], 82.0),
    ("85-89 years", [85, 86, 87, 88, 89], 87.0),
    ("90-94 years", [90, 91, 92, 93, 94], 92.0),
    ("95+ years", [95, 96, 97, 98, 99], 97.0),
]


def build_gbd_hiv_mortality_profile(
    csv_path=GBD_HIV_RATE_DATA_PATH,
    location_name="South Africa",
    sex_name="Both",
    year=2023,
    num_ages=100,
    min_rate=1e-12,
):
    """
    Build a smooth exact-age HIV/AIDS mortality profile from GBD.

    The GBD results CSV includes overlapping summary ages. This loader uses
    only the finest non-overlapping groups defined in ``GBD_HIV_RATE_GROUPS``
    and returns an exact-age mortality profile aligned to the model ages.

    Args:
        csv_path (str): path to the GBD results CSV
        location_name (str): location name in the CSV
        sex_name (str): sex name in the CSV
        year (int): year in the CSV
        num_ages (int): number of model ages
        min_rate (float): lower bound used on the log-rate scale

    Returns:
        np.array: smooth exact-age HIV mortality profile

    """
    covered_ages = np.zeros(num_ages, dtype=bool)
    required_age_labels = []
    anchor_ages = []
    for age_label, model_ages, anchor_age in GBD_HIV_RATE_GROUPS:
        required_age_labels.append(age_label)
        anchor_ages.append(anchor_age)
        for age in model_ages:
            if age < 0 or age >= num_ages:
                raise ValueError(
                    f"Invalid model age {age} in GBD group {age_label}."
                )
            if covered_ages[age]:
                raise ValueError(
                    f"GBD HIV age groups overlap at model age {age}."
                )
            covered_ages[age] = True

    if not covered_ages.all():
        missing_ages = np.where(~covered_ages)[0]
        raise ValueError(
            "GBD HIV age groups do not cover model ages "
            f"{missing_ages}."
        )

    gbd_group_rates = {}
    with open(csv_path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            if (
                row["location_name"] == location_name
                and row["sex_name"] == sex_name
                and int(row["year"]) == year
                and row["measure_name"] == "Deaths"
                and row["cause_name"] == "HIV/AIDS"
                and row["metric_name"] == "Rate"
                and row["age_name"] in required_age_labels
            ):
                gbd_group_rates[row["age_name"]] = (
                    float(row["val"]) / 100_000.0
                )

    missing_labels = sorted(set(required_age_labels) - gbd_group_rates.keys())
    if missing_labels:
        raise ValueError(
            "Missing GBD HIV rate rows for age groups "
            f"{missing_labels}."
        )

    anchor_rates = np.array(
        [
            max(gbd_group_rates[age_label], min_rate)
            for age_label in required_age_labels
        ],
        dtype=float,
    )
    log_rate_interpolator = PchipInterpolator(
        np.array(anchor_ages, dtype=float),
        np.log(anchor_rates),
        extrapolate=False,
    )

    exact_age_hiv_rates = np.zeros(num_ages)
    max_anchor_age = int(anchor_ages[-1])
    interior_ages = np.arange(max_anchor_age + 1, dtype=float)
    exact_age_hiv_rates[: max_anchor_age + 1] = np.exp(
        log_rate_interpolator(interior_ages)
    )
    exact_age_hiv_rates[max_anchor_age + 1 :] = anchor_rates[-1]

    return np.maximum(exact_age_hiv_rates, 0.0)


def main():
    hiv_mortality_profile = build_gbd_hiv_mortality_profile()
    output_path = get_pop_data.HIV_MORTALITY_PROFILE_PATH
    np.savetxt(
        output_path,
        hiv_mortality_profile,
        delimiter=",",
    )
    print(f"Saved HIV mortality profile to {output_path}")


if __name__ == "__main__":
    main()

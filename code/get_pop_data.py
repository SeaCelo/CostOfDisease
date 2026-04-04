import numpy as np
from scipy.optimize import brentq, minimize
from ogcore import demographics
from ogcore.parameters import Specifications
import os

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
DEMOG_PATH = os.path.join(CUR_DIR, "demographic_data")


def baseline_pop(p, un_country_code="710", download=False):
    """
    Returns objects for the baseline population

    Args:
        p (Specifications): baseline parameters
        un_country_code (str): UN country code
        download (bool): whether to download the data or use that
            checked into repo

    Returns:
        dict: population objects
        fert_rates (np.array): fertility rates
        mort_rates (np.array): mortality rates
        infmort_rates (np.array): infant mortality rates
        imm_rates (np.array): immigration rates
    """
    if download:
        # get initial population objects
        pop_dist, pre_pop_dist = demographics.get_pop(
            p.E,
            p.S,
            0,
            99,
            country_id=un_country_code,
            start_year=p.start_year,
            end_year=p.start_year + 1,
            download_path=DEMOG_PATH,
        )
        fert_rates = demographics.get_fert(
            p.E + p.S,
            0,
            99,
            country_id=un_country_code,
            start_year=p.start_year,
            end_year=p.start_year + 1,
            graph=False,
            download_path=DEMOG_PATH,
        )
        mort_rates, infmort_rates = demographics.get_mort(
            p.E + p.S,
            0,
            99,
            country_id=un_country_code,
            start_year=p.start_year,
            end_year=p.start_year + 1,
            graph=False,
            download_path=DEMOG_PATH,
        )
        imm_rates = demographics.get_imm_rates(
            p.E + p.S,
            0,
            99,
            country_id=un_country_code,
            fert_rates=fert_rates,
            mort_rates=mort_rates,
            infmort_rates=infmort_rates,
            pop_dist=pop_dist,
            start_year=p.start_year,
            end_year=p.start_year + 1,
            graph=False,
            download_path=DEMOG_PATH,
        )
    else:
        pop_dist = np.loadtxt(
            os.path.join(DEMOG_PATH, "population_distribution.csv"),
            delimiter=",",
        )
        pre_pop_dist = np.loadtxt(
            os.path.join(DEMOG_PATH, "pre_period_population_distribution.csv"),
            delimiter=",",
        )
        fert_rates = np.loadtxt(
            os.path.join(DEMOG_PATH, "fert_rates.csv"), delimiter=","
        )
        mort_rates = np.loadtxt(
            os.path.join(DEMOG_PATH, "mort_rates.csv"), delimiter=","
        )
        infmort_rates = np.loadtxt(
            os.path.join(DEMOG_PATH, "infmort_rates.csv"), delimiter=","
        )
        imm_rates = np.loadtxt(
            os.path.join(DEMOG_PATH, "immigration_rates.csv"), delimiter=","
        )

    deaths = total_deaths(
        pop_dist,
        fert_rates,
        mort_rates,
        infmort_rates,
        imm_rates,
        num_years=200,
    )

    pop_dict = demographics.get_pop_objs(
        p.E,
        p.S,
        p.T,
        0,
        99,
        country_id=un_country_code,
        fert_rates=fert_rates,
        mort_rates=mort_rates,
        infmort_rates=infmort_rates,
        imm_rates=imm_rates,
        infer_pop=True,
        pop_dist=pop_dist[:1, :],
        pre_pop_dist=pre_pop_dist,
        initial_data_year=p.start_year,
        final_data_year=p.start_year + 1,
        GraphDiag=False,
    )

    return (
        pop_dict,
        pop_dist,
        pre_pop_dist,
        fert_rates,
        mort_rates,
        infmort_rates,
        imm_rates,
        deaths,
    )


def excess_death_dist(
    scale_factor, pop_dist, mort_rates, excess_deaths=202_693
):
    """
    Find the scale factor to apply to the mortality rates to achieve the
    desired number of excess deaths.

    Args:
        scale_factor (float): factor to apply to the mortality rates
        pop_dist (np.array): population distribution
        mort_rates (np.array): mortality rates
        excess_deaths (int): number of excess deaths to achieve

    Returns:
        float: distance between predicted excess deaths and desired

    """
    current_deaths = np.sum(pop_dist * mort_rates)
    alt_mort_rates = np.minimum(mort_rates * (1 + scale_factor), 1.0)
    new_deaths = np.sum(pop_dist * alt_mort_rates)
    predicted_execess_deaths = new_deaths - current_deaths
    dist = (predicted_execess_deaths - excess_deaths) ** 2

    return dist


def excess_death_gap(
    scale_factor, pop_dist, mort_rates, excess_deaths=202_693
):
    """
    Compute the signed excess-deaths gap for a mortality scale factor.

    Args:
        scale_factor (float): factor to apply to the mortality rates
        pop_dist (np.array): population distribution in the target ages
        mort_rates (np.array): mortality rates in the target ages
        excess_deaths (float): target number of excess deaths

    Returns:
        float: predicted excess deaths minus target excess deaths

    """
    current_deaths = np.sum(pop_dist * mort_rates)
    alt_mort_rates = np.minimum(mort_rates * (1 + scale_factor), 1.0)
    new_deaths = np.sum(pop_dist * alt_mort_rates)
    return new_deaths - current_deaths - excess_deaths


def validate_age_bins(age_bins, age_shares, num_ages=100):
    """
    Validate that bins cover each age exactly once and shares sum to one.

    Args:
        age_bins (list): age bins as ``[(start, end), ...]``
        age_shares (np.array): target excess-death shares by bin
        num_ages (int): number of model ages covered by bins

    """
    if len(age_bins) != len(age_shares):
        raise ValueError("age_bins and age_shares must have the same length.")

    if not np.isclose(age_shares.sum(), 1.0):
        raise ValueError(
            "age_shares must sum to one, but sum to "
            f"{age_shares.sum()}."
        )

    covered_ages = np.zeros(num_ages, dtype=bool)
    for age_start, age_end in age_bins:
        if age_start < 0 or age_end > num_ages or age_start >= age_end:
            raise ValueError(
                f"Invalid age bin ({age_start}, {age_end}) for "
                f"{num_ages} ages."
            )
        if covered_ages[age_start:age_end].any():
            raise ValueError(
                f"Age bins overlap in ({age_start}, {age_end})."
            )
        covered_ages[age_start:age_end] = True

    if not covered_ages.all():
        missing_ages = np.where(~covered_ages)[0]
        raise ValueError(f"Age bins do not cover ages {missing_ages}.")


def solve_scale_factor_brentq(pop_dist, mort_rates, excess_deaths):
    """
    Solve for one mortality multiplier using a bracketed root finder.

    Args:
        pop_dist (np.array): population distribution in the target ages
        mort_rates (np.array): baseline mortality rates in the target ages
        excess_deaths (float): target excess deaths in the target ages

    Returns:
        float: multiplicative mortality increment ``phi``

    """
    if excess_deaths == 0:
        return 0.0

    max_excess_deaths = np.sum(pop_dist * (1.0 - mort_rates))
    if excess_deaths > max_excess_deaths + 1e-8:
        raise ValueError(
            "Target excess deaths exceed the bin capacity under "
            f"the mortality cap: target={excess_deaths}, "
            f"capacity={max_excess_deaths}."
        )

    lower_bound = 0.0
    upper_bound = 0.1
    while (
        excess_death_gap(
            upper_bound, pop_dist, mort_rates, excess_deaths
        )
        < 0
    ):
        upper_bound *= 2
        if upper_bound > 1e6:
            raise RuntimeError(
                "Could not bracket the mortality scale factor root."
            )

    return brentq(
        excess_death_gap,
        lower_bound,
        upper_bound,
        args=(pop_dist, mort_rates, excess_deaths),
    )


def calibrate_age_bin_scale_factors(
    pop_dist,
    mort_rates,
    excess_deaths,
    age_bins,
    age_shares,
):
    """
    Solve one mortality scale factor per age bin.

    Args:
        pop_dist (np.array): population distribution
        mort_rates (np.array): full-intensity baseline mortality row
        excess_deaths (float): total target annual excess deaths
        age_bins (list): age bins as ``[(start, end), ...]``
        age_shares (np.array): target shares by bin

    Returns:
        np.array: scale factors by age bin

    """
    validate_age_bins(age_bins, np.asarray(age_shares), num_ages=len(mort_rates))

    scale_factors = np.zeros(len(age_bins))
    for bin_idx, (age_start, age_end) in enumerate(age_bins):
        target_bin_excess = age_shares[bin_idx] * excess_deaths
        scale_factors[bin_idx] = solve_scale_factor_brentq(
            pop_dist[age_start:age_end],
            mort_rates[age_start:age_end],
            target_bin_excess,
        )

    return scale_factors


def phase_in_age_bin_mortality_rates(
    mort_rates,
    age_bins,
    scale_factors,
    phase_in_years,
):
    """
    Build a phase-in path for age-bin-specific mortality shocks.

    Args:
        mort_rates (np.array): baseline mortality path
        age_bins (list): age bins as ``[(start, end), ...]``
        scale_factors (np.array): mortality scale factors by bin
        phase_in_years (int): number of periods to phase in the shock

    Returns:
        np.array: mortality path with bin-specific shocks applied

    """
    alt_mort_rates = mort_rates.copy()
    for i in range(phase_in_years):
        phase_in_weight = (i + 1) / phase_in_years
        for bin_idx, (age_start, age_end) in enumerate(age_bins):
            alt_mort_rates[i, age_start:age_end] = np.minimum(
                mort_rates[i, age_start:age_end]
                * (1 + scale_factors[bin_idx] * phase_in_weight),
                1.0,
            )

    return alt_mort_rates


def final_year_age_bin_excess_deaths(
    pop_dist,
    fert_rates,
    mort_rates,
    infmort_rates,
    imm_rates,
    baseline_final_year_deaths,
    age_bins,
    scale_factors,
):
    """
    Compute realized final phase-in year excess deaths by age bin.

    Args:
        pop_dist (np.array): initial population distribution
        fert_rates (np.array): baseline fertility-rate path
        mort_rates (np.array): baseline mortality-rate path
        infmort_rates (np.array): baseline infant-mortality-rate path
        imm_rates (np.array): baseline immigration-rate path
        baseline_final_year_deaths (np.array): baseline final-year deaths by age
        age_bins (list): age bins as ``[(start, end), ...]``
        scale_factors (np.array): mortality scale factors by bin

    Returns:
        np.array: realized final-year excess deaths by bin

    """
    alt_mort_rates = phase_in_age_bin_mortality_rates(
        mort_rates,
        age_bins,
        scale_factors,
        mort_rates.shape[0],
    )
    alt_deaths = total_deaths(
        pop_dist,
        fert_rates,
        alt_mort_rates,
        infmort_rates,
        imm_rates,
        num_years=mort_rates.shape[0],
    )
    final_year_excess_deaths = (
        alt_deaths[mort_rates.shape[0] - 1, :] - baseline_final_year_deaths
    )

    return np.array(
        [
            final_year_excess_deaths[age_start:age_end].sum()
            for age_start, age_end in age_bins
        ]
    )


def calibrate_dynamic_age_bin_scale_factors(
    pop_dist,
    fert_rates,
    mort_rates,
    infmort_rates,
    imm_rates,
    baseline_final_year_deaths,
    excess_deaths,
    age_bins,
    age_shares,
    tolerance=1.0,
    max_iter=100,
    update_weight=0.5,
):
    """
    Iteratively match realized final-year excess deaths by age bin.

    The static one-period calibration is used as the initial guess. The update
    loop then rescales each bin's ``phi_b`` based on realized final-year excess
    deaths after the full phase-in path and endogenous demographic transitions.

    Args:
        pop_dist (np.array): initial population distribution
        fert_rates (np.array): baseline fertility-rate path
        mort_rates (np.array): baseline mortality-rate path
        infmort_rates (np.array): baseline infant-mortality-rate path
        imm_rates (np.array): baseline immigration-rate path
        baseline_final_year_deaths (np.array): baseline final-year deaths by age
        excess_deaths (float): total target annual excess deaths
        age_bins (list): age bins as ``[(start, end), ...]``
        age_shares (np.array): target excess-death shares by bin
        tolerance (float): max absolute bin gap in deaths
        max_iter (int): maximum calibration iterations
        update_weight (float): damping for multiplicative updates

    Returns:
        np.array: mortality scale factors by age bin

    """
    age_shares = np.asarray(age_shares)
    validate_age_bins(age_bins, age_shares, num_ages=mort_rates.shape[1])

    target_bin_excess_deaths = age_shares * excess_deaths
    scale_factors = calibrate_age_bin_scale_factors(
        pop_dist[0, :],
        mort_rates[-1, :],
        excess_deaths,
        age_bins,
        age_shares,
    )

    for _ in range(max_iter):
        realized_bin_excess_deaths = final_year_age_bin_excess_deaths(
            pop_dist,
            fert_rates,
            mort_rates,
            infmort_rates,
            imm_rates,
            baseline_final_year_deaths,
            age_bins,
            scale_factors,
        )
        bin_gaps = (
            realized_bin_excess_deaths - target_bin_excess_deaths
        )
        if np.max(np.abs(bin_gaps)) <= tolerance:
            return scale_factors

        safe_realized = np.where(
            np.abs(realized_bin_excess_deaths) > 1e-12,
            realized_bin_excess_deaths,
            1e-12,
        )
        update_ratio = target_bin_excess_deaths / safe_realized
        scale_factors = np.maximum(
            scale_factors
            * (1 + update_weight * (update_ratio - 1)),
            0.0,
        )

    raise RuntimeError(
        "Dynamic age-bin mortality calibration did not converge. "
        f"Final absolute bin gaps: {np.abs(bin_gaps)}."
    )


def extrapolate_demographics(rates, num_years):
    """
    Extend demographic inputs by repeating the final observed row.

    Args:
        rates (np.array): demographic rate path with year in axis 0
        num_years (int): number of years to return

    Returns:
        np.array: rate path with length num_years along axis 0

    """
    if rates.shape[0] >= num_years:
        return rates[:num_years].copy()

    extra_rows = np.repeat(rates[-1:], num_years - rates.shape[0], axis=0)
    return np.concatenate((rates, extra_rows), axis=0)


def disease_pop(
    p,
    pop_dist,
    pre_pop_dist,
    fert_rates,
    mort_rates,
    infmort_rates,
    imm_rates,
    un_country_code="710",
    excess_deaths=202_693,
    age_bins=None,
    age_shares=None,
    phase_in_years=5,
):
    """
    Modify mortality and then get new population objects
    Estimated lives saved per year in South Africa: 202,693
    Source: https://www.cgdev.org/blog/how-many-lives-does-us-foreign-aid-save
    If age bins and shares are provided, apply bin-specific mortality
    multipliers calibrated to those shares while preserving the baseline
    within-bin age profile. Otherwise, fall back to one proportional
    all-age mortality multiplier.

    Args:
        p (Specifications): baseline parameters
        fert_rates (np.array): fertility rates
        mort_rates (np.array): mortality rates
        infmort_rates (np.array): infant mortality rates
        imm_rates (np.array): immigration rates
        excess_deaths (int): number of excess deaths to achieve
        age_bins (list): optional age bins for bin-specific mortality shocks
        age_shares (np.array): optional excess-death shares by age bin
        phase_in_years (int): number of years to phase in mortality changes

    Returns:
        dict: population objects

    """

    # Preserve the baseline demographic path as 2025, 2026, 2026, ...
    num_years = phase_in_years
    fert_rates = extrapolate_demographics(fert_rates, num_years)
    mort_rates = extrapolate_demographics(mort_rates, num_years)
    imm_rates = extrapolate_demographics(imm_rates, num_years)
    alt_infmort_rates = extrapolate_demographics(infmort_rates, num_years)

    # Phase in the mortality changes over the requested number of years.
    alt_mort_rates = mort_rates.copy()
    if age_bins is None and age_shares is None:
        if excess_deaths == 0:
            scale_factor = 0.0
        else:
            # use the scipy minimize function to find the scale factor
            scale_factor_guess = 0.01
            res = minimize(
                excess_death_dist,
                scale_factor_guess,
                args=(pop_dist[0, :], mort_rates[-1, :], excess_deaths),
            )
            scale_factor = res.x[0]

        for i in range(num_years):
            alt_mort_rates[i, :] = np.minimum(
                mort_rates[i, :]
                * (1 + scale_factor * (i + 1) / num_years),
                1.0,
            )
    elif age_bins is not None and age_shares is not None:
        baseline_final_year_deaths = total_deaths(
            pop_dist,
            fert_rates,
            mort_rates,
            alt_infmort_rates,
            imm_rates,
            num_years=num_years,
        )[num_years - 1, :]
        scale_factors = calibrate_dynamic_age_bin_scale_factors(
            pop_dist,
            fert_rates,
            mort_rates,
            alt_infmort_rates,
            imm_rates,
            baseline_final_year_deaths,
            excess_deaths,
            age_bins,
            age_shares,
        )
        alt_mort_rates = phase_in_age_bin_mortality_rates(
            mort_rates,
            age_bins,
            scale_factors,
            num_years,
        )
    else:
        raise ValueError(
            "age_bins and age_shares must either both be provided or "
            "both be None."
        )

    deaths = total_deaths(
        pop_dist,
        fert_rates,
        alt_mort_rates,
        alt_infmort_rates,
        imm_rates,
        num_years=200,
    )

    pop_dict = demographics.get_pop_objs(
        p.E,
        p.S,
        p.T,
        0,
        99,
        country_id=un_country_code,
        fert_rates=fert_rates,
        mort_rates=alt_mort_rates,
        infmort_rates=alt_infmort_rates,
        imm_rates=imm_rates,
        infer_pop=True,
        pop_dist=pop_dist[:1, :],
        pre_pop_dist=pre_pop_dist,
        initial_data_year=p.start_year,
        final_data_year=p.start_year + num_years - 1,
        GraphDiag=False,
    )

    return pop_dict, deaths


def total_deaths(
    pop_dist, fert_rates, mort_rates, infmort_rates, imm_rates, num_years=200
):
    """
    This function computes total deaths each year for num_years.

    Args:
        pop_dist (NumPy array): number of people of each age
        fert_rates (NumPy array): fertility rates at each age
        mort_rates (NumPy array): mortality rates at each age
        infmort_rates (NumPy array): infant mortality rates by year
        imm_rates (NumPy array): immigration rates at each age
        num_years (int): number of years for death forecast

    Returns
        deaths fert_rates (NumPy array): number deaths for each year and age
    """
    # start by looping over years in population objects
    initial_years = mort_rates.shape[0]
    # initialize death array
    deaths = np.zeros((num_years, mort_rates.shape[1]))
    # Loop over years in pop data passed in
    pop_t = pop_dist[0, :]
    for y in range(initial_years):
        deaths[y, :] = pop_t * mort_rates[y, :]
        pop_tp1 = np.zeros_like(pop_t)
        pop_tp1[1:] = (
            pop_t[:-1] * (1 - mort_rates[y, :-1])
            + pop_t[:-1] * imm_rates[y, :-1]
        )
        pop_tp1[0] = (pop_t * fert_rates[y, :]).sum() * (1 - infmort_rates[y])
        pop_t = pop_tp1
    # now loop over all years for which pop data fixed
    for yy in range(initial_years, num_years):
        deaths[yy, :] = pop_t * mort_rates[-1, :]
        pop_tp1 = np.zeros_like(pop_t)
        pop_tp1[1:] = (
            pop_t[:-1] * (1 - mort_rates[-1, :-1])
            + pop_t[:-1] * imm_rates[-1, :-1]
        )
        pop_tp1[0] = (pop_t * fert_rates[-1, :]).sum() * (
            1 - infmort_rates[-1]
        )
        pop_t = pop_tp1

    return deaths

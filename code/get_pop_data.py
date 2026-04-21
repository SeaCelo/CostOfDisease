import numpy as np
from scipy.optimize import brentq, minimize
from ogcore import demographics
from ogcore.parameters import Specifications
import os

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
DEMOG_PATH = os.path.join(CUR_DIR, "demographic_data")
HIV_MORTALITY_PROFILE_PATH = os.path.join(
    DEMOG_PATH, "hiv_mortality_profile_gbd_sa_2023.csv"
)


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


def load_hiv_mortality_profile(
    profile_path=HIV_MORTALITY_PROFILE_PATH,
    num_ages=100,
):
    """
    Load the checked-in age-specific HIV mortality profile used in the paper.

    Args:
        profile_path (str): path to the frozen exact-age profile
        num_ages (int): expected number of model ages

    Returns:
        np.array: age-specific HIV mortality profile

    """
    hiv_mortality_profile = np.loadtxt(profile_path, delimiter=",")
    if hiv_mortality_profile.ndim != 1:
        raise ValueError(
            "The HIV mortality profile must be a one-dimensional age profile."
        )
    if hiv_mortality_profile.shape[0] != num_ages:
        raise ValueError(
            "The HIV mortality profile length does not match the model ages: "
            f"{hiv_mortality_profile.shape[0]} vs {num_ages}."
        )
    if np.any(hiv_mortality_profile < 0):
        raise ValueError("The HIV mortality profile must be nonnegative.")
    if not np.all(np.isfinite(hiv_mortality_profile)):
        raise ValueError("The HIV mortality profile must be finite.")

    return hiv_mortality_profile


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
    hiv_mortality_profile=None,
    phase_in_years=5,
):
    """
    Modify mortality and then get new population objects
    Estimated lives saved per year in South Africa: 202,693
    Source: https://www.cgdev.org/blog/how-many-lives-does-us-foreign-aid-save
    If ``hiv_mortality_profile`` is provided, use the age-specific HIV
    mortality profile as the shape for an additive mortality shock and solve
    for one scalar so year-5 realized excess deaths match the target.
    Otherwise, fall back to one proportional all-age mortality multiplier.

    Args:
        p (Specifications): baseline parameters
        fert_rates (np.array): fertility rates
        mort_rates (np.array): mortality rates
        infmort_rates (np.array): infant mortality rates
        imm_rates (np.array): immigration rates
        excess_deaths (int): number of excess deaths to achieve
        hiv_mortality_profile (np.array or None): optional 1-D age-specific
            HIV mortality profile used to shape the additive mortality shock
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
    if hiv_mortality_profile is not None:
        if (
            hiv_mortality_profile.ndim != 1
            or hiv_mortality_profile.shape[0] != mort_rates.shape[1]
            or np.any(hiv_mortality_profile < 0)
            or not np.all(np.isfinite(hiv_mortality_profile))
        ):
            raise ValueError(
                "hiv_mortality_profile must be a finite, nonnegative, "
                f"1-D array of length {mort_rates.shape[1]}."
            )
        baseline_final_year_total_deaths = total_deaths(
            pop_dist,
            fert_rates,
            mort_rates,
            alt_infmort_rates,
            imm_rates,
            num_years=num_years,
        )[num_years - 1, :].sum()

        def build_hiv_shock_path(shock_scale):
            shocked_mortality = mort_rates.copy()
            for year_idx in range(num_years):
                phase_in_weight = (year_idx + 1) / num_years
                shocked_mortality[year_idx, :] = np.minimum(
                    mort_rates[year_idx, :]
                    + shock_scale * hiv_mortality_profile * phase_in_weight,
                    1.0,
                )

            return shocked_mortality

        def year5_excess_gap(shock_scale):
            shocked_deaths = total_deaths(
                pop_dist,
                fert_rates,
                build_hiv_shock_path(shock_scale),
                alt_infmort_rates,
                imm_rates,
                num_years=num_years,
            )
            return (
                shocked_deaths[num_years - 1, :].sum()
                - baseline_final_year_total_deaths
                - excess_deaths
            )

        if excess_deaths == 0:
            shock_scale = 0.0
        else:
            lower_bound = 0.0
            upper_bound = 1.0
            while year5_excess_gap(upper_bound) < 0:
                upper_bound *= 2
                if upper_bound > 1e6:
                    raise RuntimeError(
                        "Could not bracket the HIV mortality shock scale root."
                    )
            shock_scale = brentq(
                year5_excess_gap,
                lower_bound,
                upper_bound,
            )

        alt_mort_rates = build_hiv_shock_path(shock_scale)
    else:
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
                mort_rates[i, :] * (1 + scale_factor * (i + 1) / num_years),
                1.0,
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

# This script creates the plots and tables used in the paper on the
# costs of excess deaths due to the lack of foreign aid.

# ipmort the necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.ticker as mticker
from ogcore import parameter_plots as pp
from ogcore import output_plots as op
from ogcore.utils import unstationarize_vars

# Use a custom matplotlib style file for plots
plt.style.use("ogcore.OGcorePlots")

# Set constants
BASELINE_YEAR_TO_PLOT = 2040
TIME_SERIES_PLOT_END_YEAR = 2100
NUM_YEARS_NPV = 100
YEAR_RANGE_MIN = 2040
N_YEARS = 20
YEAR_RANGE_MAX = YEAR_RANGE_MIN + 9


def plots(
    base_tpi, base_params, baseline_deaths, reform_dict, forecast, plot_path
):
    """
    Create plots and tables for the cost of disease

    Args:
        base_tpi (dict): baseline TPI output
        base_params (Specifications): baseline parameters
        baseline_deaths (NumPy array): number deaths by year and age in the basline
        reform_dict (dict): dictionary of reform results
        forecast (pd.DataFrame): forecast of macro variables
        plot_path (str): path to save plots

    """
    # Plot mortality rates in the baseline and reform(s)
    tpi_list = [base_tpi] + [
        reform_dict[k]["tpi_vars"] for k in reform_dict.keys()
    ]
    param_list = [base_params] + [
        reform_dict[k]["params"] for k in reform_dict.keys()
    ]
    labels_list = ["Jan. 2025 Baseline"] + [k for k in reform_dict.keys()]
    # Plot mort rates in different scenarios
    years = [BASELINE_YEAR_TO_PLOT]
    p0 = param_list[0]
    age_per = np.linspace(p0.E, p0.E + p0.S, p0.S)
    fig, ax = plt.subplots()
    for y in years:
        t = y - p0.start_year
        for i, p in enumerate(param_list):
            plt.plot(
                age_per[:-1],
                p.rho[t, :-1],
                label=labels_list[i],
            )
    plt.xlabel(r"Age $s$ (model periods)")
    plt.ylabel(r"Mortality Rates $\rho_{s}$")
    plt.legend(loc="upper left")
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels(["{:,.0%}".format(x) for x in ticks_loc])
    fig_path = os.path.join(plot_path, "mortality_rates")
    plt.savefig(fig_path, dpi=300)

    # Print differences in mortality rates at 2025, 2040, 2100
    print("Mortality Rates at 2025, 2040, 2100")
    for year in [2025, 2040, 2100]:
        print(
            year,
            " diff: ",
            reform_dict["Brink, et al. (2025)"]["params"].rho[
                year - base_params.start_year, 50
            ]
            - base_params.rho[year - base_params.start_year, 50],
        )

    # plot survival rates
    p0 = param_list[0]
    age_per = np.linspace(p0.E, p0.E + p0.S, p0.S)
    fig, ax = plt.subplots()
    for y in years:
        t = y - p0.start_year
        for i, p in enumerate(param_list):
            plt.plot(
                age_per,
                np.cumprod(1 - p.rho[t, :]),
                label=labels_list[i] + " " + str(y),
            )
    plt.xlabel(r"Age $s$ (model periods)")
    plt.ylabel(r"Cumulative Survival Rates")
    plt.legend(loc="lower left")
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels(["{:,.0%}".format(x) for x in ticks_loc])
    fig_path = os.path.join(plot_path, "survival_rates")
    plt.savefig(fig_path, dpi=300)

    # Plot differences in productivity profiles
    fig = pp.plot_ability_profiles(
        base_params,
        p2=reform_dict["Gandhi, et al. (2025)"]["params"],
        log_scale=True,
        t=BASELINE_YEAR_TO_PLOT - base_params.start_year,
        include_title=False,
        path=plot_path,
    )

    # Plot differences in chi_n profiles
    fig = pp.plot_chi_n(
        param_list,
        labels=labels_list,
        years_to_plot=[BASELINE_YEAR_TO_PLOT],
        include_title=False,
        path=plot_path,
    )

    # Plot differences in labor supply profiles
    fig = op.tpi_profiles(
        base_tpi,
        base_params,
        reform_dict["Gandhi, et al. (2025)"]["tpi_vars"],
        reform_dict["Gandhi, et al. (2025)"]["params"],
        by_j=False,
        var="n",
        num_years=5,
        start_year=BASELINE_YEAR_TO_PLOT,
        plot_title="Labor Supply by Age",
        path=os.path.join(plot_path, "labor_supply.png"),
    )

    # Plot population distribution in current period and in 25 years
    # under different scenarios
    fig = pp.plot_population(
        base_params,
        years_to_plot=[
            int(base_params.start_year),
            int(base_params.start_year + 25),
        ],
        include_title=False,
        path=None,
    )
    fig.savefig(os.path.join(plot_path, "pop_dist_baseline.png"), dpi=300)
    for k in reform_dict.keys():
        fig2 = pp.plot_population(
            reform_dict[k]["params"],
            years_to_plot=[
                int(base_params.start_year),
                int(base_params.start_year + 25),
            ],
            include_title=False,
            path=None,
        )
        fig2.savefig(os.path.join(plot_path, f"pop_dist_{k}.png"), dpi=300)
    # plot 2050 population distribution on the same picture
    age_vec = np.arange(base_params.E, base_params.S + base_params.E)
    fig, ax = plt.subplots()
    base_pop = base_params.omega[2050 - base_params.start_year, :]
    plt.plot(age_vec, base_pop, label="Baseline pop.")
    for k in reform_dict.keys():
        pop_dist = reform_dict[k]["params"].omega[
            2050 - base_params.start_year, :
        ]
        plt.plot(age_vec, pop_dist, label=k + " pop.")
    plt.xlabel(r"Age $s$")
    plt.ylabel(r"Pop. dist'n $\omega_{s}$")
    plt.legend(loc="lower left")
    fig.savefig(os.path.join(plot_path, "pop_dist_2050.png"), dpi=300)

    # Plot cumulative excess deaths
    # Compute and plot aggregage deaths by year
    death_dict = {
        "Jan. 2025 Baseline": baseline_deaths.sum(axis=1),
    }
    for k in reform_dict.keys():
        death_dict[k] = reform_dict[k]["deaths"].sum(axis=1)
    # Put in dataframe
    death_df = pd.DataFrame.from_dict(death_dict)
    # save as csv
    death_df.to_csv(os.path.join(plot_path, "deaths.csv"), index=False)
    # plot
    num_years_plot = 100
    fig, ax = plt.subplots()
    for k in reform_dict.keys():
        plt.plot(
            np.arange(
                base_params.start_year, base_params.start_year + num_years_plot
            ),
            (death_df[k] - death_df["Jan. 2025 Baseline"])[
                :num_years_plot
            ].cumsum()
            / 1_000_000,
            label=k,
        )
    plt.xlabel("Year")
    plt.ylabel("Cumulative Excess Deaths (millions)")
    plt.legend(loc="upper left")
    fig.savefig(
        os.path.join(plot_path, "cumulative_excess_deaths.png"), dpi=300
    )

    # Create table of level changes in macro variables
    # African GDP over a long time period (or extrapolate from some shorter term forecast)
    # compute percentage changes in macro variables
    inflation_adjust = (
        313.698 * 1.025
    ) / 237.002  # to put the 2015$ into 2025$
    GDP_series = {
        "Baseline Forecast": forecast[:NUM_YEARS_NPV] * inflation_adjust,
        "Pct Changes": {},
        "Levels": {},
        "Diffs": {},
    }
    for k in reform_dict.keys():
        Y_ref = unstationarize_vars(
            "Y",
            reform_dict[k]["tpi_vars"],
            reform_dict[k]["params"],
        )
        Y_base = unstationarize_vars(
            "Y",
            base_tpi,
            base_params,
        )
        GDP_series["Pct Changes"][k] = (
            Y_ref[:NUM_YEARS_NPV] - Y_base[:NUM_YEARS_NPV]
        ) / Y_base[:NUM_YEARS_NPV]
        GDP_series["Levels"][k] = GDP_series["Baseline Forecast"][
            :NUM_YEARS_NPV
        ] * (1 + GDP_series["Pct Changes"][k][:NUM_YEARS_NPV])
        GDP_series["Diffs"][k] = (
            GDP_series["Levels"][k][:NUM_YEARS_NPV]
            - GDP_series["Baseline Forecast"][:NUM_YEARS_NPV]
        )

    # Find avg change in GDP from begin_window to end_window
    results_first_N_years = {}
    for k in reform_dict.keys():
        results_first_N_years[k] = (
            GDP_series["Diffs"][k][:N_YEARS].mean() / 1e9
        )  # convert to billions of dollars
    results_df = pd.DataFrame.from_dict(results_first_N_years, orient="index")
    # rename column
    results_df.columns = [r"$\Delta$ GDP, Billions"]
    results_df.to_latex(
        buf=os.path.join(plot_path, f"mean_gdp_change_{N_YEARS}years.tex"),
        float_format="%.2f",
    )

    # Find NPV of levels of GDP over NUM_YEARS_NPV years
    results_NPV = {"Discount Rate": [r"1\%", r"2\%", r"3\%", r"4\%", r"6\%"]}
    npv_dict = {
        "Discount Rate": [0.01, 0.02, 0.03, 0.04, 0.06],
        "Discount Rate Label": [r"1\%", r"2\%", r"3\%", r"4\%", r"6\%"],
    }
    for k in reform_dict.keys():
        results_NPV[k] = []
        for r in npv_dict["Discount Rate"]:
            results_NPV[k].append(
                (
                    GDP_series["Diffs"][k][:NUM_YEARS_NPV]
                    / (1 + r) ** np.arange(NUM_YEARS_NPV)
                ).sum()
                / 1e9  # convert to billions of dollars
            )
    # Save table to disk
    formatted_table = pd.DataFrame(results_NPV)
    formatted_table.to_latex(
        buf=os.path.join(plot_path, "npv_gdp_table.tex"),
        index=False,
        index_names=False,
        float_format="%.2f",
    )

    # Time series plots
    # log GDP
    fig, ax = plt.subplots()
    idx = TIME_SERIES_PLOT_END_YEAR - base_params.start_year
    years = np.arange(base_params.start_year, TIME_SERIES_PLOT_END_YEAR)
    plt.plot(
        years,
        np.log(GDP_series["Baseline Forecast"][:idx] / 1e9),
        label="Jan. 2025 Baseline",
    )
    for k in reform_dict.keys():
        plt.plot(
            years,
            np.log(GDP_series["Levels"][k][:idx] / 1e9),
            label=k,
        )
    plt.legend()
    plt.xlabel("Year")
    plt.ylabel("GDP, ln(Billions 2025$)")
    plt.savefig(os.path.join(plot_path, "log_GDP_paths.png"), dpi=300)

    # Differences in GDP
    fig, ax = plt.subplots()
    idx = TIME_SERIES_PLOT_END_YEAR - base_params.start_year
    years = np.arange(base_params.start_year, TIME_SERIES_PLOT_END_YEAR)
    for k in reform_dict.keys():
        plt.plot(
            years,
            GDP_series["Diffs"][k][:idx] / 1e9,
            label=k,
        )
    plt.legend()
    plt.xlabel("Year")
    plt.ylabel("Change in GDP (Billions of 2025$)")
    plt.savefig(os.path.join(plot_path, "GDP_diff_path.png"), dpi=300)

    # Plot the capital labor ratio over time
    fig, ax = plt.subplots()
    idx = 20  # TIME_SERIES_PLOT_END_YEAR - base_params.start_year
    years = np.arange(base_params.start_year, base_params.start_year + idx)
    K_base = unstationarize_vars("K", base_tpi, base_params)
    L_base = unstationarize_vars("L", base_tpi, base_params)
    KL_base = K_base[:idx] / L_base[:idx]
    plt.plot(years, KL_base, label="Jan. 2025 Baseline")
    for k in reform_dict.keys():
        K_ref = unstationarize_vars(
            "K",
            reform_dict[k]["tpi_vars"],
            reform_dict[k]["params"],
        )
        L_ref = unstationarize_vars(
            "L",
            reform_dict[k]["tpi_vars"],
            reform_dict[k]["params"],
        )
        KL_ref = K_ref[:idx] / L_ref[:idx]
        plt.plot(years, KL_ref, label=k)
    plt.legend()
    plt.xlabel("Year")
    plt.ylabel("Capital-Labor Ratio")
    plt.savefig(os.path.join(plot_path, "capital_labor_ratio.png"), dpi=300)

    # plot the percentage change in K and L over time
    # Make K solid and L dashed
    fig, ax = plt.subplots()
    idx = 20  # TIME_SERIES_PLOT_END_YEAR - base_params.start_year
    years = np.arange(base_params.start_year, base_params.start_year + idx)
    K_base = unstationarize_vars("K", base_tpi, base_params)
    L_base = unstationarize_vars("L", base_tpi, base_params)
    k = "Gandhi, et al. (2025)"
    K_ref = unstationarize_vars(
        "K",
        reform_dict[k]["tpi_vars"],
        reform_dict[k]["params"],
    )
    L_ref = unstationarize_vars(
        "L",
        reform_dict[k]["tpi_vars"],
        reform_dict[k]["params"],
    )
    plt.plot(
        years,
        (K_ref[:idx] - K_base[:idx]) / K_base[:idx] * 100,
        label=" K",
        color="black",
    )
    plt.plot(
        years,
        (L_ref[:idx] - L_base[:idx]) / L_base[:idx] * 100,
        label=" L",
        linestyle="dashed",
        color="black",
    )
    plt.legend()
    ax.set_xticks(
        np.arange(base_params.start_year, base_params.start_year + idx, 5)
    )
    plt.xlabel("Year")
    plt.ylabel("Percentage Change from Baseline")
    plt.savefig(
        os.path.join(plot_path, "K_L_percentage_change_medium.png"), dpi=300
    )


def decomposition(base_tpi, base_params, decomp_dict, forecast, plot_path):
    """
    This function decomposes the changes in GDP into the contributions from
    changes in mortality rates and changes in productivity.
    """

    # Create table of level changes in macro variables
    # African GDP over a long time period (or extrapolate from some shorter term forecast)
    # compute percentage changes in macro variables
    inflation_adjust = (
        313.698 * 1.025
    ) / 237.002  # to put the 2015$ into 2025$
    GDP_series = {
        "Baseline Forecast": forecast[:NUM_YEARS_NPV] * inflation_adjust,
        "Pct Changes": {},
        "Levels": {},
        "Diffs": {},
    }
    for k in decomp_dict.keys():
        Y_ref = unstationarize_vars(
            "Y",
            decomp_dict[k]["tpi_vars"],
            decomp_dict[k]["params"],
        )
        Y_base = unstationarize_vars(
            "Y",
            base_tpi,
            base_params,
        )
        GDP_series["Pct Changes"][k] = (
            Y_ref[:NUM_YEARS_NPV] - Y_base[:NUM_YEARS_NPV]
        ) / Y_base[:NUM_YEARS_NPV]
        GDP_series["Levels"][k] = GDP_series["Baseline Forecast"][
            :NUM_YEARS_NPV
        ] * (1 + GDP_series["Pct Changes"][k][:NUM_YEARS_NPV])
        GDP_series["Diffs"][k] = (
            GDP_series["Levels"][k][:NUM_YEARS_NPV]
            - GDP_series["Baseline Forecast"][:NUM_YEARS_NPV]
        )

    # Find NPV of levels of GDP over NUM_YEARS_NPV years
    results_NPV = {"Discount Rate": [r"1\%", r"2\%", r"3\%", r"4\%", r"6\%"]}
    npv_dict = {
        "Discount Rate": [0.01, 0.02, 0.03, 0.04, 0.06],
        "Discount Rate Label": [r"1\%", r"2\%", r"3\%", r"4\%", r"6\%"],
    }
    for k in decomp_dict.keys():

        results_NPV[k] = []
        for r in npv_dict["Discount Rate"]:
            results_NPV[k].append(
                (
                    GDP_series["Diffs"][k][:NUM_YEARS_NPV]
                    / (1 + r) ** np.arange(NUM_YEARS_NPV)
                ).sum()
                / 1e9  # convert to billions of dollars
            )
    # Save table to disk
    formatted_table = pd.DataFrame(results_NPV)
    formatted_table.to_latex(
        buf=os.path.join(plot_path, "npv_gdp_decomposition_table.tex"),
        index=False,
        index_names=False,
        float_format="%.2f",
    )

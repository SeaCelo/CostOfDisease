import numpy as np

import get_pop_data


def main():
    hiv_mortality_profile = get_pop_data.build_gbd_hiv_mortality_profile()
    np.savetxt(
        get_pop_data.HIV_MORTALITY_PROFILE_PATH,
        hiv_mortality_profile,
        delimiter=",",
    )
    print(
        "Saved HIV mortality profile to "
        f"{get_pop_data.HIV_MORTALITY_PROFILE_PATH}"
    )


if __name__ == "__main__":
    main()

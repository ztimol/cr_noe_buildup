import numpy as np
import math


# see 8.4.1.1 Advanced topic: longer mixing times from Keeler chapter 8
def calc_theoretical_noe_intensity(
    mixing_time, lambda_one, lambda_two, r_i, r_s, r, sigma
):

    iz_over_iz_zero = ((2 * sigma) / r) * (
        math.exp(-lambda_two * mixing_time) - math.exp(-lambda_one * mixing_time)
    ) + 1

    sz_over_iz_zero = (
        (
            ((r_i - r_s) / r)
            * (
                math.exp(-lambda_one * mixing_time)
                - math.exp(-lambda_two * mixing_time)
            )
        )
        + 1
        - math.exp(-lambda_one * mixing_time)
        + math.exp(-lambda_two * mixing_time)
    )

    theoretical_noe_intensity = iz_over_iz_zero - sz_over_iz_zero

    return theoretical_noe_intensity


def calc_theoretical_noe_intensities(
    mixing_times, lambda_one, lambda_two, r_i, r_s, r, sigma
):

    theoretical_noe_intensities = {}

    for mixing_time in mixing_times:
        theoretical_noe_intensity = calc_theoretical_noe_intensity(
            mixing_time, lambda_one, lambda_two, r_i, r_s, r, sigma
        )
        theoretical_noe_intensities[mixing_time] = theoretical_noe_intensity

    return theoretical_noe_intensities.values()


def calc_corellation(list_1, list_2):

    assert len(list_1) == len(
        list_2
    ), "The length of both lists have to be equal when calculation r sqaured"

    list_1 = np.array(list(list_1))
    list_2 = np.array(list(list_2))

    fit1 = np.polyfit(list_1, list_2, 2)
    v1 = np.polyval(fit1, list_1)
    mse1 = np.mean(np.square(list_2 - v1))

    r = (
        (len(list_1) * sum([x * y for x, y in zip(list_1, list_2)]))
        - (sum(list_1) * sum(list_2))
    ) / (
        math.sqrt(
            ((len(list_1) * sum([x ** 2 for x in list_1])) - (sum(list_1) ** 2))
            * ((len(list_2) * sum([y ** 2 for y in list_2])) - (sum(list_2) ** 2))
        )
    )
    return mse1, r ** 2


def write_out_params(
    out_file,
    lambda_one,
    lambda_two,
    r_i,
    r_s,
    r,
    sigma,
    theoretical_noe_intensities=[],
    mean_squared_error=None,
    r_squared=None,
):
    out_string = (
        "{lambda_one}, {lambda_two}, {r_i}, {r_s}, {r}, {sigma}, {mean_squared_error}, {r_squared}, "
    ).format(
        lambda_one=lambda_one,
        lambda_two=lambda_two,
        r_i=r_i,
        r_s=r_s,
        r=r,
        sigma=sigma,
        mean_squared_error=mean_squared_error,
        r_squared=r_squared,
    )
    out_file.write(out_string)

    for noe in theoretical_noe_intensities:
        out_file.write(str(noe))
        out_file.write(" ")
    out_file.write("\n")


def fit_curve(measured_noe_intensity_per_mixing_time, out_file):

    # lambda_one_min = 9.0
    # lambda_two_min = 4.5
    # r_min = 3.0
    # r_s_min = 7.0
    # r_i_min = 7.0
    # sigma_min = 2.0

    # lambda_one_max = 9.1
    # lambda_two_max = 4.6
    # r_max = 3.1
    # r_s_max = 7.1
    # r_i_max = 7.1
    # sigma_max = 2.1

    lambda_one_min = 5.0
    lambda_two_min = 2.5
    r_min = 1.0
    r_s_min = 3.0
    r_i_min = 3.0
    sigma_min = 1.0

    lambda_one_max = 12.0
    lambda_two_max = 6.0
    r_max = 5.0
    r_s_max = 10.0
    r_i_max = 10.0
    sigma_max = 5.0

    mixing_times = measured_noe_intensity_per_mixing_time.keys()
    measured_noe_intensities = measured_noe_intensity_per_mixing_time.values()

    for lambda_one in np.arange(lambda_one_min, lambda_one_max, 0.1):
        for lambda_two in np.arange(lambda_two_min, lambda_two_max, 0.1):
            for r_i in np.arange(r_i_min, r_i_max, 0.1):
                for r_s in np.arange(r_s_min, r_s_max, 0.1):
                    for r in np.arange(r_min, r_max, 0.1):
                        for sigma in np.arange(sigma_min, sigma_max, 0.1):
                            if (
                                lambda_two >= lambda_one
                                or lambda_one == 0.0
                                or lambda_two == 0.0
                                or r_i == 0.0
                            ):
                                pass
                            else:
                                try:
                                    theoretical_noe_intensities = calc_theoretical_noe_intensities(
                                        mixing_times,
                                        lambda_one,
                                        lambda_two,
                                        r_i,
                                        r_s,
                                        r,
                                        sigma,
                                    )
                                    mean_squared_error, r_squared = calc_corellation(
                                        measured_noe_intensities,
                                        theoretical_noe_intensities,
                                    )
                                    write_out_params(
                                        out_file,
                                        lambda_one,
                                        lambda_two,
                                        r_i,
                                        r_s,
                                        r,
                                        sigma,
                                        list(theoretical_noe_intensities),
                                        mean_squared_error,
                                        r_squared,
                                    )
                                except:
                                    write_out_params(
                                        out_file,
                                        lambda_one,
                                        lambda_two,
                                        r_i,
                                        r_s,
                                        r,
                                        sigma,
                                    )


def main():

    measured_noe_intensities_per_mixing_time = {
        0.03: 0.019,
        0.04: 0.026,
        0.05: 0.036,
        0.06: 0.043,
        0.07: 0.051,
        0.08: 0.057,
        0.09: 0.064,
        0.11: 0.076,
        0.2: 0.085,
    }

    with open("out.log", "w") as out_file:
        fit_curve(measured_noe_intensities_per_mixing_time, out_file)


main()

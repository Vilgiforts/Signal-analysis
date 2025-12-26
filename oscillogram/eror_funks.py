import scipy.stats as sps
import numpy as np
import matplotlib.pyplot as plt
import random


def confidence_interval(
    array: list[float | int] | tuple[float | int],
    confidence_probability: float = 0.95,
    absolute_error: int | float = 0,
    relative_error: int | float = 0,
):
    t_inf = sps.t.ppf(1 / 2 + confidence_probability / 2, df=10**6)
    if len(array) > 1:
        t = sps.t.ppf(1 / 2 + confidence_probability / 2, df=len(array) - 1)
        mean_valeu = sum(array) / len(array)
        S_0 = (
            sum(((valeu - mean_valeu) ** 2 for valeu in array))
            / (len(array) * (len(array) - 1))
        ) ** 0.5
        rez = (
            mean_valeu,
            (
                (t * S_0) ** 2
                + (t_inf * (absolute_error + relative_error * mean_valeu) / 3) ** 2
            )
            ** 0.5,
        )
        return rez
    else:
        return (array[0], t_inf / 3**0.5 * (absolute_error + relative_error * array[0]))


def uneven_measurements(value_array: list | tuple, error_array: list | tuple):
    if len(value_array) != len(error_array):
        raise ValueError(
            "Список значений и список погрешностей должны быть одного размера!"
        )
    sum_w = sum(map(lambda x: 1 / x**2, error_array))
    mean_value = (
        sum((1 / error_array[i] ** 2 * value_array[i] for i in range(len(value_array))))
        / sum_w
    )
    error_value = (len(value_array) / sum_w) ** 0.5
    return (mean_value, error_value)


def error_extrems(
    array_x: list[float | int] | tuple[float | int],
    array_y: list[float | int] | tuple[float | int],
    confidence_probability: float = 0.95,
    absolute_error: int | float = 0,
    relative_error: int | float = 0,
):
    t = sps.t.ppf(1 / 2 + confidence_probability / 2, df=len(array_x) - 1)
    a = np.polyfit(array_x, array_y, 2)
    # print(a)
    f_approx = np.poly1d(a)
    # plt.scatter(array_x, array_y)
    # plt.plot(
    #     np.arange(
    #         min(array_x), max(array_x) + 1, -(min(array_x) - max(array_x) + 1) / 20
    #     ),
    #     f_approx(
    #         np.arange(
    #             min(array_x), max(array_x) + 1, -(min(array_x) - max(array_x) + 1) / 20
    #         )
    #     ),
    #     ls="--",
    # )
    S_0 = (
        (sum((array_y[i] - f_approx(array_x[i])) ** 2 for i in range(len(array_x))))
        / (len(array_x) * (len(array_x) - 1))
    ) ** 0.5 * 1.5

    x = -a[1] / (2 * a[0])
    # plt.errorbar([x], [f_approx(x)], yerr=[float(t * S_0)], fmt="o", c="red")
    # plt.show()
    return (float(x), float(f_approx(x))), float(t * S_0)


if __name__ == "__main__":
    s = 0
    n = 100000
    for i in range(n):
        x = np.arange(-10, 11, 1)
        abs_error_amp, reletiv_error_amp = 1, 0.01
        abs_error, reletiv_error = 2 * abs_error_amp * (
            np.random.random(len(x)) - 0.5
        ), (1 + 2 * (np.random.random(len(x)) - 0.5) * reletiv_error_amp)
        y = x**2 * reletiv_error + abs_error
        a = error_extrems(x, y, confidence_probability=0.95)
        # print(f"\033[{"39" if -a[1] < a[0][1] < a[1] else "31"}m {a}\033[39m")

        if not (a[0][1] - a[1] <= 0 <= a[0][1] + a[1]):
            s += 1
        if i % 100 == 0:
            print(f"{i+1}:{(s/(i + 1) * 100):.2f}%")
    print(f"{(s/n * 100):.2f}%")

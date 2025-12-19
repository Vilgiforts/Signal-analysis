import scipy.stats as sps


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
    sum_w = sum(map(lambda x: 1 / x**2, error_array))
    mean_value = (
        sum((1 / error_array[i] ** 2 * value_array[i] for i in range(len(value_array))))
        / sum_w
    )
    error_value = (len(value_array) / sum_w) ** 0.5
    return (mean_value, error_value)

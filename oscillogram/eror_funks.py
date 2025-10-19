import scipy.stats as sps


def confidence_interval(
    array: list[float | int] | tuple[float | int],
    confidence_probability: float = 0.95,
    absolute_error: int | float = None,
    relative_error: int | float = None,
):
    t = sps.t.ppf(1 / 2 + confidence_probability / 2, df=len(array) - 1)
    t_inf = sps.t.ppf(1 / 2 + confidence_probability / 2, df=10**6)
    mean_valeu = sum((abs(value) for value in array)) / len(array)
    S_0 = (
        sum(((abs(valeu) - mean_valeu) ** 2 for valeu in array))
        / (len(array) * (len(array) - 1))
    ) ** 0.5
    rez = (
        (mean_valeu, ((t * S_0) ** 2 + (t_inf * absolute_error / 3) ** 2) ** 0.5)
        if absolute_error
        else (
            mean_valeu,
            ((t * S_0) ** 2 + (t_inf * relative_error * mean_valeu / 3) ** 2) ** 0.5,
        )
    )
    return rez

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import random
from .eror_funks import *
import traceback
import sys

text_error = lambda x: f"\033[31m {x} \033[39m"
text_info = lambda x: f"\033[34m {x} \033[39m"


class oscillogram:
    """
    Класс для анализа и вицуализации осцилограмм из текстовых файлов.

    Атрибуты:
        1. path(str) - путь до файла
        2. frecquency(float | list[int | float] | tuple[int | float])[Гц] = None - частота канала в герцах,
        может быть как одним числом, так и набором значений в случае различности частот для разных каналов. Ести частота не задана, ей присваеватеся значение 0.
        3. name(str) = None - имя экземпляра, используется при выводе и постраении графика, если не заданно, присваивается значение из path.
        4. chanels(int) - число каналов, определяется по первой строке файла.
        5. time_start(float)[с] - время начало измерения, используется для обнуления временной шкалы.
        6. times(list[list[float]])[с] - все значения времени для каждого канала.
        7. values(list[list[float]]) - все значения измерямой величины для каждого канала.
        8. nulls(list[list[flaot]])[c] - все значения времени в которых сигнал переходит через 0,
        для каждого кнала, определяется как среднее между ближайщими точками с разыми знаками.
        9. new_apm_nulls(list[float]) - в случае сдвига сигнала относительно оси x = 0, вычисляются новые значения,
        а nulls, amplitude и некоторые другие характеристики, пересчитываются относительно этих новых нулей. Вычесляется только в случае обнаружения хотябы 4 'обычных' нулей.
        10. amplitude(list[float]) - значения амплитуд сигналов для каждого канала, определяется во время выполнения метод find_extrems.
        взависимость от того где достигается максимальное значение может быть как положительным так и отрицательным.
        11. extrems(list[list[tuple[float, float]]]) - экстремальные значения сигнала.
        . phase_shift(list[list(float)]) - матрица сдвига фаз n x n, где n - число каналов.
    Методы:
        1. find_nuls
        2. find_extrems
    """

    def __init__(
        self,
        path: str,
        frequency: float | list[int | float] | tuple[int | float] = None,
        name: str = None,
        relative_erorr_x: int | float = 0,
        relative_erorr_y: int | float = 0.001,
        absolute_error_x: int | float = 0,
        absolute_error_y: int | float = 0,
        confidence_probability: float = 0.95,
    ):
        self.path = path  # Инициализируетм путь до файла
        dot_index = path.find(".")
        self.name = (
            name if name else path[:dot_index] if dot_index != -1 else path
        )  # Инициализируем имя
        self.relative_erorr_x = relative_erorr_x
        self.relative_erorr_y = relative_erorr_y
        self.absolute_error_x = absolute_error_x
        self.absolute_error_y = absolute_error_y
        self.confidence_probability = confidence_probability
        self.t_inf = sps.t.ppf(1 / 2 + confidence_probability / 2, df=10**6)

        # Начало чтения файла
        with open(path, "r") as file:  # Открываем файл
            for i, line in enumerate(file):  # Пробегаемся по строкам файла
                a = line.rstrip(
                    "\n"
                ).split()  # Отбрасываем символ перехода на новую строку и разбиваем строку на элементы массива по пробелу
                if i == 0:  # По первой строке файла определяем количество каналов
                    self.channels = len(re.findall(r"\[(.+?)\]", line))
                    self.time_start = []
                    self.times = [[] for i in range(self.channels)]
                    self.values = [[] for i in range(self.channels)]
                for j in range(
                    self.channels
                ):  # Пробегаемся по элементам массива исча нужные нам данные
                    if i == 1:
                        self.time_start.append(
                            float(  # Преобразуем результат в float
                                np.dot(  # Берем сумму произведений соотвественных элементов(Как скалярное произведение numpy для сокращения записи), для перевода времени из часов, минут, секунд, в секунды
                                    (3600, 60, 1),
                                    list(
                                        map(
                                            lambda x: float(
                                                x.replace(",", ".")
                                            ),  # Меняем ',' на '.'
                                            a[1:][2 * j + 1].split(
                                                ":"
                                            ),  # Разбиваем подстроку по ':' на массив
                                        )
                                    ),
                                )
                            )
                        )
                    elif i == 2:
                        self.delta_time = tuple(
                            float(delta_t.replace(",", ".")) for delta_t in a[2:]
                        )
                        print(self.delta_time)
                    elif i > 4:  # Аналогично с предыдущим
                        self.times[j].append(
                            float(
                                np.dot(
                                    (3600, 60, 1),
                                    list(
                                        map(
                                            lambda x: float(x.replace(",", ".")),
                                            a[3 * j + 1].split(":"),
                                        )
                                    ),
                                )
                            )
                            - self.time_start[
                                j
                            ]  # Вычетаем из результата время начала, чтобы перенести временную шкалу в 0
                        )
                        self.values[j].append(float(a[3 * j + 2].replace(",", ".")))
        # Конец четния файла

        self.nulls = self.find_nulls()  # Определяем нули на осцилограмме нули
        for i in range(self.channels):
            if len(self.nulls[i]) > 4:  # Для каждого канала определяем наличие
                self.new_amp_nulls = [
                    (max(self.values[i]) + min(self.values[i])) / 2
                    for i in range(self.channels)
                ]  # Если есть хотя бы 4 нуля, то мы уверены в наличие хоты бы одного положительного и отрицательного пиков гармонического сигнала и можем вычислить новые нули(выровнять сигнал относительно оси x)
                self.nulls = self.find_nulls()
        self.extrems, self.extrems_error = self.find_extrems()
        self.frequency_error = [0 for i in range(self.channels)]
        if frequency == None:
            self.frequency, self.frequency_error = self.find_frequency()
        elif hasattr(frequency, "__iter__") and len(frequency) == self.channels:
            self.frequency = frequency
        elif not (hasattr(frequency, "__iter__")):
            self.frequency = [frequency for i in range(self.channels)]
        else:
            self.frequency = [0 for i in range(self.channels)]
            print(text_error("Что-то пошло не так при передаче частот каналов!"))

        simplitfy = [self.find_simplify_sig(i) for i in range(self.channels)]
        self.simplify_sig = [i[0] for i in simplitfy]
        self.simplify_sig_errors = [i[1] for i in simplitfy]
        self.phase_shift = [
            [0 for j in range(self.channels)] for i in range(self.channels)
        ]  # Задаем матрицу для сдвигов фаз
        self.phase_shift_error = [
            [0 for j in range(self.channels)] for i in range(self.channels)
        ]  # Задаем матрицу для сдвигов фаз
        for i in range(self.channels):
            for j in range(self.channels):
                if i < j:
                    a = self.find_shift(i, j)
                    self.phase_shift[i][j], self.phase_shift[j][i] = a, a
        self.amplitude = self.find_amplitude()

    def find_nulls(self):  # Проверить на работоспособность
        """Метод нахождения нулей сигнала для каждого канала.
        Returns:
            list[list[float]] - список, списков временных координат, в которых сигнал переходит через 0, вычесляется как среднее от двух ближайших точек разного знака.
        """
        try:
            self.new_amp_nulls
            for k in range(self.channels):
                self.values[k] = list(np.array(self.values[k]) - self.new_amp_nulls[k])
        except:
            self.new_amp_nulls = [0 for i in range(self.channels)]
        rez = [[] for i in range(self.channels)]
        rez_1 = [[] for i in range(self.channels)]
        for i in range(self.channels):
            for j in range(len(self.values[i]) - 1):
                value = self.values[i][j]
                if (value > 0 and self.values[i][j + 1] < 0) or (
                    value < 0 and self.values[i][j + 1] > 0
                ):
                    rez[i].append(self.times[i][j])
            max_distance = 0
            for j in range(len(rez[i]) - 1):
                if abs(rez[i][j + 1] - rez[i][j]) > max_distance:
                    max_distance = abs(rez[i][j + 1] - rez[i][j])
            j = 0
            while j < len(rez[i]):
                interim_nulls = []
                k = j + 1
                count = 0
                while k < len(rez[i]) and abs(rez[i][j] - rez[i][k]) < max_distance / 4:
                    interim_nulls.append(rez[i][k])
                    k += 1
                try:
                    rez_1[i].append(sum(interim_nulls) / len(interim_nulls))
                except:
                    rez_1[i].append(rez[i][j])
                j = k
        return rez_1

    def find_extrems(self):
        """Находит все экстремальные точки сигнала, а так же задает атрибут amplitude для каждого канала.

        Returns:
            list[list[tuple[float, float]]]: список в котором для каждого канала находятся списки кортежей с парами значений (время, значение).
        """
        rez = []
        new_rez = []
        error_rez = []
        for i in range(self.channels):
            rez.append([])
            new_rez.append([])
            error_rez.append([])
            flag = False
            max_value = max(abs(max(self.values[i])), abs(min(self.values[i])))
            if len(self.nulls[i]) >= 2:
                singl = 1 if self.values[i][0] > 0 else -1
                actual_max_value = (0, self.values[i][0])
                for j in range(len(self.values[i])):
                    if singl * self.values[i][j] >= 0 and abs(self.values[i][j]) > abs(
                        actual_max_value[1]
                    ):
                        actual_max_value = (j, self.values[i][j])
                    elif singl * self.values[i][j] < 0:
                        singl *= -1
                        if flag and abs((actual_max_value[1]) / (max_value)) >= 0.7:
                            rez[i].append(actual_max_value)
                        actual_max_value = (0, self.values[i][j])
                        flag = True
            if len(self.nulls[i]) < 2:
                if abs(max(self.values[i])) > abs(min(self.values[i])):
                    rez[i].append(
                        (self.times[i][self.values[i].index(max_value)], max_value)
                    )
                else:
                    rez[i].append(
                        (self.times[i][self.values[i].index(-max_value)], -max_value)
                    )
            for j in range(len(rez[i])):
                k = rez[i][j][0]
                p = int(0.15 * len(self.times[i]) / len(self.nulls[i]))
                idx = (max((k - p, 0)), min(k + p + 1, len(self.times[i]) - 1))

                value, error = error_extrems(
                    self.times[i][idx[0] : idx[1]],
                    self.values[i][idx[0] : idx[1]],
                )
                new_rez[i].append(value)
                error_rez[i].append(error)
        return new_rez, error_rez

    def find_amplitude(self):
        rez = []
        self.amplitude_error = []
        for i in range(self.channels):
            try:
                # mean_value, eror = confidence_interval(
                #     [abs(self.extrems[i][j][1]) for j in range(len(self.extrems[i]))],
                #     relative_error=self.relative_erorr_y,
                #     absolute_error=self.absolute_error_y,
                #     confidence_probability=self.confidence_probability,
                # )
                mean_value, error = uneven_measurements(
                    [abs(self.extrems[i][j][1]) for j in range(len(self.extrems[i]))],
                    self.extrems_error[i],
                )
                print(
                    text_info(
                        f"Амплитуда {i + 1}-го канала определена по {len(self.extrems[i])} экстремумам"
                    )
                )
                self.amplitude_error.append(error)
                rez.append(mean_value)
            except Exception as e:
                # exc_type, exc_value, exc_traceback = sys.exc_info()

                # Форматируем traceback в строку
                error_message = traceback.format_exc()
                print(
                    text_error(
                        f"Не удалось определить амлитуду для {i}-го канала, ошибка:\n{e}\n{error_message}"
                    )
                )
                self.amplitude_error.append(0)
                rez.append(0)
        return rez

    def find_frequency(self):
        """Возвращает массив с частотами для каждого канала.

        Returns:
            list[float]: Массив с частотами для каждого канала.
        """
        rez = []
        frequency_error = []
        self.frequency_error = []
        for i in range(self.channels):
            flag_nulls, flag_extrem = False, False
            try:  # Определение частоты по нулям
                period, eror_period = confidence_interval(
                    [
                        (2 * abs(self.nulls[i][j + 1] - self.nulls[i][j]))
                        for j in range(len(self.nulls[i]) - 1)
                    ],
                    confidence_probability=self.confidence_probability,
                    absolute_error=2 * 2**0.5 * self.delta_time[i],
                )
                freq_nulls = 1 / period
                flag_nulls = True
                eror_freq_nulls = eror_period / period**2
            except Exception as e:
                print(
                    text_error(
                        f"Не удалось определить частоту {i+1}-того канала по нулям сигнала. {e}"
                    )
                )
            try:  # Определение частоты по экстремумам
                period, eror_period = confidence_interval(
                    [
                        (2 * abs(self.extrems[i][j + 1][0] - self.extrems[i][j][0]))
                        for j in range(len(self.extrems[i]) - 1)
                    ],
                    confidence_probability=self.confidence_probability,
                    absolute_error=2 * 2**0.5 * self.delta_time[i],
                )
                freq_extrem = 1 / period
                eror_freq_extrems = eror_period / period**2
                flag_extrem = True
            except Exception as e:
                print(
                    text_error(
                        f"Не удалось опеределить частоту {i+1}-того канала по экстремумам. {e}"
                    )
                )
            if flag_nulls and flag_extrem:
                freq, freq_error = uneven_measurements(
                    (freq_nulls, freq_extrem), (eror_freq_nulls, eror_freq_extrems)
                )
                rez.append(freq)
                frequency_error.append(freq_error)
                print(
                    text_info(
                        f"Частота {i}-того канала найдена как среднее между частотой найденой по экстремумам и по нулям сигнала."
                    )
                )
            elif flag_nulls:
                rez.append(freq_nulls)
                frequency_error.append(eror_freq_nulls)
            elif flag_extrem:
                rez.append(freq_extrem)
                frequency_error.append(eror_freq_extrems)
            else:
                rez.append(0)
                frequency_error.append(0)
        return (rez, frequency_error)

    def find_simplify_sig(self, number_chenel):
        """Возвращает массив с коэффициентами k и b прямых, проходящих через соседние экстремумы

        Returns:
            list[lsit[k: float, b: float]]: Массив с коэффициентами для каждого канала.
        """
        rez = []
        rez_error = []
        for i in range(len(self.extrems[number_chenel]) - 1):
            x_1, x_2 = (
                self.extrems[number_chenel][i][0],
                self.extrems[number_chenel][i + 1][0],
            )
            y_1, y_2 = (
                self.extrems[number_chenel][i][1],
                self.extrems[number_chenel][i + 1][1],
            )
            k = (y_2 - y_1) / (x_2 - x_1)
            delta_k = (
                (y_1 / (x_2 - x_1) * self.relative_erorr_y * y_2) ** 2
                + (y_2 / (x_2 - x_1) * self.relative_erorr_y * y_1) ** 2
                + 2
                * (((y_2 - y_1) / (x_2 - x_1) * self.delta_time[number_chenel])) ** 2
            ) ** 0.5
            b = y_1 - k * x_1 + self.new_amp_nulls[number_chenel]
            delta_b = (
                (self.relative_erorr_y * y_1) ** 2
                + (delta_k * x_1) ** 2
                + (self.delta_time[number_chenel] * k) ** 2
            ) ** 0.5
            rez.append((k, b))
            rez_error.append((delta_k, delta_b))
        return (rez, rez_error)

    def find_shift(self, number_chenel_1, number_chenel_2):
        flag = False
        rez = 0
        for i in range(len(self.simplify_sig[number_chenel_1])):
            for j in range(len(self.simplify_sig[number_chenel_2])):
                if (
                    self.simplify_sig[number_chenel_1][i][0]
                    * self.simplify_sig[number_chenel_2][j][0]
                    > 0
                    and flag
                    and abs(
                        self.simplify_sig[number_chenel_1][i][1]
                        / self.simplify_sig[number_chenel_1][i][0]
                        - self.simplify_sig[number_chenel_2][j][1]
                        / self.simplify_sig[number_chenel_2][j][0]
                    )
                    < rez
                ):
                    rez = abs(
                        self.simplify_sig[number_chenel_1][i][1]
                        / self.simplify_sig[number_chenel_1][i][0]
                        - self.simplify_sig[number_chenel_2][j][1]
                        / self.simplify_sig[number_chenel_2][j][0]
                    )
                    delta_rez = (
                        (
                            self.simplify_sig_errors[number_chenel_1][i][1]
                            / self.simplify_sig[number_chenel_1][i][0]
                        )
                        ** 2
                        + (
                            self.simplify_sig_errors[number_chenel_2][j][1]
                            / self.simplify_sig[number_chenel_2][j][0]
                        )
                        ** 2
                        + (
                            self.simplify_sig[number_chenel_1][i][1]
                            * self.simplify_sig_errors[number_chenel_1][i][0]
                            / self.simplify_sig[number_chenel_1][i][0] ** 2
                        )
                        ** 2
                        + (
                            self.simplify_sig[number_chenel_2][j][1]
                            * self.simplify_sig_errors[number_chenel_2][j][0]
                            / self.simplify_sig[number_chenel_2][j][0] ** 2
                        )
                        ** 2
                    ) ** 0.5
                elif self.simplify_sig[number_chenel_1][i][0] * self.simplify_sig[
                    number_chenel_2
                ][j][0] > 0 and not (flag):
                    rez = abs(
                        self.simplify_sig[number_chenel_1][i][1]
                        / self.simplify_sig[number_chenel_1][i][0]
                        - self.simplify_sig[number_chenel_2][j][1]
                        / self.simplify_sig[number_chenel_2][j][0]
                    )
                    delta_rez = (
                        (
                            self.simplify_sig_errors[number_chenel_1][i][1]
                            / self.simplify_sig[number_chenel_1][i][0]
                        )
                        ** 2
                        + (
                            self.simplify_sig_errors[number_chenel_2][i][1]
                            / self.simplify_sig[number_chenel_2][j][0]
                        )
                        ** 2
                        + (
                            self.simplify_sig[number_chenel_1][i][1]
                            * self.simplify_sig_errors[number_chenel_1][i][0]
                            / self.simplify_sig[number_chenel_1][i][0] ** 2
                        )
                        ** 2
                        + (
                            self.simplify_sig[number_chenel_2][j][1]
                            * self.simplify_sig_errors[number_chenel_2][j][0]
                            / self.simplify_sig[number_chenel_2][j][0] ** 2
                        )
                        ** 2
                    ) ** 0.5
                    flag = True
        if not (flag):
            try:
                for i in range(len(self.extrems[number_chenel_1])):
                    for j in range(len(self.extrems[number_chenel_2])):
                        if (
                            (
                                self.extrems[number_chenel_1][i][1]
                                + self.new_amp_nulls[number_chenel_1]
                            )
                            * (
                                self.extrems[number_chenel_2][j][1]
                                + self.new_amp_nulls[number_chenel_2]
                            )
                            > 0
                            and abs(
                                self.extrems[number_chenel_1][i][0]
                                - self.extrems[number_chenel_2][j][0]
                            )
                            < rez
                            and flag
                        ):
                            rez = abs(
                                self.extrems[number_chenel_1][i][0]
                                - self.extrems[number_chenel_2][j][0]
                            )
                            delta_rez = (
                                self.delta_time[number_chenel_1] ** 2
                                + self.delta_time[number_chenel_2] ** 2
                            ) ** 0.5
                        elif (
                            self.extrems[number_chenel_1][i][1]
                            + self.new_amp_nulls[number_chenel_1]
                        ) * (
                            self.extrems[number_chenel_2][j][1]
                            + self.new_amp_nulls[number_chenel_2]
                        ) > 0 and not (
                            flag
                        ):
                            rez = abs(
                                self.extrems[number_chenel_1][i][0]
                                - self.extrems[number_chenel_2][j][0]
                            )
                            delta_rez = (
                                self.delta_time[number_chenel_1] ** 2
                                + self.delta_time[number_chenel_2] ** 2
                            ) ** 0.5
                            flag = True
            except Exception as e:
                print(text_error(f"Что-то пошло не так в функции find_shift, {e}"))
        try:
            (
                self.phase_shift_error[number_chenel_1][number_chenel_2],
                self.phase_shift_error[number_chenel_2][number_chenel_1],
            ) = (delta_rez, delta_rez)
        except Exception as e:
            print(
                text_error(
                    f"При расчте погрешности в функции find_shift что-то пшло не так! В качестве погрешности был возвращен 0. Ошибка: {e}"
                )
            )
            (
                self.phase_shift_error[number_chenel_1][number_chenel_2],
                self.phase_shift_error[number_chenel_2][number_chenel_1],
            ) = (0, 0)
        return rez

    def ploter(self, save: bool = False, show: bool = True, info_bar: bool = True):
        """Выводи осцилограму со всеми каналами, сдвигами фаз, амплитудами, экстремумами и смещениями сигнала."""
        shifts_collor = tuple(
            tuple(
                tuple(random.uniform(0.2, 0.5) for k in range(3))
                for j in range(self.channels - i)
            )
            for i in range(self.channels)
        )
        plt.title(self.name)
        for i in range(self.channels):
            color = [random.uniform(0.3, 0.45) for i in range(3)]
            index_max_color = color.index(max(color))
            color[index_max_color] = random.uniform(0.6, 0.85) * (
                1.1 if (index_max_color in (0, 2)) else 0.85
            )
            revers_color = tuple((1 - i) for i in color)
            plt.scatter(
                self.times[i],
                self.values[i],
                s=2,
                color=color,
            )
            plt.plot(
                self.times[i],
                self.values[i],
                lw=1,
                color=color,
                label=f"""Канал {i + 1}:
Частота: {self.frequency[i]:.2e}
Амплитуда: {self.amplitude[i]:.2e}""",
            )
            for j in range(self.channels):
                if i != j and j < i and len(self.nulls[i]) > 0:
                    plt.plot(
                        (
                            self.nulls[i][int(len(self.nulls[i]) / 2)],
                            self.nulls[i][int(len(self.nulls[i]) / 2)]
                            + self.phase_shift[i][j],
                        ),
                        (self.new_amp_nulls[i], self.new_amp_nulls[i]),
                        color=shifts_collor[i][j],
                        label=f"Сдвиг фаз к.{i + 1}, к.{j + 1}: {self.phase_shift[i][j]:.2e}",
                    )
                elif i != j and len(self.nulls[i]) > 0:
                    plt.plot(
                        [],
                        [],
                        color=shifts_collor[j][i],
                        label=f"Сдвиг фаз к.{i + 1}, к.{j + 1}: {self.phase_shift[i][j]:.2e}",
                    )
            try:
                extrem_x, extrem_y = zip(*self.extrems[i])
                plt.scatter(extrem_x, extrem_y, c=revers_color, zorder=10)
                # plt.errorbar(extrem_x, extrem_y, yerr=self.amplitude_error[i])
            except:
                pass
            plt.plot(
                (self.times[i][0], self.times[i][-1]),
                (
                    self.amplitude[i],
                    self.amplitude[i],
                ),
                ls="--",
                lw=0.5,
            )
        plt.grid()
        if info_bar:
            plt.legend(
                bbox_to_anchor=(0.5, -0.3), loc="lower center", ncol=self.channels
            )
        plt.tight_layout()
        if save:
            plt.savefig(f"{self.path}.png")
            print(
                text_info(f"Файл сохранен в этой дериктории с именем: {self.path}.png")
            )
        if show:
            plt.show()
        else:
            plt.close("all")

    def __str__(self):
        frequency_line = ""
        init_line = ""
        amplitudes_line = ""
        extrems_count_line = ""
        nulls_count_line = ""
        phase_shift_line = ""
        for i in range(self.channels):
            frequency_line += f"канал {i + 1}: ({(self.frequency[i] if hasattr(self.frequency, '__iter__') else self.frequency):.4e} +- {self.frequency_error[i]:.4e})Гц{',\t'if i !=  self.channels - 1 else ''}"
            init_line += (
                f"{len(self.values[i])}{', 'if i !=  self.channels - 1 else ''}"
            )
            amplitudes_line += f"канал {i + 1}: ({self.amplitude[i]:.4e} +- {self.amplitude_error[i]:.4e}){',\t'if i !=  self.channels - 1 else ''}"
            extrems_count_line += f"канал {i + 1}: {len(self.extrems[i])}{',\t'if i !=  self.channels - 1 else ''}"
            nulls_count_line += f"канал {i + 1}: {len(self.nulls[i])}{',\t'if i !=  self.channels - 1 else ''}"
            for j in range(i + 1):
                if i != j:
                    phase_shift_line += f"каналы {j + 1}, {i + 1}: ({self.phase_shift[i][j]:.4e} +- {self.phase_shift_error[i][j]:.4e})"
        return f"""
    Это экземпляр класса '{self.__class__.__name__}', с именем '{self.name}{',\t'if j !=  self.channels - 1 else ''}'.
    Путь до файла: '{self.path}'
    Осцилограмма имеет {self.channels} {'канал' if self.channels == 1 else 'канала'}, в которых по {init_line} значений x(времени) и y(измеряймой величины).
    Частота сигнала: {frequency_line}
    Амплитуда сигнала: {amplitudes_line}
    Число обнаруженых экстремумов сигнала: {extrems_count_line}
    Число обнаруженых нулей сигнала: {nulls_count_line}
    Рзность фаз каналов: {phase_shift_line}
               """

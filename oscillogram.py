import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os


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
        11. extrems(list[list[float]]) - экстремальные значения сигнала.
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
    ):
        self.path = path  # Инициализируетм путь до файла
        self.frequency = frequency if frequency else 0  # Инициализируем частоту
        self.name = name if name else path  # Инициализируем имя

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
                    if i > 4:  # Аналогично с предыдущим
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
        self.extrems = self.find_extrems()
        self.frequency = self.find_frequency()
        self.phase_shift = [
            [0 for j in range(self.channels)] for i in range(self.channels)
        ]  # Задаем матрицу для сдвигов фаз
        for i in range(self.channels):
            for j in range(self.channels):
                if i < j:
                    a = self.find_shift(i, j)
                    self.phase_shift[i][j], self.phase_shift[j][i] = a, a

    def find_nulls(self):  # Проверить на работоспособность
        """Метод нахождения нулей сигнала для каждого канала.
        Returns:
            list[list[float]] - список, списков временных координат, в которых сигнал переходит через 0, вычесляется как среднее от двух ближайших точек разного знака.
        """
        try:
            self.new_amp_nulls
        except:
            self.new_amp_nulls = [0 for i in range(self.channels)]
        rez = [[] for i in range(self.channels)]
        rez_1 = [[] for i in range(self.channels)]
        for i in range(self.channels):
            for j in range(0, len(self.values[i]) - 1):
                value = self.values[i][j]
                if (
                    value > self.new_amp_nulls[i]
                    and self.values[i][j + 1] < self.new_amp_nulls[i]
                ) or (
                    value < self.new_amp_nulls[i]
                    and self.values[i][j + 1] > self.new_amp_nulls[i]
                ):
                    rez[i].append(self.times[i][j])
            max_distance = 0
            for j in range(len(rez[i]) - 1):
                if abs(rez[i][j + 1] - rez[i][j]) > max_distance:
                    max_distance = abs(rez[i][j + 1] - rez[i][j])
            for j in range(len(rez[i]) - 1):
                sum = 0
                count = 0
                for k in range(j, len(rez[i]) - 1):
                    if abs(rez[i][j] - rez[i][k]) < max_distance / 2:
                        sum += rez[i][k]
                        count += 1
                        b = j + count
                if j == b - 1:
                    rez_1[i].append(sum / count)
        return rez_1

    def find_extrems(self):
        """Находит все экстремальные точки сигнала, а так же задает атрибут amplitude для каждого канала.

        Returns:
            list[list[tuple[float, float]]]: список в котором для каждого канала находятся списки кортежей с парами значений (время, значение).
        """
        rez = []
        self.amplitude = []
        for i in range(self.channels):
            max_value = max(abs(max(self.values[i])), abs(min(self.values[i])))
            rez.append([])
            singl = 1 if self.values[i][0] > 0 else -1
            actual_max_value = (self.times[i][0], self.values[i][0])
            flag = False
            self.amplitude.append(0)
            count = 0
            for j in range(len(self.values[i])):
                self.amplitude[i] += abs(actual_max_value[1]) if flag else 0
                count += 1 if flag else 0
                flag = False
                if singl * self.values[i][j] >= 0 and abs(self.values[i][j]) > abs(
                    actual_max_value[1]
                ):
                    actual_max_value = (self.times[i][j], self.values[i][j])
                elif singl * self.values[i][j] < 0:
                    singl *= -1
                    if abs(actual_max_value[1] / max_value) >= 0.95:
                        rez[-1].append(actual_max_value)
                        actual_max_value = (self.times[i][j], self.values[i][j])
                        flag = True
            try:
                self.amplitude[i] /= count
                print(
                    f"Для {i + 1}-того канала амплитуда найдена как среднее от {count}-х экстремумов."
                )
            except ZeroDivisionError:
                self.amplitude[i] = max_value
                print(
                    f"Для {i + 1}-того канала амплитуда найдена как максимумальное по модулю значение сигнала."
                )
        return rez

    def find_frequency(self):
        """Возвращает массив с частотами для каждого канала.

        Returns:
            list[float]: Массив с частотами для каждого канала.
        """
        rez = []
        for i in range(self.channels):
            flag_nulls, flag_extrem = False, False
            try:  # Определение частоты по нулям
                freq_nulls = (len(self.nulls[i]) - 1) / (
                    2 * abs(self.nulls[i][-1] - self.nulls[i][0])
                )
                flag_nulls = True
            except Exception as e:
                print(
                    f"Не удалось определить частоту {i+1}-того канала по нулям сигнала. {e}"
                )
            try:  # Определение частоты по экстремумам
                freq_extrem = (len(self.extrems[i]) - 1) / (
                    2 * abs(self.extrems[i][-1][0] - self.extrems[i][0][0])
                )
                flag_extrem = True
            except Exception as e:
                print(
                    f"Не удалось опеределить частоту {i+1}-того канала по экстремумам. {e}"
                )
            if flag_nulls and flag_extrem:
                rez.append((freq_nulls + freq_extrem) / 2)
                print(
                    f"Частота {i}-того канала найдена как среднее между частотой найденой по экстремумам и по нулям сигнала."
                )
            elif flag_nulls:
                rez.append(freq_nulls)
            elif flag_extrem:
                rez.append(freq_extrem)
            else:
                rez.append(0)
        return rez

    def find_shift(self, number_chenel_1, number_chenel_2):
        shift_extrems = 0
        shift_nulls = 0
        flag_extrems = False
        flag_nulls = False
        try:
            for i, value_1 in enumerate(self.extrems[number_chenel_1]):
                for j, value_2 in enumerate(self.extrems[number_chenel_2]):
                    if value_1[1] * value_2[1] > 0:
                        shift_actual = abs(
                            abs(
                                self.extrems[number_chenel_1][i][0]
                                - self.extrems[number_chenel_2][j][0]
                            )
                            - abs(i - j) / self.frequency[number_chenel_1] / 2
                        )
                    else:
                        shift_actual = abs(
                            abs(
                                self.extrems[number_chenel_1][i][0]
                                - self.extrems[number_chenel_2][j][0]
                            )
                            - (abs(i - j) / 2) / self.frequency[number_chenel_1]
                        )
                    shift_extrems += (
                        shift_actual
                        if 1 / self.frequency[number_chenel_1] / 4 > shift_actual
                        else 1 / self.frequency[number_chenel_1] / 2 - shift_actual
                    )
            shift_extrems /= len(self.extrems[number_chenel_1]) * len(
                self.extrems[number_chenel_2]
            )
            flag_extrems = True
        except Exception as e:
            print(
                f"Не удалось найти сдвиг фазы по экстремумам. Для каналов {number_chenel_1 + 1}, {number_chenel_2 + 1}. {e}"
            )
        try:
            for i, value_1 in enumerate(self.nulls[number_chenel_1]):
                for j, value_2 in enumerate(self.nulls[number_chenel_2]):

                    shift_actual = abs(
                        abs(
                            self.nulls[number_chenel_1][i]
                            - self.nulls[number_chenel_2][j]
                        )
                        - (abs(i - j)) / self.frequency[number_chenel_1] / 2
                    )
                    shift_nulls += (
                        shift_actual
                        if 1 / self.frequency[number_chenel_1] / 4 > shift_actual
                        else 1 / self.frequency[number_chenel_1] / 2 - shift_actual
                    )
            shift_nulls /= len(self.nulls[number_chenel_1]) * len(
                self.nulls[number_chenel_2]
            )
            flag_nulls = True

        except Exception as e:
            print(
                f"Не удалось найти сдвиг фазы по нулям. Для каналов {number_chenel_1 + 1}, {number_chenel_2 + 1}. {e}"
            )
        if flag_extrems and flag_nulls:
            shift = (shift_nulls + shift_extrems) / 2
            print(
                f"Разность фаз между каналами {number_chenel_1 + 1} и {number_chenel_2 + 1} была вычеслена как среднее между разностью полученой из экстремумов и нулей. Значение полученое по нулям: {shift_nulls}, значение полученое по экстремумам: {shift_extrems}"
            )
        elif flag_extrems:
            shift = shift_extrems
        elif flag_nulls:
            shift = shift_nulls
        else:
            print(
                f"При вычислении разности фаз между каналами {number_chenel_1 + 1} и {number_chenel_2 + 1}, что-то пошло не так, был возвращен 0."
            )
            shift = 0
        return shift

    def __str__(self):
        frequency_line = ""
        init_line = ""
        amplitudes_line = ""
        extrems_count_line = ""
        nulls_count_line = ""
        phase_shift_line = ""
        for i in range(self.channels):
            frequency_line += f"канал {i + 1}: {(self.frequency[i] if hasattr(self.frequency, '__iter__') else self.frequency):.4e}Гц{',\t'if i !=  self.channels - 1 else ''}"
            init_line += (
                f"{len(self.values[i])}{', 'if i !=  self.channels - 1 else ''}"
            )
            amplitudes_line += f"канал {i + 1}: {self.amplitude[i]:.4e}{',\t'if i !=  self.channels - 1 else ''}"
            extrems_count_line += f"канал {i + 1}: {len(self.extrems[i])}{',\t'if i !=  self.channels - 1 else ''}"
            nulls_count_line += f"канал {i + 1}: {len(self.nulls[i])}{',\t'if i !=  self.channels - 1 else ''}"
            for j in range(i + 1):
                if i != j:
                    phase_shift_line += (
                        f"каналы {j + 1}, {i + 1}: {self.phase_shift[i][j]:.4e}"
                    )
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


if __name__ == "__main__":
    print(
        oscillogram(
            path="C:/Users/Феодор/Documents/Конспекты/Теоретические основы электро и радиотехники/Лабы/Реакция простых цепей на гармоническое и импульсное воздействие/Логи/CR 0,1/1500.txt",
            name="Тест",
        )
    )

import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
import tabulate
import pip
pip.main(["install", "openpyxl"])


wave_lengths = [440*10**-9, 505*10**-9, 525*10**-9, 580*10**-9, 595*10**-9]

cols = [1, 2, 3]
data = pd.read_excel('./lab4.xlsx', usecols = (cols))
data.head()

I_gen = data['i'].tolist()
alphaplus_gen = data['alphaplus'].tolist()
alphaminus_gen = data['alphaminus'].tolist()

I_gen = list(map(float, I_gen))
alphaplus_gen = list(map(float, alphaplus_gen))
alphaminus_gen = list(map(float, alphaminus_gen))

I = I_gen[:11]
Verde = []
delta_Verde = []
d = 0.01



class Wave:
    def __init__(self, number):

        self.delta_alpha = []
        for val in range(0+number*11, 11+number*11, 1):
            self.delta_alpha.append(abs(alphaplus_gen[val] - alphaminus_gen[val])*60)  #заполнение дельта альфа
        print(self.delta_alpha)

        self.wave_length = wave_lengths[number]
        self.omega = 2 * np.pi * 3 * 10**8 / self.wave_length       #осознание длины волны

        self.B = []
        for val in I:
            self.B.append(val*77.8/10**3)                       #заполнение магнитной индукции

        self.psi = []
        for val in self.delta_alpha:                            #заполнение угла поворота плоскости
            self.psi.append((val / 2)*np.pi/60/180)

    def graph_analysis(self):

        self.x = np.array(self.B)
        self.y = np.array(self.psi)
                                                                                        #наименьшие квадраты для V
        self.A = np.vstack([self.x, np.ones(len(self.x))]).T
        self.m, self.c = np.linalg.lstsq(self.A, self.y, rcond=None)[0]

        sigmax = sum(np.square(self.x))/len(self.x) - np.square(np.mean(self.x))
        sigmay = sum(np.square(self.y))/len(self.y) - np.square(np.mean(self.y))

        self.sigma_m = np.sqrt((sigmax/sigmay - self.m*self.m)/9)


        plt.figure(1, figsize = (11.69, 8.27))
        plt.plot(self.x, self.x * self.m + self.c, 'r')
        plt.scatter(self.B, self.psi)
        plt.title(f"Зависимость \u03C8(B),рад. для \u03BB={self.wave_length * 10**9:.0f} нм")
        plt.show()

        print(self.sigma_m)

        delta_Verde.append(self.m*1.1*self.sigma_m / d)
        Verde.append(self.m/d)


waves = []

for val in range(0, 5, 1):
    temp = Wave(val)
    waves.append(temp)
    waves[val].graph_analysis()

omegas_squared = []

for val in waves:
    omegas_squared.append(val.omega)


for val in range(0,5,1):
    print(f"Постоянная Верде для длины волны \u03BB={wave_lengths[val]*10**9:.0f}нм: V = {Verde[val]:.0f} +- {delta_Verde[val]:.0f} рад/Тл*м")


plt.figure(1, figsize = (11.69, 8.27))

x = np.array(omegas_squared)
y = np.array(Verde)
A = np.vstack([x, np.ones(len(x))]).T
Faradey_const, c = np.linalg.lstsq(A, y, rcond=None)[0]

plt.scatter(omegas_squared, Verde)
plt.plot(x, Faradey_const*x + c, 'r')
plt.title(f"Зависимость V(omega),рад.")
plt.show()

delta_Faradey = sum(delta_Verde)/len(delta_Verde)/np.mean(omegas_squared)

print(f"Постоянная Фарадея Ф={Faradey_const} +- {delta_Faradey}")

N = 10**29
n = 1.6
m = 10**-30
c = 3 * 10**8
e = 1.6 * 10**-19
epsilon = 8.85 * 10**-12

omega_zero = np.sqrt(np.sqrt(N * e**3 / (2 * epsilon * n * m * m * c * Faradey_const)))

print(f"Собственная частота колебаний эклектрона = {omega_zero:.0f} Гц")
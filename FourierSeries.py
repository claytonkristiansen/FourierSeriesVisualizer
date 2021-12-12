import scipy
from scipy.integrate import quad
import numpy as np
import math
import matplotlib.pyplot as plt
import time

NUM_COEF = 200
FUNDAMENTAL_PERIOD = 20

def complex_quadrature(func, a, b):
    def real_func(t):
        return scipy.real(func(t))
    def imag_func(t):
        return scipy.imag(func(t))
    real_integral = quad(real_func, a, b)
    imag_integral = quad(imag_func, a, b)
    return (real_integral[0] + 1j*imag_integral[0])#, real_integral[1:], imag_integral[1:])

def FourierCoefficient(func, k, T):
    def integrand(t):
        return func(t) * np.exp(-1j * k * (2*math.pi / T) * t)
    res = complex_quadrature(integrand, 0, T)
    return (1/T) * res

def SynthesisEquation(fourierCoefficients, T, t, numCoef):
    sum = 0
    for k in range(0, numCoef):
        sum += fourierCoefficients[0][k] * np.exp(1j * -1 * k * (2*math.pi / T) * t)
    for k in range(1, numCoef):
        sum += fourierCoefficients[1][k] * np.exp(1j * k * (2*math.pi / T) * t)
    return sum

def Rect(t):
    if(abs(t % FUNDAMENTAL_PERIOD) < 10):
        return 1
    else:
        return 0

def Triangle(t):
    if(abs(t % FUNDAMENTAL_PERIOD) < 10):
        return (t % FUNDAMENTAL_PERIOD)
    else:
        return 0

# print(FourierCoefficient(Rect, 5, FUNDAMENTAL_PERIOD))


positiveCoefficients = []
negativeCoefficients = []
for i in range(0, NUM_COEF):
    positiveCoefficients.append(FourierCoefficient(Triangle, i, FUNDAMENTAL_PERIOD))
for i in range(0, NUM_COEF):
    negativeCoefficients.append(FourierCoefficient(Triangle, -1 * i, FUNDAMENTAL_PERIOD))

coefficientsPair = [negativeCoefficients, positiveCoefficients]

print("Val:", SynthesisEquation(coefficientsPair, FUNDAMENTAL_PERIOD, 5, 2))

x = np.linspace(0, 100, 400)
y1 = []
y2 = []
for xVal in x:
    y1.append(SynthesisEquation(coefficientsPair, FUNDAMENTAL_PERIOD, xVal, 1))
    y2.append(Rect(xVal))

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set(xlim =(-1, 100), ylim =(-0.1, 25))
line1, = ax.plot(x, y1, 'b-')

for numC in range(2, NUM_COEF):
    y1.clear()
    y2.clear()
    for xVal in x:
        y1.append(SynthesisEquation(coefficientsPair, FUNDAMENTAL_PERIOD, xVal, numC))
        y2.append(Rect(xVal))
    line1.set_ydata(y1)
    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(0.25)


# plt.plot(x, y1)
# plt.plot(x, y2)
# plt.show()
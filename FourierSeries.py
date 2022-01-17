import scipy
from scipy.integrate import quad
import numpy as np
import math
import matplotlib.pyplot as plt
import time

#These are the function definitions
#Feel free to change them or add additional ones as you see fit
#The FUNDAMENTAL_PERIOD must match that of the function you choose
#The functions must also be periodic because we are computing Fourier Series not Transforms
#You should use numbers larger than 10 so that Python can do the float math accurately
#The maximum number of coefficients is also around 200

#Global
FUNDAMENTAL_PERIOD = 40

def Rect(t):
    FUNDAMENTAL_PERIOD = 10
    relT = abs((t-5) % FUNDAMENTAL_PERIOD)
    if(relT < 5):
        return 10
    else:
        return 0

def RectSaw(t):
    FUNDAMENTAL_PERIOD = 40
    relT = abs(t % FUNDAMENTAL_PERIOD)
    if(relT <= 20):
        return (1)
    elif(relT > 20):
        return (-1)

def Triangle(t):
    FUNDAMENTAL_PERIOD = 20
    relT = abs(t % FUNDAMENTAL_PERIOD)
    if(relT < 10):
        return (t % FUNDAMENTAL_PERIOD)
    else:
        return 0

def IsocilesTriangle(t):
    FUNDAMENTAL_PERIOD = 40
    relT = abs(t % FUNDAMENTAL_PERIOD)
    if(relT >= 10 and relT <= 20):
        return (relT - 10)
    elif(relT > 20 and relT <= 30):
        return (30 - relT)
    else:
        return 0

def Tangent(t):
    return np.tan(t)

def AbsSin(t):
    FUNDAMENTAL_PERIOD = 10
    return np.abs(np.sin(t * (math.pi / 20)))

#Choose function from definitions above. Feel free to add your own
FUNC = Rect     
#Choose the speed at which the graph is updated
UPDATE_PERIOD = 0.25
SKIP_AMOUNT = 1
#Choose the number of Fourier coefficient calculated (this is the number to both the left and right of zero)
NUM_COEF = 200
#Choose graph bounds
X_L = -100
X_R = 100
Y_T = 20
Y_B = -20
#Choose the number of poitns 
NUM_POINTS = 400

#Calculated the definite complex integral of func from a to b
def complex_quadrature(func, a, b):
    def real_func(t):
        return np.real(func(t))
    def imag_func(t):
        return np.imag(func(t))
    real_integral = quad(real_func, a, b)
    imag_integral = quad(imag_func, a, b)
    return (real_integral[0] + 1j*imag_integral[0])

#Calulated the Fourier Coefficients using the Analysis Equation (4.13 in the notes)
def FourierCoefficient(func, k, T):
    def integrand(t):
        return func(t) * np.exp(-1j * k * (2*math.pi / T) * t)
    res = complex_quadrature(integrand, 0, T)
    return (1/T) * res

#Calcuates the equation in the time domain using the Synthesis Equation (4.10 in the notes)
def SynthesisEquation(fourierCoefficients, T, t, numCoef):
    sum = 0
    for k in range(0, numCoef):
        sum += fourierCoefficients[0][k] * np.exp(1j * -1 * k * (2*math.pi / T) * t)
    for k in range(1, numCoef):
        sum += fourierCoefficients[1][k] * np.exp(1j * k * (2*math.pi / T) * t)
    return sum


positiveCoefficients = []
negativeCoefficients = []
for i in range(0, NUM_COEF):
    positiveCoefficients.append(FourierCoefficient(FUNC, i, FUNDAMENTAL_PERIOD))
for i in range(0, NUM_COEF):
    negativeCoefficients.append(FourierCoefficient(FUNC, -1 * i, FUNDAMENTAL_PERIOD))

#Collection containing two lists, the first are the coefficients from k:(-oo, 0] and the second are the coefficients k:[0, oo)
coefficientsPair = [negativeCoefficients, positiveCoefficients]


#Everything from here down is plotting and the logic for updating the graph
x = np.linspace(X_L, X_R, NUM_POINTS)
y1 = []
y2 = []
for xVal in x:
    y1.append(SynthesisEquation(coefficientsPair, FUNDAMENTAL_PERIOD, xVal, 1))
    y2.append(FUNC(xVal))

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set(xlim =(X_L, X_R), ylim =(Y_B, Y_T))
line1, = ax.plot(x, np.real(y1), 'b-')

for numC in range(2, NUM_COEF, SKIP_AMOUNT):
    y1.clear()
    y2.clear()
    for xVal in x:
        y1.append(SynthesisEquation(coefficientsPair, FUNDAMENTAL_PERIOD, xVal, numC))
        y2.append(FUNC(xVal))
    line1.set_ydata(np.real(y1))
    text = (str(numC) + " coefficients used in sum")
    ax.set_title(text)
    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(UPDATE_PERIOD)

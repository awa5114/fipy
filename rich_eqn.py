import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math

class Layer:
    pass

class Image(Layer):
    pass

class CurveSet(Layer):
    def __init__(self, R_star, curves):
        self.R_star = R_star
        self.curves = curves  # A list of Curve objects
        self.figure, self.ax = plt.subplots(1)
        self.add_curves()

    def add_curves(self):
        for curve in self.curves:
            Theta, z_star = curve.calculate(curve.t_star, curve.C, self.R_star)
            self.ax.plot(Theta, z_star, color='blue')

    def show(self):
        plt.show()

class Curve:
    def __init__(self, t_star, C):
        self.t_star = t_star
        self.C = C

    def du_deta_func(self, eta, tau, rho):
        f = lambda x: np.exp(x ** 2) * math.erfc(x)

        A = rho * np.exp(-eta ** 2 * tau ** -1)
        B = 2 * np.exp(((eta + rho * tau) / tau ** 0.5) ** 2)
        C = (1 + rho ** -1) ** 0.5
        D = (eta - rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
        E = (eta + rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
        F = (eta - rho * tau) / tau ** 0.5
        G = (eta + rho * tau) / tau ** 0.5

        du_deta = A * (B - C * (f(D) - f(E)) + f(F) - f(G))
        return du_deta

    def u_func(self, eta, tau, rho):
        f = lambda x: np.exp(x ** 2) * math.erfc(x)

        A = 0.5 * np.exp(-eta ** 2 * tau ** -1)
        B = 2 * np.exp(((eta + rho * tau) / tau ** 0.5) ** 2)
        C = (eta - rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
        D = (eta - rho * tau) / tau ** 0.5
        E = (eta + rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
        F = (eta + rho * tau) / tau ** 0.5
        u = A * (B + f(C) - f(D) + f(E) - f(F))

        return u

    def calculate(self, t_star, C, R_star):
        m = 4 * C * (C - 1)
        rho = R_star / m
        tau = m * t_star

        z_star_list = []
        Theta_list = []
        eta_list = np.linspace(0, 1, 100)  # parameter

        for eta in eta_list:
            u = self.u_func(eta, tau, rho)
            du_deta = self.du_deta_func(eta, tau, rho)
            z_star = C ** -1 * (rho ** 2 * (1 + rho ** -1) * tau + rho * (2 + rho ** -1) * eta - np.log(u))
            Theta = C * (1 - (2 * rho + 1 - u ** -1 * du_deta) ** -1)
            z_star_list.append(z_star)
            Theta_list.append(Theta)

        return Theta_list, z_star_list

class Plot:
    def __init__(self, layers):
        self.layers = layers

curves=[]
for t_star in [0.5, 2, 4, 8, 16]:
    curves.append(Curve(t_star, 1.02))
fig3_top = CurveSet(0.2,  curves)
fig3_top.show()

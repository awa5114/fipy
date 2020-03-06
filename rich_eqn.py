import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math

f = lambda x: np.exp(x**2)*math.erfc(x)

def du_deta_func(eta, tau, rho):
    A = rho * np.exp(-eta ** 2 * tau ** -1)
    B = 2 * np.exp(((eta + rho * tau) / tau ** 0.5) ** 2)
    C = (1 + rho ** -1) ** 0.5
    D = (eta - rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
    E = (eta + rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
    F = (eta - rho * tau) / tau ** 0.5
    G = (eta + rho * tau) / tau ** 0.5

    du_deta = A * (B - C * (f(D) - f(E)) + f(F) - f(G))
    return du_deta

def u_func(eta, tau, rho):
    A = 0.5 * np.exp(-eta ** 2 * tau ** -1)
    B = 2 * np.exp(((eta + rho * tau) / tau ** 0.5) ** 2)
    C = (eta - rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
    D = (eta - rho * tau) / tau ** 0.5
    E = (eta + rho * (1 + rho ** -1) ** 0.5 * tau) / tau ** 0.5
    F = (eta + rho * tau) / tau ** 0.5
    u = A * (B + f(C) - f(D) + f(E) - f(F))

    return u

C=1.02
R_star=0.2
t_star = 8
eta_list=np.linspace(0,1,100) # parameter

m = 4*C*(C-1)
rho = R_star/m
tau = m*t_star

z_star_list=[]
Theta_list=[]

for eta in eta_list:
    u = u_func(eta, tau, rho)
    du_deta = du_deta_func(eta, tau, rho)
    z_star = C**-1*(rho**2*(1+rho**-1)*tau+rho*(2+rho**-1)*eta-np.log(u))
    Theta=C*(1-(2*rho+1-u**-1*du_deta)**-1)
    z_star_list.append(z_star)
    Theta_list.append(Theta)

img = mpl.image.imread(r'lit/R_0_2.JPG')
fig, ax = plt.subplots()
ax.set_xlim(0,1)
ax.set_ylim(5,0)
ax.imshow(img, extent=(0, 1.0, 5, 0), aspect=0.2)
ax.plot(Theta_list, z_star_list, dashes=(5,5,5,5))
plt.show()
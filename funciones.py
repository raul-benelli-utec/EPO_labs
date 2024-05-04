import numpy as np
from scipy.optimize import newton

# Definir la función que queremos encontrar la raíz
def equacion_media_onda(alpha, theta, omega, R, C):  
    return np.sin(alpha) - (np.sin(theta) * np.exp(-((2*np.pi + alpha - theta) / (omega * R * C))))

def equacion_onda_comleta(alpha, theta, omega, R, C):  
    return (np.sin(theta) * np.exp(-((np.pi + alpha - theta) / (omega * R * C))))-np.sin(alpha)

# Definir los valores conocidos

# Encontrar la raíz de la ecuación con el método de Newton-Raphson
initial_guess = 0  # Supongamos un valor inicial para alpha
def hallar_alpha(theta, omega, R, C,initial_guess = 0,media_onda=True):
    if media_onda:
        alpha = newton(equacion_media_onda, initial_guess, args=(theta, omega, R, C))
        return alpha
    else:
        alpha = newton(equacion_onda_comleta, initial_guess, args=(theta, omega, R, C))
        return alpha
    
def imprimir_valor(valor,unidad="",cifras=3):
    if unidad=="":
        return f"{round(valor,cifras)}"
    if unidad=="rad":
        return f"{round(valor,cifras)} rad"
    if valor<1:
        return f"{round(valor*1000,cifras)} m{unidad}"
    else:
        return f"{round(valor,cifras)} {unidad}"
    

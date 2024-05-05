import numpy as np
import matplotlib.pyplot as plt
import datos_circuito as circuito
from scipy.integrate import quad
from funciones import *
import numpy as np

import scipy.integrate as integrate

dos_pi=np.pi*2
Vm=circuito.Vrms*np.sqrt(2)
DVo=0.2*Vm

print("\n\n")
#Valor de C para obtener DVo=0.2*Vm
C_teorica=1/(circuito.frec*circuito.r*DVo/Vm)
print(f"El valor teorico obtenido para C= {C_teorica*10e5} pico faradios")

print("\n\n")
#Ingresar el valor de C a utilizar considerando tolerancia y valores comerciales
C_real=input("Ingrese el valore de C a utilizar en pico faradios:")
C_real=float(C_real)
print(f"Valor de C real: {C_real} pico faracdios")

C_real=C_real*1e-6
w=dos_pi*circuito.frec
wrc=w*int(circuito.r)*C_real
print(f"wrc {wrc}")
print(f"wcvm {w*C_real*Vm}")
theta=-np.arctan(wrc)+np.pi
print(f"theta {theta}")
alpha=hallar_alpha(theta,w,circuito.r,C_real,0)
d_vo=Vm*(1-np.sin(alpha))
# Función para la tensión de salida Vo(t)
def Vo(t):
    wt_mod_period = w*t % dos_pi  # Obtener el valor de wt dentro del período
    
    if (dos_pi + alpha) <= (wt_mod_period+dos_pi) < (dos_pi + theta):
        return Vm*np.sin(w*t)
    else:
        if 0<=wt_mod_period<alpha:
            wt_mod_period=wt_mod_period+dos_pi
        return Vm*np.sin(theta)*np.exp(-((wt_mod_period)-theta)/(wrc))

# Función para el cuadrado de Vo(t)
def Vo_squared(t):
    return Vo(t) ** 2

# Función para calcular el valor RMS de Vo(t)
def Vo_rms():
    # Calcula la integral numérica para encontrar el valor RMS
    integral, _ = quad(Vo_squared, 0, 2 * np.pi / w, )
    
    return np.sqrt(integral / (2 * np.pi / w))

def Vo_rms_dos():
    #Función secundaria para calcular Vo rms, comprobamos que devuelve los mismos resultados
    f = lambda t : (Vm*np.sin(theta)*np.exp(-((w*t)-theta)/(wrc)))** 2
    g= lambda t :(Vm*np.sin(w*t))** 2
    a = theta/w
    b = dos_pi_m_alpha/w
    c=dos_pi_m_alpha/w
    d=(dos_pi+theta)/w

    I = integrate.quad(f,a,b)
    I2= integrate.quad(g,c,d)
    
    Vo_rms= np.sqrt((1/0.02)*(I[0]+I2[0]))
    print(f"Vo rms dos: {Vo_rms}")
    return Vo_rms
#Vo_rms_dos()

vo_rms= Vo_rms()

def Ir(t):
    return Vo(t)/circuito.r

def Ic(t):
    wt_mod_period = w*t % dos_pi  # Obtener el valor de wt dentro del período
    
    if (dos_pi + alpha) <= (wt_mod_period+dos_pi) < (dos_pi + theta):
        return w*C_real*Vm*np.cos(w*t)
    else:
        if 0<=wt_mod_period<alpha:
            wt_mod_period=wt_mod_period+dos_pi
        return -(Vm*np.sin(theta)/circuito.r)*np.exp(-((wt_mod_period)-theta)/(wrc))
    
T_rad=dos_pi
#periodo en segundos 
T_s=dos_pi/w
#medio periodo en segundos 
medioT=T_s/2
#cuarto periodo en segundos 
cuartoT=medioT/2

t_div_4=dos_pi/(w*4)
dos_pi_m_alpha=dos_pi+alpha
tc=(((dos_pi/(w)))+t_div_4)-(dos_pi_m_alpha/w)

def Ic_rms():
    #potencia suministrada por la fuente al capacitor en el ciclo de conduccion
    #calculo rms en el periodo
    g= lambda wt :(w*C_real*Vm*np.cos(wt))** 2
    c=alpha
    d=theta

    I2= integrate.quad(g,c,d)
    
    ic_rms= np.sqrt((1/dos_pi)*(I2[0]))
    print(f"Ic rms: {ic_rms}")
    return ic_rms

ic_rms=Ic_rms()

Ir_rms=vo_rms/circuito.r
print(f"Ir rms {Ir_rms}")
Ir_p=Vm*(np.sin(alpha)/circuito.r)
Ic_p=Vm*np.cos(alpha)*w*C_real
Id_p=Ic_p+Ir_p 
Is_rms=Ir_rms+ic_rms
def Vs(t):
    return Vm * np.sin(w*t)

def Fp():
    pr=(vo_rms ** 2)/(circuito.r)
    ps=circuito.Vrms*Is_rms
    return pr/ps


def imprimir_resultados():
    print("Rectificador de media onda con filtro capacitivo.")
    print("\n")
    print("Datos:")
    print("Vs rms: "+imprimir_valor(circuito.Vrms,"V"))
    print("Frecuencia : ", imprimir_valor(circuito.frec,"Hz"))
    print("Resistencia: ", imprimir_valor(circuito.r, "Ω"))
    print("Periodo",imprimir_valor(T_s,"s"))
    print("\n")
    print("Valores calculados:")
    print("Alpha: "+imprimir_valor(alpha,"rad"))
    print("Theta: ",imprimir_valor(theta,"rad"))
    print("Delta Vo calculado: "+ imprimir_valor(d_vo,"V"))
    print(f"Vm {round(Vm,2)} V")
    print("Vs(theta): "+ imprimir_valor(Vm*np.sin(theta),"V"))
    print("Vo rms: "+imprimir_valor(vo_rms,"V"))
    print("Is_rms: "+ imprimir_valor(Is_rms,"A"))
    print("Ic_rms: "+ imprimir_valor(ic_rms,"A"))
    print("Ir_rms: "+ imprimir_valor(Ir_rms,"A"))
    print("Ir pico: "+ imprimir_valor(Ir_p,"A"))
    print("Ic pico: "+ imprimir_valor(Ic_p,"A"))
    print("Id pico: "+ imprimir_valor(Id_p,"A"))
    print(f"T inicio carga capacitor: "+ imprimir_valor((alpha/w),"s"))
    print(f"T carga capacitor: "+ imprimir_valor(tc,"s"))
    print(f"T descarga capacitor: "+ imprimir_valor((T_s-tc),"s"))
    print("Factor de potencia: "+ imprimir_valor(Fp()))


imprimir_resultados()


def graficar(periodos=1,bloquear=False,nombre="False",tansparente=False):
    puntos=5000
    rango=T_s*periodos*puntos
    x = list(range(0, int(rango+1), 1))

    fig, (graf_i, graf_v) = plt.subplots(2)
    graf_v.plot([i/puntos for i in x], [Vs(i/puntos) for i in x],label ='Vs(t)')
    graf_v.plot([i/puntos for i in x], [Vo(i/puntos) for i in x],label ='Vo(t)')
    graf_i.plot([i/puntos for i in x], [Ic(i/puntos) for i in x],label ='Ic(t)')
    graf_i.plot([i/puntos for i in x], [Ir(i/puntos) for i in x],label ='Ir(t)')

    for i in range(0,4):
        graf_i.axvline(x=((dos_pi/(w))*i)+t_div_4, color='black', linestyle='--',linewidth=0.5)
    graf_i.axvline(x=(dos_pi_m_alpha/w), color='black', linestyle='--',linewidth=0.5)
    for i in range(0,4):
        graf_v.axvline(x=((dos_pi/(w))*i)+t_div_4, color='black', linestyle='--',linewidth=0.5)

    #graf_v.axvline(x=0.04, color='green', linestyle='--', label='Línea vertical en x=0.4')
    #t inicia carga del capacitor
    graf_v.axvline(x=(alpha/w), color='black', linestyle='--',linewidth=0.6)
    graf_v.text((alpha/w), -27, 'α', rotation=0, va='bottom')

    graf_i.axvline(x=(alpha/w), color='black', linestyle='--',linewidth=0.6)
    graf_i.text((alpha/w), -0.1, 'α', rotation=0, va='bottom')

    graf_v.axhline(y=0,color='black', linestyle='--',linewidth=0.5)
    graf_v.axhline(y=Vm, color='black', linestyle='--',linewidth=0.5)
    #delta Vo
    ajuste=T_s*periodos/50
    graf_v.text(-ajuste, Vm, f'{round(Vm,2)}', rotation=0, va='bottom', ha='right')
    graf_v.axhline(y=Vm-d_vo, color='black', linestyle='--',linewidth=0.5)
    graf_v.text(-ajuste, Vm-d_vo-3, f'{round(Vm-d_vo,2)}', rotation=0, va='bottom', ha='right')
    #delimitar valores de los ejes
    graf_i.set_xlim([0, (T_s*periodos)])
    graf_v.set_xlim([0, (T_s*periodos)])
    #sub titulos
    graf_i.set_title('Subgráfico 1: Corriente')
    graf_v.set_title('Subgráfico 2: Voltajes')
    #Unidades de los ejes
    graf_i.set_xlabel('Tiempo (s)')
    graf_i.set_ylabel('Corriente (A)')
    graf_v.set_xlabel('Tiempo (s)')
    graf_v.set_ylabel('Voltaje (V)')
    # Añadir un título al gráfico completo
    fig.suptitle('Rectificador de media onda con filtro capacitivo', fontsize=16)

    plt.legend()
    graf_i.legend()
    graf_v.legend()
    plt.tight_layout()
    if nombre:
        plt.savefig(f'{nombre}.png', dpi=500, bbox_inches='tight', transparent=tansparente)
    plt.show(block=bloquear) 

#graficar(6.5,nombre="rect_media_onda_6 periodos")
#graficar(1,True,nombre="rect_media_onda_ periodo")
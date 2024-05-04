import numpy as np
import matplotlib.pyplot as plt
import datos_circuito as circuito
from scipy.integrate import quad
from funciones import *

dos_pi=np.pi*2
Vm=circuito.Vrms*np.sqrt(2)
DVo=0.2*Vm

#Valor de C para obtener DVo=0.2*Vm
C_teorica=1/(2*circuito.frec*circuito.r*DVo/Vm)
print("\n\n")
print(f"El valor teorico obtenido para C= {C_teorica*10e5} pico faradios")

#Ingresar el valor de C a utilizar considerando tolerancia y valores comerciales
print("\n\n")
C_real=input("Ingrese el valore de C a utilizar en pico faradios:")
C_real=float(C_real)
print(f"El valor ingresado es {C_real} pico faracdios")
#expresar en faradios
C_real=C_real*1e-6

#cálculos del circuito
w=dos_pi*circuito.frec
wrc=w*int(circuito.r)*C_real
theta=-np.arctan(wrc)+np.pi
alpha=hallar_alpha(theta,w,circuito.r,C_real,0,media_onda=False)
d_vo=Vm*(1-np.sin(alpha))


def Vo(t):
    #funcion a trozos en base a las caracteristicas del circuito
    wt_mod_period = w*t % np.pi  # Obtener el valor de wt dentro del período
    if (np.pi + alpha) <= (wt_mod_period+np.pi) < (np.pi + theta):
        return abs(Vm*np.sin(w*t))
    else:
        if 0<=wt_mod_period<alpha:
            wt_mod_period=wt_mod_period+np.pi
        return Vm*np.sin(theta)*np.exp(-((wt_mod_period)-theta)/(wrc))
    
def Vo_rms(w):
    # Define la función para el cuadrado de Vo(t)
    def Vo_cuadrado(t):
        return Vo(t) ** 2
    
    # Calcula la integral numérica para encontrar el valor RMS
    integral, _ = quad(Vo_cuadrado, 0, 2 * np.pi / w)
    vo_rms = np.sqrt(integral / (2 * np.pi / w))
    return vo_rms

def Vs(t):
    #Valor en función del tiempo
    return Vm * np.sin(w*t)

def I_fuente(t):
    return Vs(t)/circuito.r

vo_rms = Vo_rms(w)

def Ir(t):
    return Vo(t)/circuito.r

def Ic(t):
    #funcion a trozos en base a las caracteristicas del circuito
    wt_mod_period = w*t % dos_pi  # Obtener el valor de wt dentro del período
    
    if (dos_pi + alpha) <= (wt_mod_period+dos_pi) < (dos_pi + theta):
        return w*C_real*Vm*np.cos(w*t)
    else:
        if 0<=wt_mod_period<alpha:
            wt_mod_period=wt_mod_period+dos_pi
        return -(Vm*np.sin(theta)/circuito.r)*np.exp(-((wt_mod_period)-theta)/(wrc))   

def I_fuente(t):
    return Vs(t)/circuito.r




Ir_p=Vm*(np.sin(alpha)/circuito.r)
Ic_p=Vm*np.cos(alpha)*w*C_real
Id_p=Ic_p+Ir_p 
Is_rms=vo_rms/circuito.r

dos_pi_m_alpha=dos_pi+alpha

pr=(circuito.Vrms ** 2)/(circuito.r)
ps=circuito.Vrms*Is_rms

T_rad=dos_pi
t_div_4=dos_pi/(w*4)
dos_pi_m_alpha=dos_pi+alpha
tc=(((dos_pi/(w)))+t_div_4)-(dos_pi_m_alpha/w)

#periodo en segundos 
T_s=dos_pi/w
#medio periodo en segundos 
medioT=T_s/2
#cuarto periodo en segundos 
cuartoT=medioT/2


def Fp():
    return pr/ps
    

def imprimir_resultados():
    print("Rectificador de onda completa con filtro capacitivo.")
    print("\n")
    print("Datos:")
    print("Vs rms: "+imprimir_valor(circuito.Vrms,"V"))
    print("Frecuencia : ", imprimir_valor(circuito.frec,"Hz"))
    print("Resistencia: ", imprimir_valor(circuito.r, "Ω"))
    print("Periodo",imprimir_valor(T_s,"s"))
    print("\n")
    print("Valores calculados:")
    print("Alpha: "+imprimir_valor(alpha,"rad"))
    print("Delta Vo calculado: "+ imprimir_valor(d_vo,"V"))
    print(f"Vm {round(Vm,2)} V")
    print("Vs(theta): "+ imprimir_valor(Vm*np.sin(theta),"V"))
    
    print("Vo rms: "+imprimir_valor(vo_rms,"V"))
    print("Is_rms: "+ imprimir_valor(Is_rms,"A"))
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
    graf_i.plot([i/puntos for i in x], [I_fuente(i/puntos) for i in x],label ='Is(t)')

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
    fig.suptitle('Rectificador de onda completa con filtro capacitivo', fontsize=16)

    plt.legend()
    graf_i.legend()
    graf_v.legend()
    plt.tight_layout()
    if nombre:
        plt.savefig(f'{nombre}.png', dpi=500, bbox_inches='tight', transparent=tansparente)
    plt.show(block=bloquear) 

graficar(6.5,nombre="rectificador onda completa 6 periodos")
graficar(1,True,nombre="rectificador onda completa 1 periodo")


"""
print(f"cuartoT {cuartoT}, medioT {medioT}")
for i in range(0,4):
    ax1.axvline(x=((T)*i)+cuartoT, color='pink', linestyle=':')
ax1.axvline(x=(dos_pi_m_alpha/w), color='blue', linestyle='--')
for i in range(0,4):
    ax2.axvline(x=((dos_pi/(w))*i)+cuartoT, color='green', linestyle=':')

ax2.axvline(x=0.04, color='red', linestyle='--',linewidth=0.5)
ax2.axvline(x=0.02, color='red', linestyle='--',linewidth=0.5)
ax2.axvline(x=(dos_pi_m_alpha/w), color='blue', linestyle=':',linewidth=0.5)
ax2.axvline(x=((1.5*dos_pi/w)+(alpha/w)), color='blue', linestyle=':',linewidth=1.5)
ax2.axvline(x=(((5*dos_pi)/4*w)), color='blue', linestyle=':',linewidth=1.5)
ax2.axhline(y=0, color='blue',linestyle='-')
ax2.axhline(y=Vm, color='green', linestyle=':')
ax2.axhline(y=Vm-d_vo, color='green', linestyle=':')
ax1.set_xlim([0, 0.09])
ax2.set_xlim([0, 0.09])

print(f"W: {w}, dos_pi/w {dos_pi/w}, dos_pi_m_alpha/w {dos_pi_m_alpha/w}")
bbox_props = dict(boxstyle='round,pad=0.5', fc='white', ec='black', lw=1, alpha=0.5)
ax1.set_title('Subgráfico 1: Corriente')
ax2.set_title('Subgráfico 2: Voltajes')

# Añadir un título al gráfico completo
fig.suptitle('Rectificador de onda completa con filtro capacitivo', fontsize=16)

plt.legend()
ax1.legend()
ax2.legend()
plt.tight_layout()
plt.show()
#Gráficar
x = list(range(0, 401, 1))

fig, (ax1, ax2) = plt.subplots(2)

ax2.plot([i/5000 for i in x], [Vs(i/5000) for i in x],label ='Vs(t)')
ax2.plot([i/5000 for i in x], [Vo(i/5000) for i in x],label ='Vo(t)')
ax1.plot([i/5000 for i in x], [Ic(i/5000) for i in x],label ='Ic(t)')
ax1.plot([i/5000 for i in x], [Ir(i/5000) for i in x],label ='Ir(t)')
ax1.plot([i/5000 for i in x], [I_fuente(i/5000) for i in x],label ='Is(t)')

#ax.axvline(x=0.02, color='red', linestyle='--', label='Línea vertical en x=0.2')
#lineas verticales 
#periodo en segundos 
T=dos_pi/w
#medio periodo en segundos 
medioT=T/2
#cuarto periodo en segundos 
cuartoT=medioT/2


"""
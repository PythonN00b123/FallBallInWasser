import matplotlib.pyplot as plt
import numpy as np
from celluloid import Camera



#### TOPF U.Ä. ####
hightPot = 65
widthPot = 50
#hightBod = 8                #Höhe Körper
hightWat = 45                #Wasserhöhe

##### WERTE #####
fps = 30
t = 1/fps                   #Zeitintervall in secs

g = -9.81                      #g in m/s^2
h = 100 
r = 4
V_k = np.pi * r**2
roh_k = 0.4                   #Dichte Körper
roh_w = 1                   #Dichte Wasser
m_k = roh_k * V_k           #Masse des Körpers

######## ENERGIE ###########
E_max = abs(m_k * g * h)
E_plac = [0.02 * E_max,0.07 * E_max,0.12 * E_max]
x_tit = ["E_pot","E_kin","E_ges"]

########## GRAPH ##########
fig, ax = plt.subplots(1,2)
camera = Camera(fig)

ax[0].set_ylim(0,1.15*h)
ax[0].set_xlim(-((1.3*h/2)-widthPot),1.3*h/2)
ax[1].set_title("Energie in kJ")
ax[1].set_ylim(0,1.2*E_max/1000)

V_sub = []
#a_res = -9.81
s=[0]
v = 0
for i in range(0,4* round(np.sqrt((h-r)*2/-g) * fps)):
    
    P = [0.5*widthPot,h]

    #### WASSERSTAND ####
    if h-r > hightWat:
        V_sub.append(0) #V = Volumen submerged

    elif h > hightWat and h-r <= hightWat:
        
        a = [np.sqrt(-((P[1]-hightWat)**2 - r**2)) + P[0], hightWat]
        b = [2*P[0]-a[0], hightWat]

        x = np.arange(b[0],a[0],0.01)
        y = np.array(-np.sqrt(r**2-(x-P[0])**2)+P[1])
        V_sub.append(round(-(np.trapz(y,x)-((abs(a[0]-b[0]))*hightWat)),4))

    elif h <= hightWat and h+r > hightWat:
        a = [np.sqrt(-((P[1]-hightWat)**2 - r**2)) + P[0], hightWat]
        b = [2*P[0]-a[0], hightWat]

        x = np.arange(b[0],a[0],0.01)
        y = np.array(np.sqrt(r**2-(x-P[0])**2)+P[1])
        V_sub.append(round((np.pi*r**2)-((np.trapz(y,x)-((abs(a[0]-b[0]))*hightWat))),4))
        
    else:
        V_sub.append(np.pi*r**2)


    #print(V_sub[i]/V_k)
    ######### F_g - F_A #########

    m_w = roh_w * V_sub[i] #Masse des verdrängten Wassers

    F_res = (m_k - m_w) * g
    a_res = (F_res / m_k)
    #print(F_res)
    #############################
    v += a_res*(1/fps)
    #print(v)
    if V_sub[len(V_sub)-1] > 0:
        v *= 0.97

    E_kin = 1/2 * m_k * v**2
    E_pot = abs(m_k * g * h)

    s.append(s[len(s)-1]+v*1/fps)
    h += (s[len(s)-1]-s[len(s)-2])
    
    circ = plt.Circle((P[0],P[1]), r, color="r")
    ax[0].add_patch(circ)
    ax[0].plot([0,widthPot],[hightWat,hightWat], "blue")  #Wasser

    ax[1].bar(x_tit,[E_pot/1000,E_kin/1000,(E_pot + E_kin)/1000],color="C0")

    camera.snap()
    hightWat += (V_sub[len(V_sub)-1]-V_sub[len(V_sub)-2])/widthPot
    t += 1/fps


ax[0].plot([0,0],[0,hightPot],"black")                #Wand Links
ax[0].plot([widthPot,widthPot],[0,hightPot],"black")  #Wand Rechts
ax[0].plot([0,widthPot],[0,0],"black")                #Boden


animate = camera.animate(interval=1/fps,repeat=True,repeat_delay=500)
animate.save("C:/Users/moelj/downloads/FreierFallInsWasser.gif")
plt.show()






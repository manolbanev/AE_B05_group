import numpy as np
import matplotlib.pyplot as plt



#contsants 
W_pl =  8531 # kg payload for MTOW
W_S = 6201.56
AR = 7
n_ult = 3.75
M_dd = 0.964
M_crux = 0.935
lamda_4 = 0.756216 #rad
taper_ratio = 0.248757
t_c = 0.12
Tto = 85998 #lbf for two engines. Assume not changing the engine

#converting functions
def lbs_to_kg(lbs):
    return lbs*0.45359237

def kg_to_lbs(kg):
    return kg/0.45359237

def sq_m_to_sq_ft(S):
    return S*10.7639104


#arrays
MTOW_arr= [106227.1] #kg
OEW = [] #kg 56693.42
W_f_used = [] #kg fuel used
W_tfo = [] #kg trapped fuel and oils
W_f_res = [] #kg reserve fuel
W_f = [] # total fuel weight 17741.6
Ww = [] #kg wing
W_emp = [] #kg empanage
Wf_s= [] # kg fuel system weight
W_fc = [] # kg flight control system weight
W_iae = [] # kg instrumentation system weight
W_fur = [] # kg furnishing system weight
W_hs = [] # kg hydraulic system weight
W_apu = [] #kg apu system weight
W_paint = [] #kg paint weight
W_g_nose = [] #kg landing gear nose
W_g_main = [] #kg landing gear main
W_struc = [] # kg structural weight
W_p = [] #kg powerplant weight
W_feq = [] # kg fixed equipment weight


#Structure weight estimation

## Fuselage
W_fuselage = 15748 # kg


## Wing

def wing_weight_estimation(S,ZFW):
    cos_lamda_2 = 0.76289
    C_r = (2*S)/((1+taper_ratio)*np.sqrt(AR*S))
    t_r = t_c*C_r
    b = np.sqrt(S * AR) #ft
    Ww_new = 0.0017*ZFW*(b/cos_lamda_2)**(0.75)*(1+np.sqrt((6.3*cos_lamda_2)/b))*n_ult**(0.55)*((b*S)/(t_r*ZFW*cos_lamda_2))**(0.3)
    Ww.append(lbs_to_kg(Ww_new))


## Landing Gear

def landing_gear_weight_estimation(MTOW):
    k_g = 1.00
    A_g = [40,20]
    B_g = [0.16,0.10]
    C_g = [0.019,0]
    D_g=[1.5*10**(-5),2*10**(-6)]
    W_g_new_main = k_g*(A_g[0]+B_g[0]*MTOW**(3/4)+C_g[0]*MTOW+D_g[0]*MTOW**(3/2))
    W_g_new_nose = k_g*(A_g[1]+B_g[1]*MTOW**(3/4)+C_g[1]*MTOW+D_g[1]*MTOW**(3/2))
    W_g_nose.append(lbs_to_kg(W_g_new_nose))
    W_g_main.append(lbs_to_kg(W_g_new_main))

## Nacelle 
W_n=lbs_to_kg(0.065*Tto)

## Empennage 
def empennage_weight_estiamtion(S):
    V_d = 625.9
    K_h = 1.1
    K_v =1
    S_h = 0.1686*S+13.242
    S_v = 0.0922*S + 13.847
    cos_lamda_2_h = 0.7881
    cos_lamda_2_v = 0.65763
    W_h = K_h*S_h*((3.81*S_h**(0.2)*V_d)/(1000*np.sqrt(cos_lamda_2_h))-0.287)
    W_v = K_v*S_v*((3.81*S_v**(0.2)*V_d)/(1000*np.sqrt(cos_lamda_2_v))-0.287)
    W_emp_new = W_h+ W_v
    W_emp.append(lbs_to_kg(W_emp_new))

# Powerplant Weight Estimation

## Engines
W_eng = 6622 #kg both engines 

## Air Induction System 
None

## Fuel System 
def fuel_system_weight_estimation(Wf):
    N_e = 2
    N_t = 3
    K_fsp = 5.87
    Wfs_new = 80*(N_e +N_t - 1)+15*(N_t)**(0.5)*(Wf/K_fsp)**(0.333)
    Wf_s.append(lbs_to_kg(Wfs_new))

## Propulsion System
W_apsi = 148 #kg constant
W_tr = 0.18*W_eng #kg

#Fixed Equipment Weight Estimation

## Flight Control System
def flight_control_system(MTOW):
    W_fc_new = 0.64*MTOW**(2/3)
    W_fc.append(lbs_to_kg(W_fc_new))

## Instrumentation system
def instrumentation_system_weight_estiamtion(W_e):
    W_iae_new = 0.575*W_e**(0.556)*6326.13391**(0.25)
    W_iae.append(lbs_to_kg(W_iae_new))

## Furnishing system
def furnishing_system_weight(MTOW,Wf):
    W_fur_new = 0.211*(MTOW - Wf)**(0.91)
    W_fur.append(lbs_to_kg(W_fur_new))

## Hydraulic system 
def hydraulic_system_weight_estimation(MTOW):
    W_hs_new = 0.009*MTOW
    W_hs.append(lbs_to_kg(W_hs_new))

## APU system
def apu_system_weight_estimation(MTOW):
    W_apu_new =0.0085* MTOW
    W_apu.append(lbs_to_kg(W_apu_new))

## Paint
def paint_weight_estimation(MTOW):
    W_paint_new = 0.0045* MTOW
    W_paint.append(lbs_to_kg(W_paint_new))

## Air-conditioning, pressurization, anti- and de-icing system
W_api = 1444.461641 # kg

## Oxygen System
W_ox = 233.6907892 #kg

## Electrical System
W_els = 1494 #kg


# Fuel Weight Estimation

## Used Fuel
def used_fuel_weight_estimations(MTOW):
    c_j = 1.64854*10**(-5) #kg/s.N
    V = 257.64 #m/s
    D_L = 1/14.76
    R = 11716000 #m
    g = 9.81 #m/s^2
    M_ff = 1/(np.exp((R*g*D_L*c_j)/V))
    print(M_ff)
    W_f_used_new = MTOW*(1-M_ff)
    W_f_used.append(lbs_to_kg(W_f_used_new))

## Trapped Fuel 
def trapped_fuel_weight_estimation(MTOW):
    W_tfo_new = 0.005*MTOW
    W_tfo.append(lbs_to_kg(W_tfo_new))

## Reserve Fuel 
def reserve_fuel_weight_estimation(W_f_used):
    W_f_res_new = 0.25*W_f_used
    W_f_res.append(lbs_to_kg(W_f_res_new))

def main(MTOW_arr):
    for i in range(30):

        #W_fuel total
        used_fuel_weight_estimations(kg_to_lbs(MTOW_arr[i]))
        trapped_fuel_weight_estimation(kg_to_lbs(MTOW_arr[i]))
        reserve_fuel_weight_estimation(kg_to_lbs(W_f_used[i]))
        W_f.append(W_f_used[i]+W_f_res[i]+W_tfo[i])

        #Wing
        S = MTOW_arr[i]*9.81/W_S
        ZFW = MTOW_arr[i] - W_f[i]
        wing_weight_estimation(sq_m_to_sq_ft(S),kg_to_lbs(ZFW))

        #Landing Gear
        landing_gear_weight_estimation(kg_to_lbs(MTOW_arr[i]))

        #Empennage
        empennage_weight_estiamtion(sq_m_to_sq_ft(S))

        #Fuel System
        fuel_system_weight_estimation(kg_to_lbs(W_f[i]))

        #Flight Control System
        flight_control_system(kg_to_lbs(MTOW_arr[i]))

        #Instrumentation system
        W_e = MTOW_arr[i]-W_f[i]-W_pl
        instrumentation_system_weight_estiamtion(kg_to_lbs(W_e))

        #Furnishing system
        furnishing_system_weight(kg_to_lbs(MTOW_arr[i]),kg_to_lbs(W_f[i]))

        #Hydraulic system 
        hydraulic_system_weight_estimation(kg_to_lbs(MTOW_arr[i]))

        #APU system
        apu_system_weight_estimation(kg_to_lbs(MTOW_arr[i]))

        #Paint
        paint_weight_estimation(kg_to_lbs(MTOW_arr[i]))

        #Structural Weight sumation
        W_struc.append(W_fuselage+Ww[i]+W_g_main[i]+W_g_nose[i]+W_n+W_emp[i]) #kg

        #Powerplant Weight sumation
        W_p.append(W_eng+Wf_s[i]+W_apsi+W_tr) #kg

        #Fixed Equipment sumation
        W_feq.append(W_fc[i]+W_iae[i]+ W_fur[i]+ W_api+ W_ox+ W_els+ W_hs[i]+W_apu[i]+W_paint[i])

        #OEW sumation
        OEW.append(W_p[i]+W_feq[i]+W_struc[i])

        #MTOW sumation
        MTOW_arr.append(OEW[i]+ W_f[i]+W_pl)

main(MTOW_arr)
plt.plot(MTOW_arr,'o',label='MTOW')
plt.plot(OEW,'o',label='OEW')
plt.plot(W_f,'o',label = 'Fuel weight')
plt.xlabel('Recursive Cycle')
plt.ylabel('[kg]')
plt.title('Iteration Cycle')
plt.legend(loc= 'upper left')
plt.show()  



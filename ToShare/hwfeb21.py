import numpy as np
'''
e1 = 0.7 #Vm^-1
e2 = 3.2 #Vm^-1
phi1 = 20*(np.pi/180)
phi2 = 37*(np.pi/180)
chi = 1.0208
beta = 0.94073
e0 = e1*np.cos(phi1)/(np.cos(beta)*np.cos(chi))
i = e1**2 + e2**2
v = 2*e1*e2*np.sin(phi1-phi2)
q = e1**2 - e2**2
u = 2*e1*e2*np.cos(phi1-phi2)
polar = np.sqrt(q**2 + u**2 + v**2)/i
i = 513.6
q = 117.2
u = 227.8
v = 445.3
polar = np.sqrt(q**2 + u**2 + v**2)/i


import matplotlib.pyplot as plt
from scipy import optimize
def solve_keplers_equation(M, e, tol=1e-6):
    """Solve Kepler's equation M = E - e*sin(E) for the eccentric anomaly E using Newton's method."""
    E = M  # Initial guess
    while True:
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def compute_orbit(a, e, omega, M0, num_steps=1000, years=60):
    """Compute the orbit over a given period in years with num_steps time intervals."""
    
    GMs = 4 * np.pi ** 2  # GM_sun in AU^3/yr^2
    n = np.sqrt(GMs / a**3)  # Mean motion
    
    times = np.linspace(0, years, num_steps)
    rarr, Larr, nuarr,Marr = [], [], [], []
    #T = -M0/n
    for t in times:
        M = M0 + n*t 
        #M = (n*(t-T)) % 2*np.pi # Mean anomaly
        E = solve_keplers_equation(M, e) #% 2*np.pi # Eccentric anomaly
        r = a*(1-(e*np.cos(E)))  # Radius
        
        # Compute true anomaly (nu)
        nu = 2*np.arctan((((1+e)/(1-e))**(1/2))*np.tan(E/2)) #% 2*np.pi
        #nu = 2*np.arctan(tan_nu_2) #% 2*np.pi
        
        # Compute true longitude L = omega + nu
        L = (omega + nu) #% 2*np.pi
        
        # Convert to Cartesian coordinates
        x = r * np.cos(L)
        y = r * np.sin(L)
        
        rarr.append(r)
        Larr.append(L)
        nuarr.append(nu)
        Marr.append(M)
    
    return times, rarr, Larr, nuarr, Marr

def plot_orbit(xs, ys):
    """Plot the computed orbit."""
    plt.figure(figsize=(6,6))
    plt.plot(xs, ys, label='Orbit')
    plt.scatter([0], [0], color='yellow', label='Sun', s=100)
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.title("Orbital Path of Planet")
    plt.legend()
    plt.grid()
    plt.axis('equal')
    plt.show()


# Example usage
a = 1.0  # Semi-major axis in AU
e = 0.2  # Eccentricity
omega = np.radians(45)  # Longitude of pericenter in radians
M0 = np.radians(-10)  # Initial mean anomaly in radians

#planet values
a_sat = 9.53707032; e_sat = 0.05415060; om_sat = 92.43194*np.pi/180; lam_sat = 49.94432*np.pi/180
a_ur = 19.19126393; e_ur = 0.04716771; om_ur = 170.96424*np.pi/180; lam_ur = 313.23218*np.pi/180
a_nep = 30.06896348; e_nep = 0.00858587; om_nep = 44.97135*np.pi/180 ; lam_nep = 304.88003*np.pi/180

times_ur, rarr_ur, Larr_ur, nuarr_ur, Marr_ur = compute_orbit(a_ur, e_ur, om_ur, lam_ur - om_ur)
times_sat, rarr_sat, Larr_sat, nuarr_sat, Marr_sat= compute_orbit(a_sat, e_sat, om_sat, lam_sat - om_sat)
#import pdb; pdb.set_trace()
Larr = np.array(Larr_sat)*180/np.pi % 360
f_M = (np.array(nuarr_sat) % (2*np.pi) -np.array(Marr_sat) % (2*np.pi))  #% 2*np.pi
#plt.plot(times,Larr)
plt.plot(times_sat,f_M*180/np.pi) # % 2*np.pi)
plt.xlabel("t [years]")
plt.ylabel("f-M [deg]")
plt.title("Uranus: f-M vs. time")
import pdb; pdb.set_trace()
plt.show()
#plot_orbit(xs_sat, ys_sat)


j2 = 1e-3
a_R = 1.5
e = 0.01
inc = np.linspace(0,np.pi,100)
bigomegan = (-3/2)*j2*(a_R)**(-2)*np.cos(inc)/(1-e**2) 
lilomegan = 3*j2*(a_R)**(-2)*(1-((5/4)*(np.sin(inc))**2))/((1-e**2)**2)
mn = (3/2)*j2*(a_R)**(-2)*(1-((3/2)*(np.sin(inc))**2))/((1-e**2)**(3/2))
plt.plot(inc,bigomegan)
plt.xlabel("inc [rad]")
plt.ylabel("Big Omega Dot/n")
plt.show()
plt.close()
plt.plot(inc,lilomegan)
plt.xlabel("inc [rad]")
plt.ylabel("Little Omega Dot/n")
plt.show()
plt.close()
plt.plot(inc,mn)
plt.xlabel("inc [rad]")
plt.ylabel("M Dot/n - 1")
plt.show()
'''

M_sun = 1.989e30 #kg
P_moon = 27*86400 #s
P_earth = 365*86400 #s
omega_earth = 2*np.pi/(24*3600)
i_earth = 23.5 #degrees
i_moon = 5 #degrees
J2_earth = 1083e-6 
J2_sat = 0.01656
J2_ur = 0.003343
r_e = 6378137 #m
r_sat = 58232000 #m
r_ur =  25559*1000 #m
g = 6.67e-11 #m^3 kg^-1 s^-2 
a_moon = 3.84e8 #m 
a_earth = 1.496e+11 #m
a_sat = 9.58*a_earth #m
a_ur = 19.19*a_earth #m
M_earth = 5.972e+24 #kg 
M_sat = 5.683e26 #kg 
M_ur = 8.681e25 #kg
M_moon = 7.35e22 #kg 
m_mimas = 0.38 * 10**20 #kg
a_mimas = 185.52 * 10**6 #m 
m_encel = 0.65 *10**20
a_encel = 238.02  * 10**6
m_miranda = 0.659 *10**20
a_miranda = 129.8  * 10**6
m_ariel = 13.53*10**20
a_ariel = 191.2 * 10**6
a_tethys = 294.66 * 10**6
m_tethys = 6.27 * 10**20
a_rhea = 527.04 * 10**6
m_rhea = 23.1 * 10**20
domdt_j2 = (-3/2)*(2*np.pi/P_moon)*J2_earth*(r_e/a_moon)**2 * np.cos(i_earth*np.pi/180)
domdt_sun = (-3/4)*np.cos(i_moon*np.pi/180)*(2*np.pi/P_earth)**2/(2*np.pi/P_moon)
ac_earth = (((2*J2_earth*r_e**2 * a_earth**3 *M_earth)/(M_sun))**(1/5))/r_e 
ac_sat = (((2*J2_sat*r_sat**2 * a_sat**3 *M_sat)/(M_sun))**(1/5))/r_sat
ac_ur = (((2*J2_ur*r_ur**2 * a_ur**3 *M_ur)/(M_sun))**(1/5))/r_ur
k2 = 0.4
roche_s = 2.5*r_sat
roche_u = 2.5*r_ur
delt = 4.5e9*3.154e+7
q_satmimas = 1/((2*np.sqrt(M_sat))*(a_mimas**(13/2) - roche_s**(13/2))*(np.sqrt(g)*13*3*k2*m_mimas*delt*r_sat**5)**(-1))
q_satencel =  1/((2*np.sqrt(M_sat))*(a_encel**(13/2) - roche_s**(13/2))*(np.sqrt(g)*13*3*k2*m_encel*delt*r_sat**5)**(-1))
#print(q_satmimas)
#print(q_satencel)
q_urmir = 1/((2*np.sqrt(M_ur))*(a_miranda**(13/2) - roche_u**(13/2))*(np.sqrt(g)*13*3*k2*m_miranda*delt*r_ur**5)**(-1))
q_urar = 1/((2*np.sqrt(M_ur))*(a_ariel**(13/2) - roche_u**(13/2))*(np.sqrt(g)*13*3*k2*m_ariel*delt*r_ur**5)**(-1))
#print(q_urmir)
#print(q_urar)
q_satmeasured = 0.38/(1.6e-4)
#print(q_satmeasured)
t_elapsed_mimas = (a_mimas**(13/2) - roche_s**(13/2))*(2/13)*((3*(1.6e-4)*(g**(1/2))*m_mimas*r_sat**5)/(M_sat**(1/2)))**(-1)
t_elapsed_tethys = (a_tethys**(13/2) - roche_s**(13/2))*(2/13)*((3*(1.6e-4)*(g**(1/2))*m_tethys*r_sat**5)/(M_sat**(1/2)))**(-1)
t_elapsed_rhea = (a_rhea**(13/2) - roche_s**(13/2))*(2/13)*((3*(1.6e-4)*(g**(1/2))*m_rhea*r_sat**5)/(M_sat**(1/2)))**(-1)
#print(t_elapsed_mimas/(1e9*1e7*np.pi))
#print(t_elapsed_tethys/(1e9*1e7*np.pi))
#print(t_elapsed_rhea/(1e9*1e7*np.pi))
#import pdb; pdb.set_trace()

#define beta 
k2 = 0.29
qmoon = 12
qearth = 100
iearth = 0.331*M_earth*r_e**2 
beta = (3*k2*(g**(1/2))*M_moon*(r_e**5))/(qmoon*(M_earth**(1/2)))
#evolution of a 
#timearr = np.linspace(-5e9*3.154e+7,0,1000)
timearr = np.linspace(0, 5e9*3.154e+7,1000)
#a_evol = ((a_moon)**(13/2) + (13/2)*beta)**(2/13)#*time^2/13
tidaltorque = ((3/2)*k2*g*(M_sun**2)*(r_e**5)/(qearth*a_earth**6))
import matplotlib.pyplot as plt
a_arr = []
omega_arr = []
n_arr = []
for i in range(len(timearr)):
    a_evol = ((a_moon)**(13/2) + (13/2)*beta*timearr[i])**(2/13)
    n_evol = np.sqrt((4*(np.pi**2)*a_evol**3)/(g*M_earth))
    #omega_evol = (((3/2)*k2*g*(M_sun**2)*(r_e**5)/(qmoon*a_evol**6))*timearr[i])/
    

    n_arr.append(n_evol)
    a_arr.append(a_evol) 

#a_evol_syn = ((a_moon)**(13/2) + (13/2)*beta*6.75e9**3.154e+7)**(2/13)
import pdb; pdb.set_trace()

#plt.plot(timearr, tidaltorque, label = 'torque')
plt.plot(timearr/(1e9*3.154e+7),np.array(a_arr)/r_e, label = 'semimajor axis')
plt.xlabel("Time (gigayears)")
plt.ylabel("Semimajor axis of Moon (R_e)")
plt.legend()
plt.title("Evolution of Moon's Semimajor axis over time (going forward)")
plt.show()
plt.plot(timearr/(1e9*3.154e+7),np.array(n_arr)/3600, label = '2pi/n')
plt.xlabel("Time (gigayears)")
plt.ylabel("Lunar Orbital Period (hours)")
plt.title("Evolution of Lunar Period (2pi/n) going forward")
plt.legend()
plt.show()
#import pdb; pdb.set_trace()

e_moon = 0.0549
c_earth  = 0.331
Lorbmoon = M_moon*np.sqrt(g*(M_earth)*a_moon*(1-e_moon**2))
Lspinearth = (2*np.pi/(24*3600))*c_earth*M_earth*r_e**2
Lorbearth = M_earth*np.sqrt(g*(M_sun)*a_earth*(1-0.01671**2))

Lspinearth_new = (2*np.pi/(3.154e+7))*c_earth*M_earth*r_e**2
Lorbearth_new = Lspinearth_new/(Lspinearth/Lorbearth)
a_earth_new = (Lorbearth_new/(M_earth*np.sqrt(g*M_sun*(1-0.01671**2))))**2
delta_a = a_earth - a_earth_new

print(delta_a)
#print(Lspinearth/Lorbmoon)
#print(Lspinearth/Lorbearth)
'''
omega = np.linspace(0.1,4.1,100)
omega0 =1
beta = 0.1
A_omega = 1/np.sqrt((omega0-omega**2)**2 + (beta*omega)**2)
import matplotlib.pyplot as plt
plt.plot(omega,A_omega)
plt.xlabel("omega")
plt.ylabel("A(omega)")
plt.show()
import pdb; pdb.set_trace()
epsilon = beta*omega/(omega0**2 - omega**2)
plt.plot(omega, np.arctan2(beta * omega, omega0**2 - omega**2))
plt.xlabel("omega")
plt.ylabel("epsilon")
plt.show()
import pdb; pdb.set_trace()'''
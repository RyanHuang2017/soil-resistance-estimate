from archimede import compK
import numpy as np
import matplotlib.pyplot as plt

phi = 25*np.pi/180.0  # The friction angle of the Ottawa F65 sample
# R = 0.02                # unit:m, the diameter of the cone and cylinder
theta = 30*np.pi/180    # half apex angle of the cone
thos = 1537             # packing density of the sands

# hcone = R/np.tan(theta) # height of the cone
# S = np.pi*R**2          # cross section area of the cone and cylinder
k_phi = compK(phi)      # k_phi value from (Kang et al.2018)
g = 9.81                # gravitational acceleration m/s2


Radii = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
Depth = [0.10, 0.10, 0.10, 0.10, 0.10, 0.10]
Hcone = Radii/np.tan(theta)
Sii = [np.pi*x**2 for x in Radii]
V0 = Radii * Hcone * (1.0/3.0)
# print(V0)
V0 = V0 + Sii * (Depth - Hcone)
# print(V0)
Fz = k_phi*thos*g*V0
print(Fz)
fig, ax = plt.subplots(figsize=(3.0,3.0))
plt.plot(Radii,Fz/1000, '-x')
ax.set_xlim([0,0.05])
ax.set_ylim([0,0.5])
ax.set_xlabel('penetrator radii (m)')
ax.set_ylabel('penetration resistance (kN)')
fig.subplots_adjust(bottom=0.18, top=0.95, left=0.18, right=0.92, wspace=0.25, hspace=0.05)
plt.savefig('scaling law.png', dpi=600)
plt.show()
# ================================
# Side_and_brake_force_calculation
# ================================

# import modules
import argparse
import re
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

# value = input("Vehicle mass in kg: ")

# set up argument parser
parser = argparse.ArgumentParser(description="Parse variables")

with open("Readme.txt", "r", encoding="utf-8") as textfile:
    for line in textfile:
        veh_mass_match = re.search(r"Vehicle mass=(\d+)kg", line)
        slip_match = re.search(r"Slip=(\d+)°", line)
        camber_match = re.search(r"Camber=(\d+)°", line)
        my_match = re.search(r"my=(\d+)", line)

        # vehicle mass in kg
        if veh_mass_match:
            veh_mass = float(veh_mass_match.group(1))
            parser.add_argument("--veh_mass", type=float,
                                default=veh_mass, help="vehicle mass in kg")

        # slip angle in degrees
        if slip_match:
            slip = float(slip_match.group(1))
            parser.add_argument("--alpha", type=float,
                                default=slip, help="slip angle in degrees")

        # camber in degrees
        if camber_match:
            camber = float(camber_match.group(1))
            parser.add_argument("--gamma", type=float,
                                default=camber, help="camber in degrees")

        # friction coefficient
        if my_match:
            my = float(my_match.group(1))
            parser.add_argument("--my", type=float,
                                default=my, help="friction coefficient")


# vertical loads in kN
parser.add_argument("--v_load", type=float, default=[2,4,6,8],
                    help="vertical loads in kN")

# parse the variables
args=parser.parse_args()

Fz = args.veh_mass*const.g/4000     # vertical load of one wheel in kN
alpha = np.radians(args.alpha)      # slip angle in radians
gamma = np.radians(args.gamma)      # camber in radians
kappa = np.linspace(0,1,1001)       # longitudinal slip

print("veh_mass = " + str(args.veh_mass) + "kg")
print("Vertical load Fz = " + str(round(Fz*1000)) + "N")
print("Slip = " + str(args.alpha) + "°")
print("Camber = " + str(args.gamma) + "°")
print("my = " + str(args.my))

# Constant coefficients
A1_FX = -21.3
A2_FX = 1144
A3_FX = 49.6
A4_FX = 226
A5_FX = 0.069
A6_FX = -0.006
A7_FX = 0.056
A8_FX = 0.486

A1_FY = -22.1
A2_FY = 1011
A3_FY = 1078
A4_FY = 1.82
A5_FY = 0.208
A6_FY = 0
A7_FY = -0.354
A8_FY = 0.707
A9_FY = 0.028
A10_FY = 0
A11_FY = 14.8
A12_FY = 0.022
A13_FY = 0

for Fz in args.v_load:
    # ==========
    # Side force
    # ==========

    # coefficients for side force
    D_side = A1_FY*Fz**2 + A2_FY*Fz             # peak factor
    C_SIDE = 1.30                               # shape factor
    B_side = ((A3_FY*np.sin(A4_FY*np.arctan(A5_FY*Fz))) /
              (C_SIDE*D_side)) * (1-A12_FY*np.abs(gamma))     # stiffness factor
    E_side = A6_FY*Fz**2 + A7_FY*Fz + A8_FY     # curvature factor

    # horizontal and vertikal shift
    d_shift_h = A9_FY * gamma
    d_shift_v = (A10_FY*Fz**2 + A11_FY*Fz) * gamma

    # Phi
    phi_side = ((1-E_side)*(np.degrees(alpha)+d_shift_h) +
                (E_side/B_side)*np.arctan(B_side*(alpha+d_shift_h)))

    # Side force
    Fy0 = D_side*np.sin(C_SIDE*np.arctan(B_side*phi_side))+d_shift_v

    # ===========
    # Brake force
    # ===========

    # coefficients for brake force
    D_brake = A1_FX * Fz**2 + A2_FX * Fz            # peak factor
    C_BRAKE = 1.65                                  # shape factor
    B_brake = (A3_FX*Fz**2 + A4_FX*Fz)/(C_BRAKE*D_brake*np.e**(A5_FX*Fz))       # stiffness factor
    E_brake = A6_FX * Fz**2 +A7_FX * Fz + A8_FX     # curvature factor

    # Phi
    phi_brake = (1-E_brake)*kappa + (E_brake/B_brake)*np.arctan(B_brake*kappa)*180/np.pi

    # Brake force
    Fx0 = D_brake*np.sin(C_BRAKE*np.arctan(B_brake*phi_brake))

    # ===============
    # Brake and side force for combined cornering and braking
    # ===============

    # calculate Sigma
    sigma_x = -kappa/(1+kappa)
    sigma_y = -np.tan(alpha)/(1+kappa)

    sigma = np.sqrt(sigma_x**2 + sigma_y**2)

    # Brake force
    Fx = -(sigma_x/sigma)*Fx0
    # Side force
    Fy = -(sigma_y/sigma)*Fy0

    plt.plot(kappa*100,Fy,label="Fy at " + str(Fz*1000) +" N")
    plt.plot(kappa*100,Fx,label="Fx at " + str(Fz*1000) + " N")
    plt.legend(loc="best")

plt.title("Side and brake force")
plt.xlabel("Longitudinal slip " + r"$\kappa$" + " [%]")
plt.ylabel("Side force / Brake force [N]")
plt.show()

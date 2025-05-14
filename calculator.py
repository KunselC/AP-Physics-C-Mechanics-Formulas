# AP Physics C: Mechanics Formulas for TI-84 Plus CE Python
# Note: Input values should be in standard SI units (meters, kilograms, seconds, etc.)

import math

# --- Constants ---
G = 6.674e-11  # Gravitational constant (N*m^2/kg^2)
g = 9.81       # Acceleration due to gravity near Earth (m/s^2) - can be adjusted

# --- Kinematics ---
def kinematics():
    print("--- Kinematics ---\nCalculus Forms:\n  r = r0 + integral(v dt)\n  v = v0 + integral(a dt)")
    input("Press Enter for more...")
    print("Constant Acceleration:\n  dx = (v + v0)/2 * t\nVector Forms (Constant Accel):\n  v_vec = v0_vec + a_vec*t")
    input("Press Enter for more...")
    print("  dr_vec = v0_vec*t + 0.5*a_vec*t**2\nProjectile Motion (level ground, launch angle theta):\n  Max Height (y_max): (v0*sin(theta))**2 / (2*g)\n  Range (R): (v0**2 * sin(2*theta)) / g")

# --- Dynamics ---
def dynamics():
    print("--- Dynamics ---\nNewton\'s Second Law:\n  F_net = m*a\n  F_net_vec = m*a_vec")
    input("Press Enter for more...")
    print("Friction:\n  Friction (kinetic): fk = mu_k * N\nDrag Force (examples):\n  F_drag = -b*v")
    input("Press Enter for more...")
    print("  F_drag = -c*v**2\nTerminal Velocity (for F_drag = -bv): v_t = mg/b")

# --- Work, Energy, Power ---
def work_energy_power():
    print("--- Work, Energy, Power ---\nWork:\n  Work (constant force): W = F * d * cos(theta)\n  W = F_vec . d_vec")
    input("Press Enter for more...")
    print("Potential Energy:\n  Work done by spring: Ws = -delta_Us\n  Ws = -(0.5*k*xf**2 - 0.5*k*xi**2)")
    input("Press Enter for more...")
    print("Conservation of Energy:\n  E_mech = K + U\n  delta_E = W_nc (Work by non-conservative forces)\n  E_initial = E_final (if W_nc = 0)")
    input("Press Enter for more...")
    print("Power:\n  P = F . v")

# --- Momentum & Collisions ---
def momentum():
    print("--- Momentum & Collisions ---\nConservation of Momentum (isolated system):\n  p_initial_total_vec = p_final_total_vec\n  sum(m_i * v_i_initial) = sum(m_i * v_i_final)")
    input("Press Enter for more...")
    print("Center of Mass:\n  a_cm_vec = sum(m_i * a_i_vec) / sum(m_i) = F_net_ext_vec / M_total\n1D Elastic Collisions (Special Case):\n  v1f = ((m1-m2)/(m1+m2))*v1i + ((2*m2)/(m1+m2))*v2i")
    input("Press Enter for more...")
    print("  v2f = ((2*m1)/(m1+m2))*v1i + ((m2-m1)/(m1+m2))*v2i\nCoefficient of Restitution (e): |v2f - v1f| / |v1i - v2i|\n  (e=1: elastic, 0<=e<1: inelastic, e=0: perfectly inelastic)")

# --- Circular Motion & Rotation ---
def rotation():
    print("--- Circular Motion & Rotation ---\nUniform Circular Motion:\n  Centripetal Force: Fc = m * ac = m*v**2 / r\nAngular Variables:\n  theta = s / r")
    input("Press Enter for more...")
    print("Relationship between Linear and Angular (vectors):\n  a_rad_vec = omega_vec x v_vec (Radial/Centripetal acceleration)\nAngular Kinematics (Constant Angular Accel):\n  dtheta = (omega + omega0)/2 * t")
    input("Press Enter for more...")
    print("Torque:\n  Magnitude: tau = r * F * sin(theta)\nCommon Shapes (axis through CM):\n    Hoop: I = MR^2\n    Disk/Cylinder: I = 0.5*MR^2")
    input("Press Enter for more...")
    print("    Rod (center): I = (1/12)ML^2\n    Rod (end): I = (1/3)ML^2\n    Solid Sphere: I = (2/5)MR^2\n    Hollow Sphere: I = (2/3)MR^2")
    input("Press Enter for more...")
    print("Rolling Motion (No Slipping):\n  a_cm = R * alpha\n  K_total = K_trans + K_rot = 0.5*m*v_cm**2 + 0.5*I_cm*omega**2\nRotational Power: P = tau * omega")
    input("Press Enter for more...")
    print("Angular Momentum:\n  Magnitude (point particle): L = mvr*sin(theta)\nConservation of Angular Momentum:\n  L_initial_total_vec = L_final_total_vec (if tau_net_ext = 0)")

# --- Oscillations (Simple Harmonic Motion) ---
def oscillations():
    print("--- Oscillations (SHM) ---\nDefining Equation: d^2x/dt^2 = -omega**2 * x\nSolution:\n  Velocity: v(t) = -A*omega * sin(omega*t + phi)")
    input("Press Enter for more...")
    print("  Acceleration: a(t) = -A*omega**2 * cos(omega*t + phi) = -omega**2 * x\nAngular Frequency:\n  General: omega = 2*pi*f = 2*pi/T\n  Torsional pendulum: omega = sqrt(kappa / I)")
    input("Press Enter for more...")
    print("Frequency: f = omega / (2*pi) = 1/T\nEnergy in SHM:\n  Total Energy: E = K + Us = 0.5 * k * A**2 = 0.5 * m * vmax**2")

# --- Gravitation ---
def gravitation():
    print("--- Gravitation ---\nNewton\'s Law of Gravitation:\n  Force Vector: Fg_vec = - (G * m1 * m2 / r**2) * r_hat\nGravitational Field (g):\n  g_vec = Fg_vec / m = - (G * M / r**2) * r_hat\n  Magnitude: g = G * M / r**2")
    input("Press Enter for more...")
    print("Orbital Mechanics (Circular Orbits):\n  Orbital Speed: v = sqrt(G*M / r)\n  Orbital Period: T = 2*pi*r / v = 2*pi * sqrt(r**3 / (G*M))")
    input("Press Enter for more...")
    print("Kepler\'s Laws:\n  1st: Planets move in elliptical orbits with the Sun at one focus.\n  2nd: A line connecting a planet to the Sun sweeps out equal areas in equal times (dL/dt = 0).")
    input("Press Enter for more...")
    print("  3rd: T^2 / a^3 = 4*pi**2 / (G*M) (a=semi-major axis, r for circular)\nTotal Energy in Circular Orbit: E = K + Ug = -G*M*m / (2*r)\nEscape Speed: v_esc = sqrt(2*G*M / R)")

# --- Main Menu ---
def main_menu():
    while True:
        print("""
--- AP Physics C: Mechanics Formulas ---
Select a topic:
1. Kinematics
2. Dynamics
3. Work, Energy, Power
4. Momentum & Collisions
5. Rotation
6. Oscillations (SHM)
7. Gravitation
0. Exit""")

        choice = input("Enter choice: ")

        if choice == '1':
            kinematics()
            input("Press Enter to return to menu...")
        elif choice == '2':
            dynamics()
            input("Press Enter to return to menu...")
        elif choice == '3':
            work_energy_power()
            input("Press Enter to return to menu...")
        elif choice == '4':
            momentum()
            input("Press Enter to return to menu...")
        elif choice == '5':
            rotation()
            input("Press Enter to return to menu...")
        elif choice == '6':
            oscillations()
            input("Press Enter to return to menu...")
        elif choice == '7':
            gravitation()
            input("Press Enter to return to menu...")
        elif choice == '0':
            print("Exiting.")
            break
        else:
            print("Invalid choice. Please try again.")

# --- Run the program ---
main_menu()
import numpy as np
import os
import matplotlib.pyplot as plt

#Aaron Stirk - Engineering Challenge - 6/4/2020
#Aknoledgements to JoshtheEngineer.com for executing xfoil commands using python
#code certainly could be optimized to run faster, however my main goal was to submit it in a timely manner

cptxt = 'Cp.txt'
ptxt = 'Polar.txt'
xfnm = 'xfinput.txt'

re = ['3e6','10e6','15e6']

fig, axs = plt.subplots(4, 1,)
for x in re:
    if os.path.exists(ptxt):
        os.remove(ptxt)
    if os.path.exists(xfnm):
        os.remove(xfnm)

    fid = open(xfnm, "w")
    fid.write("load " + 'naca633618.dat' + "\n")
    fid.write("pane\n")
    fid.write("oper\n")
    fid.write("visc\n")
    fid.write(x + "\n")
    fid.write("iter 1000\n")
    fid.write("pacc\n")
    fid.write(ptxt + "\n\n")
    fid.write("aseq\n")
    fid.write("-10\n")
    fid.write("20\n")
    fid.write("1\n")
    fid.write("\n\n\nquit\n")
    fid.close()

    os.system("xfoil.exe < xfinput.txt")

    if os.path.exists(xfnm):
        os.remove(xfnm)

    data1 = np.loadtxt(ptxt, skiprows=12)
    alfa = data1[:, 0]
    cl = data1[:, 1]
    cd = data1[:, 2]
    cm = data1[:, 4]
    clcd = cl / cd

    axs[0].plot(alfa, cl)
    axs[0].set_ylabel('Lift Coeff (Cl)')
    fig.suptitle("Aerodynamic Polar Coefficients Over Range of AoA", fontsize=14)

    axs[1].plot(alfa, cd)
    axs[1].set_ylabel('Drag Coeff (Cd)')

    axs[2].plot(alfa, cm)
    axs[2].set_ylabel('Moment Coeff (Cm)')

    axs[3].plot(alfa, clcd)
    axs[3].set_xlabel('Angle of Attack (deg)')
    axs[3].set_ylabel('Lift to Drag Ratio (Cl/Cd)')

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)

axs[0].legend(('Re = 3e6', 'Re = 10e6', 'Re = 15e6'))

fig, axs = plt.subplots(3, 1, sharex=True)
alf = ['0','4','8','12']
col = ['r', 'b', 'g', 'c']
for i in alf:

    if os.path.exists(xfnm):
        os.remove(xfnm)
    if os.path.exists('3e6' + cptxt):
        os.remove('3e6' + cptxt)
    if os.path.exists('10e6' + cptxt):
        os.remove('10e6' + cptxt)
    if os.path.exists('15e6' + cptxt):
        os.remove('15e6' + cptxt)

    if i == '0':
        color = col[0]
    if i == '4':
        color = col[1]
    if i == '8':
        color = col[2]
    if i == '12':
        color = col[3]

    for x in re:

        fid = open(xfnm, "w")
        fid.write("load " + 'naca633618.dat' + "\n")
        fid.write("pane\n")
        fid.write("oper\n")
        fid.write("visc\n")
        fid.write(x + "\n")
        fid.write("iter 1000\n")
        fid.write("alfa\n")
        fid.write(i + "\n")
        fid.write("cpwr "+ x + cptxt + "\n")
        fid.write("\n\nquit\n")
        fid.close()

        os.system("xfoil.exe < xfinput.txt")

        if os.path.exists(xfnm):
            os.remove(xfnm)

        db = np.loadtxt(x + cptxt, skiprows=3)

        X_0 = db[:, 0]
        Y_0 = db[:, 1]
        Cp_0 = db[:, 2]

        XU = X_0[0:80]
        XL = X_0[80:]
        YU = Y_0[0:80]
        YL = Y_0[80:]

        CPU = Cp_0[0:80]
        CPL = Cp_0[80:]

        if x == '3e6':
            axs[0].plot(XU, CPU, color + '-', label='AoA = ' + i)
            axs[0].plot(XL, CPL, color + '--')
            axs[0].set_ylabel('Pressure Coefficient (Cp)')
            axs[0].set_title('Pressure Coefficient Distributions for Re = ' + x)
            axs[0].set_xlim(-.01,1)

            fig.suptitle('Pressure Coefficients with AoA Plotted For Each Re', fontsize=14)
        elif x =='10e6':
            axs[1].plot(XU, CPU, color + '-', label='AoA = ' + i)
            axs[1].plot(XL, CPL, color + '--')
            axs[1].set_ylabel('Pressure Coefficient (Cp)')
            axs[1].set_title('Pressure Coefficient Distributions for Re = ' + x)

        elif x =='15e6':
            axs[2].plot(XU, CPU, color + '-', label='AoA = ' + i)
            axs[2].plot(XL, CPL, color + '--')
            axs[2].set_ylabel('Pressure Coefficient (Cp)')
            axs[2].set_xlabel('Position on Foil (x)')
            axs[2].set_title('Pressure Coefficient Distributions for Re = ' + x)

        if os.path.exists(x + 'Cp.txt'):
            os.remove(x + 'Cp.txt')

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.7,
                        wspace=0.35)
axs[0].invert_yaxis()
axs[1].invert_yaxis()
axs[2].invert_yaxis()
axs[0].legend()

fig, axs = plt.subplots(3, 1)
fig.suptitle('Boundary Layer Properties at Trailing Edge', fontsize=14)
for x in re:

    if os.path.exists(xfnm):
        os.remove(xfnm)

    a = np.arange(-10,21,1)
    hka = np.zeros(31)
    dta = np.zeros(31)
    theta = np.zeros(31)

    for i in a:

        if os.path.exists('dt.txt'):
            os.remove('dt.txt')
        if os.path.exists('hk.txt'):
            os.remove('hk.txt')

        fid = open(xfnm, "w")
        fid.write("load " + 'naca633618.dat' + "\n")
        fid.write("pane\n")
        fid.write("oper\n")
        fid.write("visc\n")
        fid.write(x + "\n")
        fid.write("iter 1000\n")
        fid.write("pacc\n\n\n")
        fid.write("a " + str(i) + "\n")
        fid.write("vplo\n")
        fid.write("hk\n")
        fid.write("dump hk.txt\n")
        fid.write("dt\n")
        fid.write("dump dt.txt\n")
        fid.write("\n\n\nquit\n")
        fid.close()

        os.system("xfoil.exe < xfinput.txt")

        if os.path.exists(xfnm):
            os.remove(xfnm)

        data1 = np.loadtxt('hk.txt', skiprows=7)
        data2 = np.loadtxt('dt.txt', skiprows=7)

        X = data1[:, 0]
        X2 = data2[:, 0]
        hkr = data1[:, 1]
        dtr = data2[:, 1]

        hk = hkr[X == 1]
        dt = dtr[X2 == 1]

        hka[i+10] = hk[0]

        if len(dt) > 1:
            dta[i+10] = dt[0]
        elif len(dt) == 1:
            dta[i+10] = dt

        theta[i+10] = dta[i+10]/hka[i+10]

        if os.path.exists('hk.txt'):
            os.remove('hk.txt')
        if os.path.exists('dt.txt'):
            os.remove('dt.txt')

    axs[0].plot(a, dta)
    axs[0].set_ylabel('Displacement Thickness (*delta)')

    axs[1].plot(a, theta)
    axs[1].set_ylabel('Moment Thickness (Theta)')

    axs[2].plot(a, hka)
    axs[2].set_ylabel('Shape Factor (Hk)')
    axs[2].set_xlabel('Angle of Attack (deg)')

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.7,
                        wspace=0.35)
axs[0].legend(('Re = 3e6', 'Re = 10e6', 'Re = 15e6'))
plt.show()
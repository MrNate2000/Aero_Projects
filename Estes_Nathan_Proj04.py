import numpy as np
from matplotlib import pyplot as plt

print("Quasi-1D Flow SOlution (MacCormack's Technique)")
print(" Supersonic, Nonconservative Form")
print(" Nozzle Shape: A/At = 1.00 + 2.20(x - 1.00); dx: 0.1")
print(" Courant #: 0.5")
print(' Tolerance: 1.00E-9')
for x in range(1,4) :
    
    xmin = 0
    xmax = 3
    gridpoints = 30 * x + 1
    dx = (xmax - xmin)/(gridpoints - 1)
    print("\n Gridpoints:", gridpoints)
    residual = 1
    n1 = 0

    #Nozzle Parameters, initial conditions
    Temp = []
    rho = []
    vel = []
    Area = []
    xval = []
    pres = [] 
    mdot = []
    mach = []
    iterate = []
    itercnt = []
    itercnt = np.append(itercnt, 0)
    x = 0
    cour = 0.5
    a,b,c = [1, 2.2, 1.5]

    def finddt(n1, nmxm1, dx, vel, Temp, cour):
        dtmin = 1
        for i in range(n1, nmxm1):
            dtmax = dx/(vel[i]+np.sqrt(Temp[i]))
            dtmin = min(dtmin, dtmax)
        dtmin = cour*dtmin
        return dtmin



    for i in range(gridpoints):
        xval = np.append(xval, round(x, 4))
        Temp = np.append(Temp, 1-.2314*x)
        rho = np.append(rho, 1 - 0.3146*x)
        vel = np.append(vel, (0.1 + 1.09*x)*np.sqrt(Temp[i]))
        pres = np.append(pres, rho[i]*Temp[i])
        Area= np.append(Area,  a+ b*(x - c)**2)
        mach = np.append(mach, vel[i]/np.sqrt(Temp[i]))
        mdot = np.append(mdot, rho[i]*vel[i]*Area[i])
        x += dx
        iterate = np.append(iterate, (i +1))
    print('\n----------------------------------------')
    print('\t\tInitial Conditions')
    print('----------------------------------------')
    print('i: ',[round(x, 4)for x in iterate])
    print('\nx: ', [round(x,4)for x in xval])
    print('\nArea:', [round(x,4)for x in Area] )
    print('\nrho: ',[round(x,4)for x in rho])
    print('\ntemp: ',[round(x,4)for x in Temp])
    print('\nvel: ',[round(x,4)for x in vel])

    # Calculate the timestep based on internal points

    dt = finddt(1, gridpoints-1, dx, vel, Temp, cour)

    # Integrate in time; use MacCormack's technique
    throat = int(np.argwhere(Area == min(Area)))

    nmax = gridpoints - 1

    pres[0] = 1
    it = 0
    it2 = 1
    tol = 0.000000001
    itmax = 50000
    pp=[]
    ppp=[]
    mm=[]
    mmm=[]
    pp.append([rho[i]*Temp[i] for i in range(0,gridpoints)])
    ppp.append(rho[throat]*Temp[throat])
    mm.append([vel[i]/np.sqrt(Temp[i]) for i in range(0,gridpoints)])
    mmm.append(vel[throat]/np.sqrt(Temp[throat]))
    
    while  it < itmax and residual > tol :
        prho = []
        ptemp = []
        pvel = []
        drdt = []
        dudt = []
        dTdt = []
        
        prho = np.append(prho, rho[0])
        ptemp = np.append(ptemp, Temp[0])
        pvel = np.append(pvel, 0)
    
        # Predictor step
        for i in range(1, nmax):
        # Forward Difference
          fdrho = (rho[i+1]-rho[i])/dx
          fdtemp = (Temp[i+1]-Temp[i])/dx
          fdvel = (vel[i+1]-vel[i])/dx
          fdArea = (np.log(Area[i+1])-np.log(Area[i]))/dx
        
        # Compute t derivatives
          drdt = np.append(drdt, -rho[i]*fdvel - rho[i]*vel[i]*fdArea - vel[i]*fdrho)
          dudt = np.append(dudt, -vel[i]*fdvel - (1/1.4)*(fdtemp+(Temp[i]/rho[i])*fdrho))
          dTdt = np.append(dTdt, -vel[i]*fdtemp-0.4*Temp[i]*(fdvel+vel[i]*fdArea))
        
        # Compute predicted variables
          prho = np.append(prho, rho[i] + drdt[i-1]*dt)
          ptemp = np.append(ptemp, Temp[i] + dTdt[i-1]*dt)
          pvel = np.append(pvel, vel[i] + dudt[i-1]*dt)
          
        # Extrapolate velocity to outflow
        pvel[0] = 2*pvel[1] - pvel[2]
        
        # Corrector Step
        for i in range(1, nmax):
        # Backward Difference
          bdrho = (prho[i]-prho[i-1])/dx
          bdtemp = (ptemp[i]-ptemp[i-1])/dx
          bdvel = (pvel[i]-pvel[i-1])/dx
          bdArea = (np.log(Area[i])-np.log(Area[i-1]))/dx
        
        # Compute t derivatives
          cdrdt = -prho[i]*bdvel - prho[i]*pvel[i]*bdArea - pvel[i]*bdrho
          cdudt = -pvel[i]*bdvel - (1/1.4)*(bdtemp+(ptemp[i]/prho[i])*bdrho)
          cdTdt = -pvel[i]*bdtemp-0.4*ptemp[i]*(bdvel+pvel[i]*bdArea)
        
        # Average t derivatives
          avgdrdt = 0.5*(drdt[i-1] + cdrdt)
          avgdTdt = (dTdt[i-1] + cdTdt)/2
          avgdudt = (dudt[i-1] + cdudt)/2
          
        # Update variables
          rho[i] =  rho[i] + avgdrdt*dt
          vel[i] =  vel[i] + avgdudt*dt
          Temp[i] = Temp[i] + avgdTdt*dt
          
        # Define residual at throat
          if i == throat:
              residual = abs(avgdrdt)
              ppp.append(rho[i] * Temp[i])
              mmm.append(vel[i] / np.sqrt(Temp[i]))
    
        # Apply boundary conditions grid point 1
        vel[0] = 2*vel[1] - vel[2]
    
        # Apply boundary conditions grid point N (float all variables)
        vel[nmax] = 2*vel[nmax-1] - vel[nmax-2]
        rho[nmax] = 2*rho[nmax-1] - rho[nmax-2]
        Temp[nmax] = 2*Temp[nmax-1] - Temp[nmax-2]
        
        
        dt=finddt(1,gridpoints-1,dx,vel,Temp,cour)
        itercnt = np.append(itercnt, it)
        it=it+1
        if it==25 or it==50 or it==100:
            pp.append([rho[i]*Temp[i] for i in range(0,gridpoints)])
            mm.append([vel[i]/np.sqrt(Temp[i]) for i in range(0,gridpoints)])
            
    # Exact Values
    def findmexactf(g,aa,m):

        f = ((1 / m) ** 2) * (((2 / (g + 1)) * (1 + ((g - 1) / 2) * m * m)) ** ((g + 1) / (g - 1))) - aa ** 2
        return f


    def findmexactsub(g,aa):
    
        f = findmexactf
        m0 = 0.1
        m1 = 9
        tol = 1e-9
        maxn = 1000
        ncnt = 2
        p0 = m0
        p1 = m1
        q0 = f(g, aa, p0)
        q1 = f(g, aa, p1)
    
        while ncnt <= maxn:
            p = p1-q1*(p1-p0)/(q1-q0)
            dif = p-p1
            ncnt = ncnt+1
            if abs(p-p1) < tol:
                break
            p0 = p1
            q0 = q1
            p1 = p
            q1 = f(g, aa, p1)
            m = p
        return m
    
    def findmexactsup(g,aa):
    
        f = findmexactf
        m0 = 3
        m1 = 4
        tol = 0.00001
        maxn = 1000
        ncnt = 2.5
        p0 = m0
        p1 = m1
        q0 = f(g, aa, p0)
        q1 = f(g, aa, p1)
    
        while ncnt <= maxn:
            p = p1-q1*(p1-p0)/(q1-q0)
            dif = p-p1
            ncnt = ncnt+1
            if abs(p-p1) < tol:
                break
            p0 = p1
            q0 = q1
            p1 = p
            q1 = f(g, aa, p1)
            m = p
        return m
    
    def findTTo(g,m):
        TTo = 1 + ((g - 1) / 2) * m * m
        return TTo
    
    def findppo(g,m):
        ppo = (1 + m * m * (g - 1) / 2) ** (g / (g - 1))
        return ppo
    
    def findrro(g,m):
        rro = (1 + m * m * (g - 1) / 2) ** (1 / (g - 1))
        return rro


    itercnt = np.append(itercnt,it)
    pp.append([rho[i]*Temp[i] for i in range(0,gridpoints)])
    ppp.append(rho[throat]*Temp[throat])
    mm.append([vel[i]/np.sqrt(Temp[i]) for i in range(0,gridpoints)])
    mmm.append(vel[throat]/np.sqrt(Temp[throat]))
    
    print('\n------------------------------------------------------')
    print('Final Results:', it,'iterations')
    print('------------------------------------------------------')
    print('\nFinal rho: ',[round(x,3)for x in rho])
    print('\nFinal temp: ',[round(x,3)for x in Temp])
    print('\nFinal vel: ',[round(x,3)for x in vel])
    print('\nFinal mach: ',[round(x,3)for x in [vel[i]/np.sqrt(Temp[i]) for i in range(0,gridpoints)]])
    print('\nFinal mdot: ',[round(x,3)for x in [rho[i]*vel[i]*Area[i] for i in range(0,gridpoints)]])
    
    
    Mach_exact = []
    rho_exact = []
    Temp_exact= []
    vel_exact = []
    pres_exact = []
    mdot_exact = []

    i = 0
    supersonic_check = 0
    while i < gridpoints:
        if supersonic_check == 0:
            Mach_exact.append(findmexactsub(1.4,Area[i]))
            rho_exact.append(1/findrro(1.4,Mach_exact[i]))
            Temp_exact.append(1/findTTo(1.4,Mach_exact[i]))
            vel_exact.append(Mach_exact[i] * np.sqrt(Temp_exact[i]))
            pres_exact.append(1/findppo(1.4,Mach_exact[i]))
            mdot_exact.append(rho_exact[i] * vel_exact[i] * Area[i])

        if supersonic_check == 1:
            Mach_exact.append(findmexactsup(1.4,Area[i]))
            rho_exact.append(1/findrro(1.4,Mach_exact[i]))
            Temp_exact.append(1/findTTo(1.4,Mach_exact[i]))
            vel_exact.append(Mach_exact[i] * np.sqrt(Temp_exact[i]))
            pres_exact.append(1/findppo(1.4,Mach_exact[i]))
            mdot_exact.append(rho_exact[i] * vel_exact[i] * Area[i])

        if Area[i] == 1:
            supersonic_check = 1
        i = i + 1
        
    a = Mach_exact[int((gridpoints-1)/2)]
    b = vel[int((gridpoints-1)/2)]/np.sqrt(Temp[int((gridpoints-1)/2)])
    diffab = (a-b)/((a+b)/2)*100
    c = vel_exact[int((gridpoints-1)/2)]
    d = vel[int((gridpoints-1)/2)]
    diffcd = (c-d)/((c+d)/2)*100
    e = findrro(1.4, a)
    f = findrro(1.4, b)
    diffef = (e-f)/((e+f)/2)*100
    g = findTTo(1.4, a)
    h = findTTo(1.4, b)
    diffgh = (g-h)/((g+h)/2)*100
    i = findppo(1.4, a)
    j = findppo(1.4, b)
    diffij = (i-j)/((i+j)/2)*100
    k = mdot_exact[int((gridpoints-1)/2)]
    l = rho[int((gridpoints-1)/2)]*vel[int((gridpoints-1)/2)]*Area[int((gridpoints-1)/2)]
    diffkl = (k-l)/((k+l)/2)*100
    
    print('\n-----------------------------------------------')
    print('\nExact Solution @ x = 1.500 (AoA* = 1.00)')
    print('-----------------------------------------------')
    print('\nMach: Exact =', a,' Solution = ',b,' %Diff =',abs(diffab),'%')
    print('\nVel:  Exact =', c,' Solution = ',d,' %Diff =',abs(diffcd),'%')
    print('\nR/Ro: Exact =', e,' Solution = ',f,' %Diff =',abs(diffef),'%')
    print('\nT/To: Exact =', g,' Solution = ',h,' %Diff =',abs(diffgh),'%')
    print('\nP/Po: Exact =', i,' Solution = ',j,' %Diff =',abs(diffij),'%')
    print('\nmdot: Exact =', k,' Solution = ',l,' %Diff =',abs(diffkl),'%')
    
    if gridpoints == 61 : 
      for i in range(0, gridpoints):
          pres[i] = rho[i]*Temp[i]
          mdot[i] = rho[i]*vel[i]*Area[i]
          mach[i] = vel[i]/np.sqrt(Temp[i])
        
      fig1, axs = plt.subplots(3, 2, figsize=(10, 8))
      axs[0,0].plot(xval,pp[0])
      axs[0,0].plot(xval,pp[1])
      axs[0,0].plot(xval,pp[2])
      axs[0,0].plot(xval,pp[3])
      axs[0,0].plot(xval,pp[4])
      axs[0, 0].set(xlabel=r'$x/L$',ylabel=r'$P/P_0$')
      axs[0, 0].set_title("Time Evolution")
    
      
      axs[0,1].plot(xval,mm[0])
      axs[0,1].plot(xval,mm[1])
      axs[0,1].plot(xval,mm[2])
      axs[0,1].plot(xval,mm[3])
      axs[0,1].plot(xval,mm[4])
      axs[0, 1].set(xlabel=r'$x/L$',ylabel='Mach Number')
      axs[0, 1].set_title("Time Evolution")

     
      axs[1,0].plot(itercnt, ppp)
      axs[1, 0].set(xlabel=r'$Time$',ylabel=r'$P/P_0$')
      axs[1, 0].set_title("Convergence at the Throat")
     
      axs[1,1].plot(itercnt, mmm)
      axs[1, 1].set(xlabel=r'$Time$',ylabel=r'$P/P_0$')
      axs[1, 1].set_title("Convergence at the Throat")

      axs[2,0].plot(xval,pp[4])
      axs[2, 0].set(xlabel=r'$x/L$',ylabel=r'$P/P_0$')
      axs[2, 0].set_title("Steady State Solution")
      
      axs[2,1].plot(xval,mm[4])
      axs[2, 1].set(xlabel=r'$x/L$',ylabel='Mach Number')
      axs[2, 1].set_title("Steady State Solution")
      fig1.suptitle('Supersonic Quasi-1D Flow (61 Grid Points| C = 0.5)')
      
      plt.show()
      
      fig1.tight_layout() 
      fig2, axs = plt.subplots(2, 2, figsize=(10, 8))
      axs[0, 0].plot(xval,pres_exact,color='black',marker='o',markersize=3,linestyle='None')
      axs[0, 0].plot(xval,pres,'-',color='red')
      axs[0, 0].set(xlabel=r'$x/L$',ylabel=r'$P/P_0$')
      axs[0, 0].legend(['Exact','CFD'])
      axs[0, 1].plot(xval,vel_exact,color='black',marker='o',markersize=3,linestyle='None')
      axs[0, 1].plot(xval,vel,color='blue')
      axs[0, 1].set(xlabel=r'$x/L$',ylabel=r'$u/a_0$')
      axs[0, 1].legend(['Exact','CFD'])

      axs[1, 0].plot(xval,Mach_exact,color='black',marker='o',markersize=3,linestyle='None')
      axs[1, 0].plot(xval,mach,color='green')
      axs[1, 0].set(xlabel=r'$x/L$',ylabel='Mach Number')
      axs[1, 0].legend(['Exact','CFD'])

      axs[1, 1].plot(xval,mdot_exact,color='black',marker='o',markersize=3,linestyle='None')
      axs[1, 1].plot(xval,mdot,color='cyan')
      axs[1, 1].set(xlabel=r'$x/L$',ylabel=r'$\rho uA$')
      axs[1, 1].set_ylim([0,1.5])
      axs[1, 1].legend(['Exact','CFD'])
       
      fig2.suptitle('Supersonic Quasi-1D Flow (61 Grid Points| C = 0.5)')
      fig2.tight_layout()
      plt.show()
       
       
       
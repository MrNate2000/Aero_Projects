import math as mm

def radians(x):
    y = (mm.pi/180)*x
    return y


def degrees(x):
    y = (180/mm.pi)*x
    return y 


def PM(g,v,m):
    zero = degrees(mm.sqrt((g+1)/(g-1))*mm.atan(mm.sqrt(((g-1)/(g+1))*(m**2-1)))-mm.atan(mm.sqrt(m**2-1)))-v
    return zero


def ToT(g,m):
    t0_t = 1+(g-1)/2*m**2
    return t0_t


def PoP(g,m):
    p0_p = (1+(g-1)/2*m**2)**(g/(g-1))
    return p0_p


def p2p1(p01_p1,p02_p2):
    p2_p1 = p01_p1*p02_p2**-1
    return p2_p1


def t2t1(t01_t1,t02_t2):
    t2_t1 = t01_t1*t02_t2**-1
    return t2_t1


def mu(m):
    angle = degrees(mm.asin(m**-1))
    return angle


def secant(fx,a,b,tol,maxn,g,v):
    i = 2
    q_0 = fx(g,v,a)
    q_1 = fx(g,v,b)
    while i <= maxn:
        p = b-q_1*(b-a)/(q_1-q_0)
        if abs(p-b)<tol:
            return p
        i += 1
        a = b
        q_0 = q_1
        b = p
        q_1 = fx(g,v,p)
        if maxn < i:
            print("Failed to converge")
            return
    return p

def findnu(g,m):
    v = degrees(mm.sqrt((g+1)/(g-1))*mm.atan(mm.sqrt(((g-1)/(g+1))*(m**2-1)))-mm.atan(mm.sqrt(m**2-1)))
    return v


def findmach(g,v):
    mach = secant(PM,1,6,.001,15,g,v)
    return mach

print('What would you like to do?\n')
print('1: Find Nu given Mach number\n')
print('2: Find Mach number given Nu\n')
print('3: Find M2 given M1 and theta\n')
print('4: Quit\n')

i = int(input("\nWhat would you like to do? "))

while i != 4:
    if i < 2:
        g = float(input('\nWhat is your gamma value: '))
        mach = float(input('\nWhat is your mach number: '))
        nu = findnu(g,mach)
        print('\nNu(M): ',nu)
    elif i < 3:
        g = float(input('\nWhat is your gamma value: '))
        nu = float(input('\nWhat is your Nu value: '))
        mach = findmach(1.4,nu)
        print('\nMach number: ',mach)
    elif i < 4:
        g = float(input('\nWhat is your gamma value: '))
        m1 = float(input('\nWhat is your M1 value: '))
        theta = float(input('\nWhat is your theta(degrees): '))
        mu1 = mu(m1)
        nu1 = findnu(g,m1)
        p01_p1 = PoP(g,m1)
        t01_t1 = ToT(g,m1)
        nu2 = theta + nu1
        m2 = findmach(g,nu2)
        mu2 = mu(m2)
        p02_p2 = PoP(g,m2)
        t02_t2 = ToT(g,m2)
        p2_p1 = p2p1(p01_p1,p02_p2)
        t2_t1 = t2t1(t01_t1,t02_t2)
        print('\n\n*** Results ***\n')
        print('M1:         ',m1)
        print('mu1:        ', mu1)
        print('nu1:        ',nu1)
        print('po1/p1:     ',p01_p1)
        print('to1/t1:     ',t01_t1)
        print('M2:         ',m2)
        print('mu2:        ',mu2)
        print('nu2:        ',nu2)
        print('po2/p2:     ',p02_p2)
        print('to2/t2:     ',t02_t2)
        print('p2/p1:      ',p2_p1)
        print('t2/t1:      ',t2_t1)
        print('\n************** Done **************\n')
    i = int(input('\n\nEnter an option to run, or enter 4 to quit: '))

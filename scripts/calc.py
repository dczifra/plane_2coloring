import numpy as np


def get_segment(eps):
    return (eps/(2)) - (np.cos(eps/2)*np.sin(eps/2))

def calc1():
    m = 0
    ind = -1
    #for alpha in np.linspace(0.001, np.pi/3, 100000):
    for alpha in np.linspace(0.001, 0.722734248+0.00001, 100000):
        diff = np.sin(alpha)/np.tan((np.pi-alpha)/2)
        beta = np.arcsin((np.cos(alpha)-diff)/2)
        gamma = (np.pi-2*alpha-2*beta)/2

        x = np.cos(beta)-np.sin(alpha)

        segment_hor = get_segment(2*beta)
        segment_ver = get_segment(alpha)
        trapez = np.sin(alpha)*(2*np.sin(beta)+1)/2
        triangle = (0.5*x)/2-get_segment(gamma)

        covered_area = segment_hor+2*segment_ver+trapez+2*triangle
        all_area_inf = 1/(2*np.sqrt(3))

        # If we place the pattern 60 deg below th others, than x should be:
        x_should = 0.5*np.tan(np.pi/6)
        e1 = (np.pi-gamma)/2
        e2 = np.arctan(0.5/x)
        #e,f = 2*(np.sin(alpha)+x+np.cos(np.pi-e1-e2)), (2*np.sin(np.pi/2-alpha-gamma/2))*2
        e,f = 2*(2*np.sin(alpha)+x), 2*(2*np.cos(alpha)-0.5)
        new_area2 = e*f/2
        if(covered_area/new_area2 > m):
            m = covered_area/new_area2
            #ind = 180*alpha/np.pi
            ind = {
                "radius":(alpha, beta, gamma, x, np.pi-e1-e2),
                "area":(covered_area, new_area2),
                "e,f":(e,f)}

        if(np.abs(x-x_should) < 0.00001):
            new_area = np.sqrt((2*np.cos(gamma/2))**2-1)*2

            print(covered_area/new_area2, new_area2, new_area)
            print("e,f: ", e,f, new_area2)
            print(f"Alpha: {alpha} beta: {beta} gamma: {gamma}  [x:{x}]==> Final 0: {covered_area*all_area_inf:.4f} [advanced: {covered_area/new_area:.4f}]")
    print(m, ind)

def calc2():
    pass

calc1()

import numpy as np


def get_segment(eps):
    return (eps/(2)) - (np.cos(eps/2)*np.sin(eps/2))

def calc1():
    for alpha in np.linspace(0.001, np.pi/2, 100000):
        diff = np.sin(alpha)/np.tan((np.pi-alpha)/2)
        beta = np.arcsin((np.cos(alpha)-diff)/2)
        gamma = (np.pi-2*alpha-2*beta)/2

        x = np.cos(beta)-np.sin(alpha)

        segment_hor = get_segment(2*beta)
        segment_ver = get_segment(alpha)
        trapez = np.sin(alpha)*(2*np.sin(beta)+1)/2
        triangle = (0.5*x)/2-get_segment(gamma)

        covered_area = segment_hor+2*segment_ver+trapez+2*triangle
        all_area_inf = 1/(2*np.sqrt(2))

        # If we place the pattern 60 deg below th others, than x should be:
        x_should = 0.5*np.tan(np.pi/6)
        if(np.abs(x-x_should) < 0.00001):
            base = 2
            new = 1+np.cos(gamma/2)
            ratio = base/new
            print(covered_area)
            covered_area = 0.6774

            print(2*180*beta/np.pi)

            #print(np.cos(gamma/2), base/new)
            #print(segment_hor, segment_ver, trapez, triangle, "Final", ratio*covered_area*all_area_inf)
            #print(f"Alpha: {alpha:.2f} beta: {beta:.2f} gamma: {gamma:.2f}  [x:{x}]==> Final 0: {covered_area*all_area_inf:.4f} [advanced: {ratio*covered_area*all_area_inf:.4f}]")
            
            print(f"Alpha: {alpha} beta: {beta} gamma: {gamma}  [x:{x}]==> Final 0: {covered_area*all_area_inf:.4f} [advanced: {ratio*covered_area*all_area_inf:.4f}]")

def calc2():
    pass

calc1()

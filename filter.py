import numpy as np
import pandas as pd
from io import StringIO
import math
import os

def filter(V,V0, U, c0, c1, c2):
    

    eta0 = 0.95
    Vm = 20.03
    delV = 0.2
    qv = 7.1
    Um = 24.0 * (24/3600)
    delU = 1.7 * (24/3600)

    etaV = max(int(V<V0), ((1-((V-V0)/qv)**2)/(1+np.exp((V-Vm)/delV))))
    etaU = 1/(1+np.exp((Um-U)/delU))
    etaTL = min(1, c0 + c1/(U-c2))

    eta = eta0*etaTL*etaU*etaV
    print('------------------------------------------probabilities (trailing loss, velocity, mag, total)--------------------------------------------')
    print(etaTL)
    print(etaU)
    print(etaV)
    print(eta)
    rng = np.random.uniform()

    return (rng < eta)

def is_on_right_side(x, y, xy0, xy1):
    x0, y0 = xy0
    x1, y1 = xy1
    a = float(y1 - y0)
    b = float(x0 - x1)
    c = - a*x0 - b*y0
    return a*x + b*y + c >= 0

def intersect(x0, y0, x1, y1, vertices):
    for i in range(4):
        x00 = x0
        y00 = y0
        x01 = (x1-x0)
        y01 = (y1-y0)
        x10 = vertices[2*i]
        y10 = vertices[2*i+1]
        x11 = vertices[(2*(i+1))%8] - x10
        y11 = vertices[(2*(i+1)+1)%8] - y10
        d = x11*y01 - x01 * y11
        s0 = ((1/d)*((x00 - x10) * y01 - (y00 - y10) * x01)) < 1
        s1 = ((1/d)*((x00 - x10) * y01 - (y00 - y10) * x01)) > 0
        t0 = ((1/d)*(-(-(x00 - x10)* y11 + (y00 - y10)* x11))) < 1
        t1 = ((1/d)*(-(-(x00 - x10)* y11 + (y00 - y10)* x11))) > 0

        if all([s0, s1, t0, t1]):
            return True

    return False
    

def test_point(x, y, vertices):
    num_vert = len(vertices)
    is_right = [is_on_right_side(x, y, vertices[i], vertices[(i + 1) % num_vert]) for i in range(num_vert)]
    all_left = not any(is_right)
    all_right = all(is_right)
    return all_left or all_right    

totaldays = 8978
surveys = [["CATALINA", 0.19, 0.36, 0.06], ["E12", 0.26, 0.42, 0.], ["F51", 0., 0.27, 0.2], ["G96", 0.56, 0.18, 0.], ["LINEAR", 0.19, 0.45, 0.], ["LONEOS", 0., 0.47, 0.], ["NEAT", 0., 0.25, 0.3], ["SPACEWATCH", 0., 0.17, 0.39], ["WISE", 0., 0.32, 0.78]]
cat = pd.read_csv('../neopop.cat', skiprows=5, sep='\s+', usecols=['!Name', 'H'])
ast = np.genfromtxt('ast.txt', dtype = 'str')
cat = cat[cat['!Name'].isin(ast)]
for a in cat['!Name']:
    cat.replace(a, a.replace("'",""))
N = 50
day = 0
d = 212
for y in range(1998, 2024):
    d = d%365
    while d <=365:
        s = 0
        for survey in surveys:
            j = '../surveys/{}/{}{:03d}.DAT'.format(survey[0], y, d)
            if os.stat(j).st_size == 0:
                s += 1
        if s == 9:
            d +=1
            day += 1
            continue
        print('{}{:03d}'.format(y, d))
        temp = ''
        with open("vEARTHMOON.dat", 'r') as ad:
            line_numbers = [day+4]
            for i, line in enumerate(ad):
                if i in line_numbers:
                    temp += line.strip()+'\n'
                elif i > line_numbers[-1]:
                    break
        with open("vall.dat", 'r') as ad:
            line_numbers = [*range(12+(N+1)*day, 11+(N+1)*(day+1))]
            for i, line in enumerate(ad):
                if i in line_numbers:
                    if (ast[(i-12)%(N)] in cat['!Name'].values):
                        temp += line.strip()+'\n'
                elif i > line_numbers[-1]:
                    break
        temp = np.genfromtxt(StringIO(temp))
        sdist = np.sqrt(np.diag(np.matmul(temp[1:,:3],np.transpose(temp[1:,:3]))))
        temp += -temp[0]
        dist = np.sqrt(np.diag(np.matmul(temp[1:,:3],np.transpose(temp[1:,:3]))))
        vel = (365.24*math.pi/180)*np.divide(np.diag(np.matmul(temp[1:,:3], np.transpose(temp[1:,3:6]))), dist)
        dec = np.arcsin(np.divide(temp[1:,2], dist))*(180/math.pi)
        RA = []
        for i in range(len(dist)):
            if np.divide(temp[i+1,1], dist[i]) > 0:
                RA.append(np.arccos(np.divide(np.divide(temp[i+1,0], dist[i]), np.cos(dec[i]*(math.pi/180))))*(180/math.pi))
            else:
                RA.append(360-np.arccos(np.divide(np.divide(temp[i+1,0], dist[i]), np.cos(dec[i]*(math.pi/180))))*(180/math.pi))
        V = cat['H'].to_numpy()+5*np.log10(np.multiply(dist, sdist))
        cat['V'] = V
        cat['U'] = vel
        cat['RA'] = RA
        cat['dec'] = dec


        cat = cat.reset_index(drop = True)
        cat = cat.T
        for survey in surveys:
            j = '../surveys/{}/{}{:03d}.DAT'.format(survey[0], y, d)
            if os.stat(j).st_size == 0:
                continue
            skydat = np.genfromtxt(j)
            if skydat.ndim == 2:
                for skycov in skydat:
                    V0 = skycov[-1]
                    skycov = [(skycov[2*i], skycov[2*i+1]) for i in range(4)]
                    for i in cat:
                        if test_point(cat[i]['RA'], cat[i]['dec'], skycov):
                            print(cat[i]['!Name'])
                            if filter(cat[i]['V'], V0, cat[i]['U'], survey[1], survey[2], survey[3]):
                                cat.pop(i)
                                break
            else:
                V0 = skydat[-1]
                skydat = [(skydat[2*i], skydat[2*i+1]) for i in range(4)]
                for i in cat:
                    if test_point(cat[i]['RA'], cat[i]['dec'], skycov):
                        print(cat[i]['!Name'])
                        if filter(cat[i]['V'], V0, cat[i]['U'], survey[1], survey[2], survey[3]):
                            cat.pop(i)
                            break
        
        cat = cat.T
        day += 1
        d += 1
        if day == totaldays:
            break
    if day == totaldays:
            break

print(day)
with open('found.txt', 'a') as results:
    np.savetxt(results, ast[np.isin(ast, cat['!Name'].values, invert = True)], fmt='%s')



#!/usr/bin/env python3
"""verifica.py -- corre o canal de Hartmann para varios Ha e compara f*Re com a
solucao exacta. Requer o solver mhdturbFoamQS e o OpenFOAM 6 no PATH."""
import math
import os
import re
import shutil
import subprocess

H = 0.05
NU = 4.1792e-4
RHO = 1000.0
SIGMA = 70.9
G = 0.0105
DH = 4 * H
HB = math.sqrt(SIGMA / (RHO * NU))
SOLVER = os.environ.get('SOLVER', 'mhdturbFoamQS')


def exacto(Ha):
    if Ha == 0:
        return 96.0
    t = math.tanh(Ha)
    return 32 * Ha * t / (1 - t / Ha)


def main():
    print('Ha      B0[T]      Ubar(sim)   f*Re(sim)   f*Re(exacto)   desvio')
    print('-' * 72)
    piores = []
    for Ha in (0.0, 0.5, 1.0, 5.0, 10.0):
        b = Ha / (H * HB)
        s = open('0/B0').read()
        s = re.sub(r'internalField\s+uniform \(0 [0-9.eE+-]+ 0\);',
                   'internalField   uniform (0 %.8g 0);' % b, s)
        s = re.sub(r'value uniform \(0 [0-9.eE+-]+ 0\);',
                   'value uniform (0 %.8g 0);' % b, s)
        open('0/B0', 'w').write(s)
        for d in os.listdir('.'):
            if re.match(r'^\d+(\.\d+)?$', d) and d != '0':
                shutil.rmtree(d)
        shutil.rmtree('postProcessing', ignore_errors=True)
        r = subprocess.run([SOLVER], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, universal_newlines=True)
        log = r.stdout
        if 'FOAM FATAL' in log:
            print('%5.2f  %8.4f   FALHOU -- ver log.ha%s' % (Ha, b, Ha))
            open('log.ha%s' % Ha, 'w').write(log)
            continue
        m = re.findall(r'volAverage\(\) of U\s*=\s*\(\s*([0-9.eE+-]+)', log)
        if not m:
            print('%5.2f  %8.4f   sem Ubar no log' % (Ha, b))
            continue
        Ub = float(m[-1])
        fRe = (2 * DH * G / Ub ** 2) * (Ub * DH / NU)
        ex = exacto(Ha)
        d = 100 * (fRe - ex) / ex
        piores.append(abs(d))
        print('%5.2f  %8.4f   %9.6f   %9.4f   %9.4f      %+7.3f %%'
              % (Ha, b, Ub, fRe, ex, d))
    if piores:
        print()
        print('desvio absoluto maximo: %.3f %%' % max(piores))
        print('criterio de aceitacao sugerido: < 1 %% em todos os Ha')


if __name__ == '__main__':
    main()

import numpy as np
import matplotlib.pyplot as mp

###------------------ VARIAVEIS ------------------###
column = 0

L = 0.02
d1 = 0.001
d2 = 0.001
n = 5
v0 = 1
vn = 0

e0 = 8.85418782 * 10**(-12)
er1 = 2
er2 = 4

# Epsilon correspondente ao segmento em que estamos calculando
def getEpsilon(n, s, d1, d2):
    if s * ((d1+d2) / n) < d1:
        return er1
    else:
        return er2

###------------------ MATRIZZES ------------------###
matrizSL = np.zeros(shape = (n, n))
matrizSLR = np.zeros(shape = (n-2, n-2))

matrizResposta = np.zeros(shape = (n-2, 1))
matrizTensoes = np.zeros(shape = (n, 1))

minemonicosTensoes = []

###------------------ MATRIZ DE SISTEMAS LINEAR 2x2 ------------------###
for i in range(n-1):
  if i == 0:
    matrizSL[i][column] = -(e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n)
    matrizSL[i][column+1] = (e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n)
    matrizSL[i+1][column] = (e0*getEpsilon(n, i+2, d1, d2))/((d1+d2)/n)
    matrizSL[i+1][column+1] = -(e0*getEpsilon(n, i, d1, d2))/((d1+d2)/n) - (e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n)
  else:
    matrizSL[i][column+1] = (e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n)
    matrizSL[i+1][column] = (e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n)
    matrizSL[i+1][column+1] = -(e0*getEpsilon(n, i+1, d1, d2))/((d1+d2)/n) - (e0*getEpsilon(n, i+2, d1, d2))/((d1+d2)/n)
  column+=1


###------------------ MATRIZ DE SISTEMA LINEAR 2x2 REDUZIDA ------------------###
for i in range(1, n-1):
  for j in range(1, n-1):
    matrizSLR[i-1][j-1] = matrizSL[i][j]


###------------------ MATRIZ DE TENSOES ------------------###
matrizResposta[0][0] = -((e0*getEpsilon(n, 1, d1, d2))/((d1+d2)/n) * v0)
matrizResposta[n-3][0] =  -((e0*getEpsilon(n, n-1, d1, d2))/((d1+d2)/n) * vn)



###------------------ SISTEMA LINEAR ------------------###
linearSolved = np.linalg.solve(matrizSLR, matrizResposta)

# adicionando valores de V0 e Vn
matrizTensoes[0][0] = v0
matrizTensoes[n-1][0] = vn

###------------------ MATRIZ DE TENSOES ------------------###
for i in range(n-2):
  matrizTensoes[i+1][0] = linearSolved[i][0]

###------------------ MINEMONICOS DE TENSOES ------------------###
for i in range(n):
  minemonicosTensoes.append('V_'+str(i))


###------------------ PLOT GRAFICO DAS TENSOES ------------------###
print('Plotagem do gráfico de tensões')
mp.plot(minemonicosTensoes, matrizTensoes)

###------------------ CALCULO DA CAPACITANCIA ------------------###
C = e0*er2*(L*L*(matrizTensoes[n-1][0] - matrizTensoes[n-2][0]))/(matrizTensoes[n-1][0] - matrizTensoes[0][0])

print('A capacitância é: ' + str(C) + 'F')
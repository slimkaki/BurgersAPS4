import numpy as np
import matplotlib.pyplot as plt
from math import ceil, sin, pi

class Mcdonalds(object):
    def __init__(self, dx, dy, dimensions, alpha, dt, K):
        self.dx = dx
        self.dy = dy
        self.dimensions = dimensions

        self.alpha = alpha
        self.dt = dt
        self.K = K
        self.Q_ponto = Q_ponto
        self.dx2 = dx**2
        self.dy2 = dy**2

    def create_concentration_matrix(self):
        nx = ceil(self.dimensions[0]/self.dx)
        ny = ceil(self.dimensions[1]/self.dy)
        matrix = np.zeros([nx, ny])
        return matrix

    def DiferencaCentralX(self, Cmenos1, Cmais1, atual):
        return (Cmais1 + Cmenos1 - 2*atual) / self.dx2

    def DiferencaCentralY(self, Cmenos1, Cmais1, atual):
        return (Cmais1 + Cmenos1 - 2*atual) / self.dy2

    def DiferencaAdiantadaX(self, Cxmais1, Cxmenos1):
        return (Cxmais1 - Cxmenos1)/ (2*self.dx)

    def DiferencaAdiantadaY(self, Cymais1, Cymenos1):
        return (Cymais1 - Cymenos1)/ (2*self.dy)

    def newC(self, d2Cx, d2Cy, dCx, dCy, atual, v):
        return self.dt * (self.K * (d2Cx + d2Cy) - self.alpha * dCx - v * dCy ) + atual

    def CalcV(self, x):
        return self.alpha * sin((pi/5)*x)

    def passo(self, matrix):
        copy = np.copy(matrix)
        for i in range(1, matrix.shape[0]-1):
            for j in range(1,matrix.shape[1]-1):
                v = self.CalcV(i)
                # v = 0
                u = self.alpha
                atual = matrix[i][j]
                d2Cx = self.DiferencaCentralX(matrix[i - 1][j], matrix[i + 1][j], atual)
                d2Cy = self.DiferencaCentralY(matrix[i][j + 1], matrix[i][j - 1], atual)
                dCx = self.DiferencaAdiantadaX(matrix[i + 1][j], matrix[i - 1][j])
                dCy = self.DiferencaAdiantadaY(matrix[i][j - 1], matrix[i][j + 1])
                nova_c = self.newC(d2Cx, d2Cy, dCx, dCy, atual, v)
                if nova_c < 0:
                    nova_c = 0
                copy[i, j] = nova_c
        #topo/base
        for i in range(copy.shape[0]):
            copy[i][0] = copy[i][1]
            copy[i][-1] = copy[i][-2]

        #direita/esquerda
        for j in range(copy.shape[1]):
            copy[0][j] = copy[1][j]
            copy[-1][j] = copy[-2][j]
            
        return copy

    def solve(self, matrix, tMax, a, b, Tdespejo, qc):
        list_t = np.arange(0, tMax + 0 * self.dt, self.dt)
        for t in list_t:
            matrix = self.passo(matrix)
            if t < Tdespejo:
                matrix[a][b] += qc * self.dt
            # self.plotColorMap(matrix, k)
        return matrix

    def parserText(self, matriz):
        stringF = ""
        for i in range(len(matriz)):
            stringF += "[   "
            for j in range(matriz.shape[1]):
                stringF += str(round(matriz[i][j], 5)) + "    "
            stringF += "]\n"
        with open("saida.txt", "w") as txt:
            txt.write(str(stringF))


    def plotColorMap(self, matrix):
        color_map = plt.imshow(matrix)
        color_map.set_cmap("viridis")
        plt.title("Dispersão do poluente no espaço 2D")
        plt.colorbar()
        plt.gca().invert_yaxis()
        plt.plot(a, b, "ro")
        
        plt.savefig("grafico_colormap1.png")
        plt.show()
    

# grupo 11
dx = 0.5
dy = 0.5
dimensions = [30, 20]
dt = 0.05 #s
alpha = 1
k  = 1 # m^2/s
Q_ponto = 100 # kg/ms
qc = Q_ponto / (dx * dy)
Tderramamento = 3 #s

tMax = 10*Tderramamento # s

n = 11
a = int(n / 1.4 / dx) 
b = int(60 / (n + 5) / dy) 

def checkConvergencia(K, dt, dx):
    if((dt/(dx**2)) < (1/(4*K))):
        print("Converge, bola que segue")
    else:
        print("Nao converge, azedo")

# # grupo 5
# dx = 0.5
# dy = 0.5
# dimensions = [30, 30]
# dt = 0.05 #s
# alpha = 1
# k  = 1 # m^2/s
# Q_ponto = 80 # kg/ms
# qc = Q_ponto / (dx * dy)
# Tderramamento = 2 #s
# u = 1
# v = 0
# tMax = 5 # s
# a = int(15 / dx)
# b = int(15 / dy)

checkConvergencia(k, dt, dx)


xtudo = Mcdonalds(dx, dy, dimensions, alpha, dt, k)

matrix = xtudo.create_concentration_matrix()

solved = xtudo.solve(matrix, tMax, a, b, Tderramamento, qc)
# print(round(solved[40][40], 3))

xtudo.parserText(solved)
xtudo.plotColorMap(solved)
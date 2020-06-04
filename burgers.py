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
        matrix = np.zeros([ny, nx])
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
        v = 0
        for i in range(1, matrix.shape[0]-1):
            for j in range(1,matrix.shape[1]-1):
                """
                i eh linha (y)
                j eh coluna (x)
                """
                v = self.CalcV(i)
                atual = matrix[i][j]
                d2Cx = self.DiferencaCentralX(matrix[i][j - 1], matrix[i][j + 1], atual)
                d2Cy = self.DiferencaCentralY(matrix[i + 1][j], matrix[i - 1][j], atual)
                dCx = self.DiferencaAdiantadaX(matrix[i][j + 1], matrix[i][j - 1])
                dCy = self.DiferencaAdiantadaY(matrix[i - 1][j], matrix[i + 1][j])
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
        list_t = np.arange(0, tMax, self.dt)
        plt.title(f"Dispersão do poluente no espaço 2D em {list_t[-1]+self.dt}s")
        for t in list_t:
            if t <= Tdespejo:
                matrix[b][a] += qc * self.dt
            matrix = self.passo(matrix)
            self.parserText(matrix)
        return matrix

    def parserText(self, matriz):
        stringF = ""
        for i in range(len(matriz)):
            stringF += "[   "
            for j in range(matriz.shape[1]):
                stringF += str(round(matriz[i][j], 5)) + "    "
            stringF += "]\n"
        with open("saida.txt", "w") as txt:
            txt.write("")
        with open("saida.txt", "a") as txt:
            txt.write(str(stringF))
            txt.write("\n")

    def plotColorMap(self, matrix):
        color_map = plt.imshow(matrix)
        color_map.set_cmap("viridis")
        plt.colorbar()
        plt.gca().invert_yaxis()
        plt.plot(a, b, "ro")
        
        plt.savefig("grafico_colormap1.png")
        plt.show()
    
def checkConvergencia(K, dt, dx):
    if((dt/(dx**2)) < (1/(4*K))):
        print("Converge, bola que segue")
    else:
        print("Nao converge, azedo")

grupo = 11

if grupo == 11:
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
    tMax = 10 * Tderramamento # s
    a = int((grupo / 1.4) / dx)
    b = int((60 / (grupo + 5)) / dy)

elif grupo == 5:
    # grupo 5
    dx = 0.5
    dy = 0.5
    dimensions = [30, 30]
    dt = 0.05 #s
    alpha = 1
    k  = 1 # m^2/s
    Q_ponto = 80 # kg/ms
    qc = Q_ponto / (dx * dy)
    Tderramamento = 2 #s
    tMax = 5 # s
    a = int(15 / dx)
    b = int(15 / dy)

checkConvergencia(k, dt, dx)

bigmac = Mcdonalds(dx, dy, dimensions, alpha, dt, k)

matrix = bigmac.create_concentration_matrix()

solved = bigmac.solve(matrix, tMax, a, b, Tderramamento, qc)

bigmac.parserText(solved)
bigmac.plotColorMap(solved)
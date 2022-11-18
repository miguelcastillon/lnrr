from cmath import pi
import numpy as np

# Input variables
X = np.array([[0,0.1,0.2],[1,0.1,0.2],[2,0.1,0.2],[3,0.1,0.2],[4,0.1,0.2],[5,0.1,0.2],[6,0.1,0.2],[7,0.1,0.2],[8,0.1,0.2],[9,0.1,0.2]])
Y = np.array([[0,1,0],[1,1,1],[2,1,2],[3,1,3],[4,1,3],[5,1,2],[6,1,1],[7,1,0]])
w = 0.1
sigma = 1
distance_threshold = 5

# Other variables
N = X.shape[0]
M = Y.shape[0]

def exp(x, y, sigma2):
  if np.linalg.norm(x - y) < distance_threshold:
    return np.exp(-np.linalg.norm(x - y)**2 / (2*sigma2))
  else:
    return 0

def computeP(w, sigma, X, Y):
  sigma2 = sigma**2
  c = (w*M)/((1-w)*N)*(2*pi*sigma2)**(3/2)

  K = np.zeros((M, N))
  for n in range(N):
    for m in range(M):
      K[m, n] = exp(X[n], Y[m], sigma2)

  P = np.zeros((M, N))
  for n in range(N):
    for m in range(M):
      if K[m, n] > 0:
        d=0
        for l in range(M):
          d += exp(X[n], Y[l], sigma2)
        P[m, n] = K[m, n] / (c+d)
  return P

P = computeP(w, sigma, X, Y)
print("P: ", P)
# P1 = P.sum(axis=1)
# PT1 = P.sum(axis=0)
# print("P1: ", P1)
# print("PT1: ", PT1)
# print("PX: ", np.matmul(P, X))



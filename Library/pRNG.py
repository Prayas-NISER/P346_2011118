# Defining a LCG random generator to generate numbers b/w [0,1]

def LCG(r0: float, n: int):
  a = 1103515245
  c = 12345
  m = 32768
  l = []

  for i in range (0, n):
    r0 = float(((a*r0 + c) % m)/m)
    l.append(r0)
  return l

# Defining a LCG random generator to generate numbers b/w [-k,k]

def Random_bw(r0: float, n: int, k: int):
  a = 1103515245
  c = 12345
  m = 32768
  l = []

  for i in range (0, n):
    r0 = float(((a*r0 + c) % m)/m)
    p0 = 2*k*r0 - k
    l.append(p0)
  return l

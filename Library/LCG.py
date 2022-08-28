# Defining a LCG random generator to generate numbers b/w [0,1]

def LCG(r0: float, n: int):
  a = 1103515245
  c = 12345
  m = 32768
  l = []

  for i in range (0, n):
    r0 = float(((a*r0 + c) % m)/m)
    p0 = 2*r0 - 1
    l.append(p0)
  return l

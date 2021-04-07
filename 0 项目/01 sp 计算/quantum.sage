# m的量子化[m]q
# var('q')
quantum_bracket = lambda m:(q**m-q**-m)/(q-q**-1)
def quantum_factorial(m):
    '''量子化的阶乘'''
    res = 1
    for i in range(1,m+1):
        res *= quantum_bracket(i)
    return res
def quantum_comb(m,k):
    '''量子化的组合数'''
    numer = quantum_factorial(m)
    denomin = quantum_factorial(m-k) * quantum_factorial(k)
    return numer/denomin
import pyperclip
import re,sympy
from sympy import Matrix

def save_vari(v,name):
    '''文件存储'''
    import pickle
    fout = open(name,'wb')
    pickle.dump(v,fout)
    fout.close()
    
def read_vari(name):
    '''读取文件'''
    import pickle
    fin = open(name,'rb')
    a = pickle.load(fin)
    fin.close()
    return a

def symbols(name,num=1):
    """生成若干符号变量"""
    if num==1:
        return (sympy.symbols(name),)
    return sympy.symbols((name+'%d ')*num%tuple(i+1 for i in range(num)))

def matlab_equations(mat):
    '''将sympy矩阵转化为若干线性方程'''
    assert isinstance(mat,sympy.Matrix)
    txt = str(mat).split('(')[-1][:-1] #提取矩阵部分
    n,m = mat.shape
    ys = (1,)+symbols('y',m-1)
    eqs = []
    for i in range(n):
        row = mat.row(i)
        eq = sum([a*b for a,b in zip(ys,row)])
        eqs.append(str(eq))
    return eqs,ys

def symmat_to_matlab(mat,n,copy=True):
    '''sympy的矩阵转为matlab的矩阵'''
    assert isinstance(mat,sympy.Matrix)
    syms = 'syms ' + ('x%d '*n)%tuple(i for i in range(1,n+1))
    txt =str(mat)
    txt = re.sub(r'\], \[',r']; [',txt)
    txt = syms + '\n' + txt.split('(')[-1][:-1]
    if copy:pyperclip.copy(txt)
    return txt

def symmat_to_command(mat,n,copy=True):
    '''mat转化为matlab中，解线性方程组的命令'''
    eqs,ys = matlab_equations(mat)
    ys = list(ys[1:])
    syms = 'syms ' + ('x%d '*n)%tuple(i for i in range(1,n+1))
    for i in ys:
        syms = syms + str(i) +' '
    command = 'sol = solve(['
    eqs_txt = ''
    for i,line in enumerate(eqs):
        eq = 'eq%d = '%(i+1) + line + ';\n'
        eqs_txt += eq
        command += 'eq%d'%(i+1) + ', ' 
    eqs_txt = eqs_txt[:-1]
    command = command[:-2] +  '],' + str(ys) + ')'
    txt = 'tic\n' + syms + '\n' + eqs_txt + '\n' + command + ';\ntoc\n'
    if copy: pyperclip.copy(txt)
    return txt

def read_solution(n,copy=True):
    '''获取matlab中得到的数据，n为变量个数'''
    s = '['
    for i in range(n):
        s += 'sol.y%d, '%(i+1) 
    s = s[:-2] + ']'
    if copy: pyperclip.copy(s)
    return s

def solve_matrix(mat):
    n = len(mat.row(0))
    b = sympy.Matrix([0 for i in range(len(mat.col(0)))]) #零矩阵
    para = sympy.symbols(('y%d '*n)%tuple(range(1,n+1)))
    return sympy.linsolve((mat,b),para)
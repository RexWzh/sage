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

def symbols(name,num):
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

def latex_table(mat,scale=1):
    '''表格的latex代码，默认表头深色，内容灰白相间，四周有线段'''
    import pyperclip
    mat = [[str(i) for i in line] for line in mat] #转字符串
    m,n = len(mat),len(mat[0]) #行列
    #版头版尾
    beg = """\\begin{table}[h]
\\centering
\\scalebox{%.3f}{
\\rowcolors{2}{gray!25}{white}
\\begin{tabular}{|c|%s|}
\\rowcolor{gray!50}
"""%(scale,'c'*(n-1))
    end = '\n\\end{tabular}}\n\\end{table}'
    #内容部分
    content = "\\hline\n"
    for i,line in enumerate(mat):
        for s in line[:-1]:
            content = content + s + " &"
        content = content + line[-1] + '\\\\'
        content = content + ('\\hline\n' if i in [0,m-1] else '\n') #首行加个换行
    content = content[:-1] #去掉换行
    out = beg + content + end
    print(out)
    pyperclip.copy(out)
    

def latex_matrix(mat,hrows=False,hcols=False):
    '''矩阵的latex代码，默认为align环境，不带行线，不带列线'''
    import pyperclip #用于剪贴板
    #确保为字符串
    mat = [[str(i) for i in line] for line in mat] #转字符串
    m,n = len(mat),len(mat[0]) #行列
    #设置列线行线
    if not hrows:
        hrows = {}
    if not hcols:
        hcols = {}
    s = ""
    for i in range(n+1):
        s = s + ('|c' if i in hcols else 'c')
    #版头版尾
    beg="\\begin{align*}\n\\left(\\begin{array}{%s}\n"%s[:-1]
    end = "\n\\end{array}\\right)\n\\end{align*}"
    #内容部分
    content = "\\hline\n" if 0 in hrows else ""
    for i,line in enumerate(mat):
        for s in line[:-1]:
            content = content + s + ' &'
        content = content + line[-1] + '\\\\'
        content = content + ('\\hline\n' if i+1 in hrows else '\n')
    content = (content[:-9] if i+1 in hrows else content[:-3])#去掉多余末尾
    out = beg+content+end
    print(out)
    pyperclip.copy(out)
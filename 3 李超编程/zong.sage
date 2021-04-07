import sympy,time,pyperclip
from sympy.matrices import Matrix #最后处理矩阵用
#导入工具包
load('pkg/Lie.sage') #李代数相关的工具
load('pkg/spo.sage') #李超代数spo的工具
load('pkg/MatSp.sage') #处理矩阵的工具
load('pkg/AlgebraBSC.sage') #结构常数相关工具
load('pkg/pbw.sage') #李代数pbw基的工具
load('pkg/pbws.sage') #李超代数pbw基的工具
load('pkg/rex.py') #解方程等相关工具
print('''导入相关工具包：
1. Lie.sage：李代数工具，以及C,D族李代数
2. spo.sage：李超代数工具，以及spo
3. MatSp.sage：矩阵空间，求矩阵关于矩阵基的线性表示
4. AlgebraBSC：求李/李超代数的结构常数，以及由结构常数生成代数
5. pbw.sage/pbws.sage 李/李超代数的pbw基
6. rex.py 解方程工具
''')

#V初始化数据
tt = time.time()
L = SpO(m,n) #李超代数spo(2m|2n)
mats = L.negative_matrixs + L.h_matrixs + L.positive_matrixs #矩阵信息
pos_m = len(L.positive_roots)
print('李超代数：spo(%d,%d)\n空间维数为：%d\n正根数目为：%d\n'%(2*m,2*n,len(mats),pos_m))
es = Lie.symbols('e',pos_m)
fs = Lie.symbols('f',pos_m)
hs = Lie.symbols('h',m+n) 
syms = fs + hs +es #用作结构常数的符号变量
#奇根对应的符号变量
roots = L.positive_roots 
ind1 = [roots.index(r) for r in L.odd_roots if r in roots]
roots = L.simple_roots
ind2 = [roots.index(r) for r in L.odd_roots if r in roots]
odds = [fs[i] for i in ind1] + [hs[i] for i in ind2] + [es[i] for i in ind1] 

#计算结构常数
print('计算结构常数\n')
N = AlgebraBSC.lie_sup_SC(mats,L.odd_mats)
'''
Lf = {i+1:Element(element,N,syms) for i,element in enumerate(fs)}
Le = {i+1:Element(element,N,syms) for i,element in enumerate(es)}
Lh = {i+1:Element(element,N,syms) for i,element in enumerate(hs)}
#'''
#pbw基
print('计算李超代数泛包络的pbw基')
vf = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(fs)}
ve = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(es)}
vh = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(hs)}
ide = SPBWElement({tuple():1},N,syms,odds) #幺元
vectors = L.vectors
print('初始化完毕，用时 %.3fs\n'%(time.time()-tt))

#权的正根分解
tt = time.time()
weight = L.basis[s-1] + L.basis[m+t-1];print('待讨论的权：',weight) #待讨论的权（正交基形式）
root = L.vector2root(weight) #转单根形式
dec = L.positive_dec(root) #正根分解（元组构成）
negative_roots = L.negative_roots #负根
dec_root = []
for line in dec:
    new = []
    for i,element in enumerate(line):
        if element:
            new = new + [negative_roots[i]]*element
    dec_root.append(tuple(new))
dec_vector = [tuple(L.root2vector(r) for r in line) for line in dec_root]
print('权空间维数：',len(dec)) #分解总数目
print('用时：%.3fs\n'%(time.time()-tt))

#相应的pbw基
basis = [tuple(fs[negative_roots.index(i)] for i in line) for line in dec_root]
pbw_basis = [ide.tuple2element(i) for i in basis]
print('权空间对应的pbw基：')
for i in pbw_basis:print(i)

#最高权向量
x = Lie.symbols('x',n)
z = list(Lie.symbols('z',m))
z[s-1] = x[t-1]+2*n+s-t-m-1
la = sum([i*j for i,j in zip(z,L.basis[:m])]) + sum([i*j for i,j in zip(x,L.basis[m:])])
#d_s+e_t: a_s = b_t+2n-s-t
print('\n最高权：',la)

#泛包络上计算结果res
tt = time.time()
res = [[ve[i]*element for element in pbw_basis] for i in range(1,m+n+1)]
print("计算泛包络上的乘积结果，乘积总数：",(m+n)*len(pbw_basis))
print('用时：%.3fs\n'%(time.time()-tt))

#化简为Verma模的元素
tt = time.time()
res_verma = [[i.to_verma(la,L,m+n) for i in line] for line in res]
print('转化为Verma模上结果，用时：%.3fs\n'%(time.time()-tt))

#转化为矩阵
tt = time.time()
data = []
for line in res_verma:
    mat,_ = ide.linear_rep(line)
    data.append(mat)
data = [mat for mat in data if mat] #去掉空项
mat = matrix.block(len(data),1,data)
mat_sym = Matrix(mat)
print('矩阵阶数为 %d*%d  '%mat.dimensions(),end='')
print('计算用时：%.3fs'%(time.time()-tt))
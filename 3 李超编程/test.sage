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

L = SpO(m,n) #李超代数spo(2m|2n)
mats = L.negative_matrixs + L.h_matrixs + L.positive_matrixs #矩阵信息
pos_m = len(L.positive_roots)
print('李超代数：spo(%d,%d)\n空间维数为：%d\n正根数目为：%d'%(2*m,2*n,len(mats),pos_m))
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

N = AlgebraBSC.lie_sup_SC(mats,L.odd_mats)

#pbw基
vf = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(fs)}
ve = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(es)}
vh = {i+1:SPBWElement(element,N,syms,odds) for i,element in enumerate(hs)}
ide = SPBWElement({tuple():1},N,syms,odds) #幺元
vectors = L.vectors

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

basis = [tuple(fs[negative_roots.index(i)] for i in line) for line in dec_root]
pbw_basis = [ide.tuple2element(i) for i in basis]
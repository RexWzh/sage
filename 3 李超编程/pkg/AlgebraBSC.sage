
def lie_super_b(m1,m2,odd_mats):
    '''李超计算'''
    if m1 in odd_mats and m2 in odd_mats:
        return m1*m2+m2*m1
    return m1*m2-m2*m1

class AlgebraBSC(object):
    '''由结构常数生成的代数，输入为tuple结构常数矩阵'''
    def __init__(self,N,syms=False):
        dim = len(N) #代数维数
        if not syms: #初始化符号基
            syms = Lie.symbols('v',dim) 
        self.syms = syms
        basis = [tuple(0 if i-j else 1 for i in range(dim)) for j in range(dim)]
        self.basis = [Element(i,N,syms) for i in basis]
        self.N = N
        self.dim = len(N)
        self.N_sym = self.tuple2sym_SC(N,syms)
    
    def __str__(self):
        return '结构常数生成的代数，基元为：\n'+str(self.syms)
    
    def __repr__(self):
        return self.__str__()
    
    @staticmethod
    def check_lie_identity(N):
        '''检查jacobi identity'''
        n = len(N)
        f = lambda i,j,k:N[i][j][k]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for t in range(n):
                        if sum([f(i,j,s)*f(k,s,t)+f(j,k,s)*f(i,s,t)+f(i,k,s)*f(s,j,t) for s in range(n)]):
                            return True
        return True
        
    
    @staticmethod
    def check_antisymmetric(N_sym):
        '''结构常数是否反对称'''
        Ns = matrix(N_sym)
        return bool(Ns.T==-Ns)
    
    @staticmethod
    def tuple2sym_SC(N,syms=False):
        '''矩阵形式的结构常数转为元组形式'''
        n = len(N)
        if not syms:
            syms = Lie.symbols('v',n)
        N_sym = [[sum([i*j for i,j in zip(syms,ele)])for ele in line]for line in N]
        return N_sym
    
    @staticmethod
    def sym2tuple_SC(N_sym,syms=False):
        n = len(N_sym)
        if not syms:
            syms = Lie.symbols('v',n)
        N = [[Lie.coefficients(i,syms) for i in line]for line in N_sym]
        return N
    
    @staticmethod
    def lie_SC(mats):
        '''获取李代数的结构常数（顺带验证运算封闭性）'''
        V = MatSp(mats)
        n = len(mats) #代数的空间维数
        N = [[0 for i in range(n)] for j in range(n)] #结构常数
        for i in range(n):
            for j in range(n):
                mat = mats[i]*mats[j]-mats[j]*mats[i]
                N[i][j] = tuple(V(mat))
        return N
    
    @staticmethod
    def lie_sup_SC(mats,odd_mats):
        '''获取李超代数的结构常数（顺带验证运算封闭性）'''
        V = MatSp(mats)
        n = len(mats) #代数的空间维数
        N = [[0 for i in range(n)] for j in range(n)] #结构常数
        for i in range(n):
            for j in range(n):
                mat = lie_super_b(mats[i],mats[j],odd_mats)
                N[i][j] = tuple(V(mat))
        return N
    
    """@staticmethod
    def sage_SC(T):
        '''sage内置的结构常数转化'''
        pass
    """
    
    def checke_lie_closure(basis):
        '''检查基元运算封闭性'''
        V = MatSp(basis)
        for m1 in basis:
            for m2 in basis:
                m = m1*m2-m2*m1
                if not m in V:
                    #print(m1,m2,m,sep='\n\n')
                    return False
        return True
    
    @staticmethod
    def check_sup_closure(basis,odd_mats):
        '''检查基元运算封闭性，basis需线性无关且奇部分都在odd_mats中'''
        V = MatSp(basis)
        for m1 in basis:
            for m2 in basis:
                m = lie_super_b(m1,m2,odd_mats)
                if not m in V:
                    #print(m1,m2,m,sep='\n\n')
                    return False
        return True


class Element():
    '''结构常数生成代数的元素
    内容：
    存储元素的符号及元组信息
    定义元素间的运算
    输入：
    输入元组或者符号变量
    '''
    def __init__(self,element,N,syms):
        assert isinstance(element,(tuple,Expression)),'输入格式有误！'
        dim = len(N)
        self.syms,self.N,self.dim = syms,N,dim
        if isinstance(element,tuple):
            self.element = element
            self.sym = self.element2sym(element)
        else:
            self.element = self.sym2element(element)
            self.sym = element
        
    def sym2element(self,sym):
        '''将符号变量转元组'''
        return tuple(sym.diff(i) for i in self.syms)

    def element2sym(self,element):
        '''将元组转符号变量'''
        return sum([i*j for i,j in zip(element,self.syms)])
    
    def __call__(self,element):
        return Element(element,self.N,self.syms)

    def __str__(self):
        return str(self.sym)

    def __repr__(self):
        return repr(self.sym)
    
    def __eq__(self):
        assert isinstance(obj,Element)
        return bool(self.sym==obj.sym)

    def __add__(self,obj):
        '''加法，obj为0或Element'''
        assert isinstance(obj,Element)
        return self(tuple(i+j for i,j in zip(self.element,obj.element)))

    def __sub__(self,obj):
        '''减法，obj为0或Element'''
        assert isinstance(obj,Element)
        return self(tuple(i-j for i,j in zip(self.element,obj.element)))

    def __neg__(self):
        return self(tuple(-i for i in self.element))

    def __mul__(self,obj):
        '''乘积，支持数乘与元素乘
        注意：数乘不检查数据类型！！'''
        if not isinstance(obj,Element): #若obj为数
            return self(tuple(obj*i for i in self.element))
        element1 = self.element
        element2 = obj.element
        dim = self.dim
        res = [0 for i in range(dim)]
        for i,e1 in enumerate(element1):
            if e1 ==0: continue
            for j,e2 in enumerate(element2):
                if e2==0: continue
                res = [a+e1*e2*b for a,b in zip(res,N[i][j])]
        return self(tuple(res))

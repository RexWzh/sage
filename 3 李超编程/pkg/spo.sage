class SpO(Lie):
    '''
    编写李超代数spo(2m|2n)的根系，矩阵基及一些操作工具
    
    ### 根系介绍 ###
    1.根系分even,odd，用符号向量给出两种表示方式
        一种是单根分解
        一种是正交基分解
    2.每个根对应矩阵，同样分even和odd两部分
    3.奇偶排序上：
        even按sp+,sp-,o+,o-排雷
        odd 按i,j从低到高排列（具体看定义过程）
        单根分解，正交基分解，以及奇偶部分的矩阵基，按位置对应
    4.整体排序上
        按单根分解的ht分级，内部按字典序排列。
        这部分函数继承自Lie
    '''
        
    def __init__(self,m=1,n=2):
        sp,o = Sp(m),Oth(n) #李代数对象sp和o
        d,e = sp.basis,o.basis
        basis = sp.basis + o.basis #空间的正交基
        self.sp,self.o,self.basis,self.m,self.n = sp,o,basis,m,n
        
        #偶根部分
        simple_vectors = sp.simple_vectors[:-1] + [ d[-1]-e[0] ] + o.simple_vectors #单根
        even_vectors = sp.vectors+o.vectors #偶根
        self.simple_vectors,self.even_vectors = simple_vectors,even_vectors
        #偶根对应的矩阵
        sp_mats = [matrix.block_diagonal(mat,matrix(2*n)) for mat in sp.matrixs]
        o_mats = [matrix.block_diagonal(matrix(2*m),mat) for mat in o.matrixs]
        even_mats = sp_mats + o_mats
        self.even_mats = even_mats
        
        #奇根部分
        odd_vectors = [] #正交基表示
        for i in range(m):
            for j in range(n):
                odd_vectors.append(d[i]+e[j])
                odd_vectors.append(-d[i]-e[j])
                odd_vectors.append(d[i]-e[j])
                odd_vectors.append(-d[i]+e[j])
        #奇根对应的矩阵
        odd_mats = []
        for i in range(m):
            for j in range(n):
                odd_mats.append(matrix(2*m+2*n,{(2*m+j,i+m):1,(i,j+n+2*m):1}))
                odd_mats.append(matrix(2*m+2*n,{(2*m+j+n,i):1,(i+m,j+2*m):-1}))
                odd_mats.append(matrix(2*m+2*n,{(2*m+j+n,i+m):1,(i,j+2*m):1}))
                odd_mats.append(matrix(2*m+2*n,{(2*m+j,i):1,(i+m,j+n+2*m):-1}))
        self.odd_vectors,self.odd_mats = odd_vectors,odd_mats
        vectors = even_vectors + odd_vectors
        mats = even_mats + odd_mats
        
        #单根与正交基的过渡矩阵
        simple_roots = self.symbols('a',m) + o.simple_roots #单根（其中o部分与原来相同）
        tran_mat = matrix.identity(n+m) #正交基到单根的过渡矩阵
        tran_mat[-2,-1] = 1
        for i in range(n+m-1):
            tran_mat[i+1,i] = -1
        self.tran_mat,self.tran_mat_inv = tran_mat,tran_mat.inverse()
        self.simple_roots = simple_roots
        
        #单根表示下的根系
        even_roots = [self.vector2root(vec) for vec in even_vectors]
        odd_roots = [self.vector2root(vec) for vec in odd_vectors]
        roots = [self.vector2root(vec) for vec in vectors]
        self.even_roots,self.odd_roots = even_roots,odd_roots
        
        #求正根，并排序根系
        positive_roots = [root for root in roots if sum(self.coefficients(root,simple_roots))>0]
        positive_vectors = [self.root2vector(root) for root in positive_roots]
        positive_matrixs = [mats[roots.index(r)] for r in positive_roots]
        ind = self.ordering(simple_roots,positive_roots)
        self.positive_vectors = [positive_vectors[i] for i in ind]
        self.positive_matrixs = [positive_matrixs[i] for i in ind]
        self.positive_roots = [positive_roots[i] for i in ind]
        self.negative_mats = [mats[roots.index(-r)] for r in self.positive_roots]
        
        #代数生成元矩阵
        self.e_matrixs = self.positive_matrixs[:n+m]
        self.f_matrixs = self.negative_matrixs[:n+m]
        #self.h_matrixs = [i*j-j*i for i,j in zip(self.e_matrixs,self.f_matrixs)]
        #'''
        self.h_matrixs = []
        for i,j in zip(self.e_matrixs,self.f_matrixs):
            mat = i*j + (j*i if i in odd_mats and j in odd_mats else -j*i)
            self.h_matrixs.append(mat)
        self.gen_matrixs = self.e_matrixs + self.f_matrixs + self.h_matrixs
        self.all_matrixs = self.matrixs + self.h_matrixs
        #'''
        #根对角阵（测试用）
        sym = basis[:m]+tuple(-i for i in basis[:m])+basis[m:]+tuple(-i for i in basis[m:])
        self.diag = matrix.diagonal(sym)
        #其他数据
        self.dim = len(self.all_matrixs)
        
    @property
    def negative_matrixs(self):
        '''李超正负根的矩阵不是转置关系，需重定义'''
        return self.negative_mats
    
    def lie_super_b(self,m1,m2):
        '''李超括号（暂时仅支持基元运算）'''
        assert m1 in self.all_matrixs and m2 in self.all_matrixs, "仅支持矩阵基上计算！"
        if m1 in self.odd_mats and m2 in self.odd_mats:
            return m1*m2 + m2*m1
        return m1*m2-m2*m1
    
    def inner_product(self,r1,r2,is_root=False):
        '''计算内积'''
        v1,v2 = (self.root2vector(r1),self.root2vector(r2)) if is_root else (r1,r2)
        return self.sp.inner_product(v1,v2)-self.o.inner_product(v1,v2)
    
    def root_on_diag(self,alpha,h,is_root=False):
        '''根对对角阵的作用'''
        assert matrix.diagonal(h.diagonal())==h, 'h必须为对角阵'
        diag,m,n = h.diagonal(),self.m,self.n
        assert diag[:m]==[-i for i in diag[m:2*m]],'h必须在Cartan子代数上！'
        assert diag[2*m:2*m+n]==[-i for i in diag[2*m+n:]],'h必须在Cartan子代数上！'
        if is_root:
            alpha = self.root2vector(alpha)
        coef = self.coefficients(alpha,self.basis)
        return sum([i*j for i,j in zip(diag[:m]+diag[2*m:2*m+n],coef)])
    
    def super_trace(self,mat):
        '''计算矩阵的超迹'''
        m,n = self.m,self.n
        return mat[:2*m,:2*m].trace()-mat[-2*n:,-2*n:].trace()
    
    def B(self,m1,m2):
        '''双线性型'''
        return self.super_trace(m1*m2)
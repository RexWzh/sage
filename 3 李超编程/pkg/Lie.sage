class Lie():
    '''
    李代数的基础函数工具
    '''
    
    @property
    def negative_roots(self):
        return [-r for r in self.positive_roots]

    @property
    def roots(self):
        return self.positive_roots+self.negative_roots

    @property
    def negative_vectors(self):
        return [-r for r in self.positive_vectors]

    @property
    def vectors(self):
        return self.positive_vectors+self.negative_vectors
    
    @property
    def negative_matrixs(self):
        return [m.transpose() for m in self.positive_matrixs]
    
    @property
    def matrixs(self):
        return self.positive_matrixs+self.negative_matrixs
    
    @staticmethod
    def coefficients(a,basis):
        '''符号向量的系数'''
        return tuple(a.diff(i) for i in basis)
    
    @classmethod
    def transform(cls,vec,tran_mat,basis1,basis2):
        '''求向量在另一组基下的坐标
        vec为关于basis1的坐标（行向量）
        basis1 = basis2*tran_mat
        '''
        coef = cls.coefficients(vec,basis1)
        coef = matrix(coef)*tran_mat.transpose()
        return sum([i*j for i,j in zip(coef.row(0),basis2)])
    
    @staticmethod
    def symbols(name,num):
        """生成若干符号变量"""
        if num==1:
            return (var(name),)
        return var((name+'%d ')*num%tuple(i+1 for i in range(num)))
    
    def root2vector(self,root):
        '''符号根转向量（正交基表示）'''
        tran_mat,basis1,basis2 = self.tran_mat,self.simple_roots,self.basis
        return self.transform(root,tran_mat,basis1,basis2)
    
    def vector2root(self,vec):
        '''向量转符号根（单根表示）'''
        tran_mat,basis1,basis2 = self.tran_mat_inv,self.basis,self.simple_roots
        return self.transform(vec,tran_mat,basis1,basis2)
    
    def root2matrix(self,root):
        '''根转矩阵（仅对根系有效）'''
        roots = self.roots
        assert root in roots, "根必须在根系上！"
        k = roots.index(root) #定位
        return self.matrixs[k]
    
    def vector2matrix(self,vec):
        '''向量转矩阵（仅对根系有效）'''
        vectors = self.vectors
        assert vec in vectors, "向量必须在根系上"
        k = vectors.index(vec) #定位
        return self.matrixs[k]
    
    @classmethod
    def ordering(cls,basis,roots):
        '''正根排序：先按ht排序，再按坐标排'''
        coef = [cls.coefficients(r,basis) for r in roots] #转系数形式
        new = coef.copy()
        new.sort(reverse=True) #先按坐标排列
        new.sort(key=sum) #再按ht排序
        return tuple(coef.index(i) for i in new)
    
    @staticmethod
    def lie_bracket(a,b):
        '''李括号计算'''
        return a*b-b*a
    
    def inner_product(self,r1,r2,is_root=False):
        '''计算内积'''
        v1,v2 = (self.root2vector(r1),self.root2vector(r2)) if is_root else (r1,r2)
        return sum([v1.diff(i)*v2.diff(i) for i in self.basis])
    
    def cartan_integer(self,beta,alpha,is_root=False):
        '''计算Cartan整数'''
        return 2*self.inner_product(alpha,beta,is_root)/self.inner_product(alpha,alpha,is_root)
    
    def positive_dec(self,root):
        '''正根分解，即P函数'''
        coef = self.coefficients(root,self.simple_roots)
        assert all([i>=0 for i in coef]),'输入必须为正的根'
        positive = self.positive_roots #正根系
        n = len(positive) #向量总长
        coef = coef + tuple(0 for i in range(n-len(coef))) #补0凑等长
        res = [{coef}] #最长解
        while True:
            new = set()
            for root in res[-1]: #上一层解
                for i in range(n): #抓一根
                    if not root[i]: continue
                    for j in range(i,n): #再抓一根
                        if not root[j]: continue
                        if i==j and root[i]==1: continue
                        r = positive[i] + positive[j] #求和得到的根
                        if not r in positive: continue
                        k = positive.index(r)
                        new_coef = list(root)
                        new_coef[i] -= 1
                        new_coef[j] -= 1
                        new_coef[k] += 1
                        new_coef = tuple(new_coef)
                        new.add(new_coef)
            if not len(new): break #没有新的解
            res.append(new) #加入长度-1的解
        return [i for line in res for i in line]
    
    def root_on_diag(self,alpha,h,is_root=False):
        '''根对对角阵的作用'''
        assert matrix.diagonal(h.diagonal())==h, 'h必须为对角阵'
        diag,n = h.diagonal(),self.n
        assert diag[:n]==[-i for i in diag[n:]],'h必须在Cartan子代数上！'
        if is_root:
            alpha = self.root2vector(alpha)
        coef = self.coefficients(alpha,self.basis)
        return sum([i*j for i,j in zip(diag[:n],coef)])
    
    #测试内容
    @staticmethod
    def check_weight(h,mat,vec):
        '''检查mat的权是否为vec'''
        return True if vec*mat == h*mat-mat*h else False
    
    @classmethod
    def check_weights(cls,h,mats,vectors):
        for i,j in zip(mats,vectors):
            if not cls.check_weight(h,i,j):
                return False
        return True
    
    @staticmethod
    def find_weight(h,mat):
        '''求mat的权'''
        if mat.rank() == 0:return 'zero matrix'
        new = h*mat-mat*h
        for i,row in enumerate(mat):
            for j,ele in enumerate(row):
                if ele != 0:break
            else:continue #没有break继续
            break
        root = new[i,j]/ele
        return root if new == root*mat else 'not weight matrix!'
    

class Sp(Lie):
    '''Sp的根系和矩阵基'''
    def __init__(self,n,root_name='a',basis_name='delta'):
        #符号根
        self.n = n
        simple_roots = self.symbols(root_name,n) #sp的符号单根
        basis = self.symbols(basis_name,n) #sp的符号规范正交基
        #两组基之间的过渡矩阵（单根=正交基*mat）
        tran_mat = matrix.identity(n) 
        tran_mat[-1,-1] = 2
        for i in range(n-1):
            tran_mat[i+1,i] = -1
        #添加为类的属性
        self.tran_mat,self.tran_mat_inv = tran_mat,tran_mat.inverse()
        self.simple_roots,self.basis = simple_roots,basis

        #规范正交基下的正根
        positive_vectors = [2*basis[i] for i in range(n)] 
        for i in range(n-1):
            for j in range(i+1,n):
                positive_vectors.append(basis[i]-basis[j])
                positive_vectors.append(basis[i]+basis[j])
        #符号单根下的正根
        positive_roots = [self.vector2root(i) for i in positive_vectors] #等价写法
        #正根对应的矩阵
        positive_matrixs = [matrix(2*n,{(i,n+i):1}) for i in range(n)]
        for i in range(n-1):
            for j in range(i+1,n):
                positive_matrixs.append(matrix(2*n,{(i,j):1,(n+j,n+i):-1}))
                positive_matrixs.append(matrix(2*n,{(i,n+j):1,(j,n+i):1}))
        #数据排序，并添加到类的属性
        ind = self.ordering(simple_roots,positive_roots) #按单根排序
        self.positive_vectors = [positive_vectors[i] for i in ind]
        self.positive_matrixs = [positive_matrixs[i] for i in ind]
        self.positive_roots = [positive_roots[i] for i in ind]
        self.simple_vectors = self.positive_vectors[:n]
        #代数生成元矩阵
        self.e_matrixs = self.positive_matrixs[:n]
        self.f_matrixs = self.negative_matrixs[:n]
        self.h_matrixs = [self.lie_bracket(i,j) for i,j in zip(self.e_matrixs,self.f_matrixs)]
        self.gen_matrixs = self.e_matrixs + self.f_matrixs + self.h_matrixs
        
        #根对角阵（测试用）
        sym = self.basis + tuple(-i for i in self.basis)
        self.diag = matrix.diagonal(sym)
        self.all_matrixs = self.matrixs+self.h_matrixs
        self.dim = len(self.all_matrixs)

class Oth(Lie):
    '''
    正交李代数，D族李的根系和矩阵基(n>=2)
    命名上，sp单根和规范正交基为alpha，delta，用a,d表示
    o的单根和规范正交基为beta,epsilon，用b,e表示
    '''
    def __init__(self,n,root_name='b',basis_name='epsilon'):
        #符号根
        self.n = n
        simple_roots = self.symbols(root_name,n) #sp的符号单根beta
        basis = self.symbols(basis_name,n) #sp的符号规范正交基epsilon
        #两组基之间的过渡矩阵（单根=正交基*mat）
        tran_mat = matrix.identity(n) 
        tran_mat[-2,-1] = 1
        for i in range(n-1):
            tran_mat[i+1,i] = -1
        #添加为类的属性
        self.tran_mat,self.tran_mat_inv = tran_mat,tran_mat.inverse()
        self.simple_roots,self.basis = simple_roots,basis

        #规范正交基下的正根
        positive_vectors = [] 
        for i in range(n-1):
            for j in range(i+1,n):
                positive_vectors.append(basis[i]-basis[j])
                positive_vectors.append(basis[i]+basis[j])
        
        #符号单根下的正根
        positive_roots = [self.vector2root(i) for i in positive_vectors] #等价写法
        #正根对应的矩阵
        positive_matrixs = []
        for i in range(n-1):
            for j in range(i+1,n):
                positive_matrixs.append(matrix(2*n,{(i,j):1,(n+j,n+i):-1}))
                positive_matrixs.append(matrix(2*n,{(i,n+j):1,(j,n+i):-1}))
        #数据排序，并添加到类的属性
        ind = self.ordering(simple_roots,positive_roots) 
        self.positive_vectors = [positive_vectors[i] for i in ind]
        self.positive_matrixs = [positive_matrixs[i] for i in ind]
        self.positive_roots = [positive_roots[i] for i in ind]
        self.simple_vectors = self.positive_vectors[:n]
        #代数生成元矩阵
        self.e_matrixs = self.positive_matrixs[:n]
        self.f_matrixs = self.negative_matrixs[:n]
        self.h_matrixs = [self.lie_bracket(i,j) for i,j in zip(self.e_matrixs,self.f_matrixs)]
        self.gen_matrixs = self.e_matrixs + self.f_matrixs + self.h_matrixs
        
        #根对角阵（测试用）
        sym = self.basis + tuple(-i for i in self.basis)
        self.diag = matrix.diagonal(sym)
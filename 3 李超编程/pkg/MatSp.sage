class MatSp():
    '''定义矩阵空间'''
    def __init__(self,basis,field=QQ):
        assert len(basis),'基元个数不能为空'
        self.basis = [matrix(mat) for mat in basis]
        
        m,n = basis[0].dimensions()
        basis = [self.mat2vector(mat) for mat in basis] #转向量
        self.nrows,self.ncols = m,n
        
        V = VectorSpace(field,m*n) #生成向量空间
        U = V.subspace_with_basis(basis) #子空间
        self.V,self.U,self.dimension = V,U,U.dimension()
    
    def __str__(self):
        txt = '%dx%d的矩阵空间，维数为%d\n'%(self.nrows,self.ncols,self.dimension)
        return txt
    
    def __repr__(self):
        return self.__str__()
    
    def __contains__(self,mat):
        '''判断包含'''
        vec = self.V(self.mat2vector(mat)) #转向量空间元素
        return vec in self.U
    
    def __call__(self,mat):
        '''快速调用矩阵'''
        return self.coordinate_vector(mat)
        
    def coordinate_vector(self,mat):
        '''求mat的坐标'''
        U,V = self.U,self.V
        vec = V(self.mat2vector(mat))
        assert vec in U,"矩阵不在空间上！"
        return U.coordinate_vector(vec)
    
    @staticmethod
    def mat2vector(mat):
        return [i for line in mat for i in line]
    
    @staticmethod
    def vector2mat(vec,nrows=False):
        '''向量转矩阵(默认方阵)'''
        if nrows:
            m,n = nrows,len(vec)/nrows
        else:
            n = sqrt(len(vec))
            m = n
        return matrix(QQ,[[vec[j*n+i] for i in range(n)] for j in range(m)])
class BilinearForm():
    '''双线性型
    输入
        sym 基元变量名
        gram 格拉姆矩阵
        ID 用于运算时检验数据类型
        display 数据显示方式
    属性
        .syms 基元变量名
        .gram 格拉姆矩阵
        .dim 空间维数
        .ID 运算ID
        .disaplay 显示方式（默认为符号变量）
        .zero 空间的零向量
    方法
        .basis 空间的一组基
        .vec2element 向量转当前数据类型
        .product_matrix 求给定向量集的乘积矩阵（向量应为当前数据类型，不检查）
    内置类
        Element 用于向量的基本操作
    注意事项
        内积通过幂次符号 ^ 或 ** 运算
        双线性型上自带偏序：v1>=v2 <=> v1-v2 所有系数非负
    '''
    def __init__(self, dim, sym=None, gram=None, ID=None, display='sym'):
        '''sym 基元名称 e.g. alpha
        gram 默认单位矩阵
        ID 用于运算时检验数据类型'''
        assert dim>0,'维数必须为正'
        sym = 'v' if sym is None else sym
        sym = var((sym+"%d ")*dim%tuple(range(1,dim+1))) # 标准基的符号
        self.syms = (sym,) if dim ==1 else sym # 一阶需确保为元组
        self.gram = matrix.identity(dim) if gram is None else matrix(gram)
        self._basis = None    
        self.dim = dim
        self.ID = ID
        self.display = display
        self.zero = self.Element(self,vector(0 for i in range(dim)))
        
    def basis(self):
        '''双线性型的标准基'''
        if self._basis is None:
            self._basis = []
            for line in matrix.identity(self.dim):
                self._basis.append(self.Element(self,line))
        return self._basis
    
    @staticmethod
    def product_matrix(vectors):
        '''计算vectors双线性型下的矩阵'''
        return matrix([[v1^v2 for v1 in vectors] for v2 in vectors])
    
    def vec2element(self,vec):
        '''向量转元素'''
        assert len(vec) == self.dim,'维数不一致！'
        return self.Element(self,vec)
    
    def __str__(self):
        return 'bilinear form on a vector space of dimension %d'%dim
    
    class Element():
        '''用于数据运算
        
        属性及方法
            .parent 访问父类数据
            .vector 向量
            .sym_vector 向量的符号变量形式（用于打印）
            .vec2sym() 向量转变量
            
        支持运算
            __pow__ 内积运算，通过幂次调用
            __add__ 右加，检查数据类型
            __radd__ 左加，用于累和，调用__add__
            __neg__ 取负
            __sub__ 作差，调用__add__
            __mul__ 数乘，或矩阵右乘
            __rmul__ 数乘，或矩阵左乘
            __truediv__ 右除数值
            __eq__ 判断相等
            __le__,__ge__ 偏序下的小于等于，大于等于
            __lt__,__gt__ 偏序下的严格小于，严格大于
        
        待添加
            向量换基（参考之前写的MatSpace）
        '''
        def __init__(self,parent,vec):
            self.parent = parent
            self.vector = vector(vec)
            self.sym_vector = self.vec2sym(vec)
        
        def __add__(self,vec):
            '''向量加法'''
            assert vec.parent.ID == self.parent.ID,'数据类型不匹配'
            return self.parent.vec2element(self.vector + vec.vector)
        
        def __radd__(self,other):
            '''右加函数，用于被sum调用'''
            if other == 0: # sum 函数从0开始
                return self
            return self.__add__(other)
        
        def __neg__(self):
            '''向量取负'''
            return self.parent.vec2element(-self.vector)
        
        def __sub__(self,vec):
            '''向量减法'''
            return self.__add__(-vec)
        
        def __mul__(self,obj):
            '''数乘或矩阵右乘'''
            return self.parent.vec2element(self.vector*obj)
        
        def __rmul__(self,obj):
            '''数乘或矩阵左乘'''
            return self.parent.vec2element(obj*self.vector)
        
        def __truediv__(self,obj):
            '''数量除法'''
            return self.parent.vec2element(self.vector/obj)
        
        def __pow__(self,vec):
            '''gram矩阵下的向量内积'''
            assert vec.parent.ID == self.parent.ID,'数据类型不匹配'
            gram = self.parent.gram
            return self.vector*gram*vec.vector
        
        def __eq__(self,obj):
            '''判别相等'''
            if obj == 0: obj = self.parent.zero
            assert self.parent.ID == obj.parent.ID
            return set(self.vector-obj.vector)=={0}
        
        def __ge__(self,obj):
            '''判别大于等于'''
            if obj == 0: obj = self.parent.zero
            assert self.parent.ID == obj.parent.ID
            vec = (self.vector-obj.vector)
            return not any([i<0 for i in vec])
        
        def __le__(self,obj):
            '''判别小于等于'''
            if obj == 0: obj = self.parent.zero
            assert self.parent.ID == obj.parent.ID
            vec = (self.vector-obj.vector)
            return not any([i>0 for i in vec])
        
        def __gt__(self,obj):
            '''判别严格大于'''
            return self.__ge__(obj) and not self.__eq__(obj)
        
        def __lt__(self,obj):
            '''判别严格小于'''
            return self.__le__(obj) and not self.__eq__(obj)
        
        def __str__(self):
            return str(self.sym_vector) if self.parent.display=='sym' else str(self.vector)
        def __repr__(self):
            return self.__str__()
        
        def vec2sym(self,vec):
            '''向量转符号变量形式'''
            syms = self.parent.syms
            return sum([i*j for i,j in zip(vec,syms)])
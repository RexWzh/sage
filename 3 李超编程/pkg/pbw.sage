class PBWElement():
    '''pbw基元素运算
    输入：
    元组或者符号变量
    element为字典类型'''
    def __init__(self,element,N,syms):
        '''初始元素在U1上'''
        assert isinstance(element,(dict,Expression)),'输入基元格式不对'
        L_dim = len(N) #李代数的维数
        self.syms = syms #李代数的符号基，用作pbw序
        if isinstance(element,Expression):
            element = Lie.coefficients(element,syms)
            element = {(syms[i],):element[i] for i in range(L_dim) if element[i]} #字典形式
        else:
            element = {key:val for key,val in zip(element.keys(),element.values()) if val}
        self.element = element #统一化字典形式，并消去零元
        self.L_dim,self.N = L_dim,N
        self.keys = element.keys()
    
    def __call__(self,obj):
        '''元素转pbw类'''
        return PBWElement(obj,self.N,self.syms)

    def __add__(self,obj):
        '''pbw元素加法'''
        assert isinstance(obj,PBWElement),'加法必须为pbw类'
        a,b = self.element,obj.element
        return self(self.add_dict(a,b))

    def __sub__(self,obj):
        '''pbw元素减法'''
        assert isinstance(obj,PBWElement),'减法必须为pbw类'
        a,b = self.element,obj.element
        b = {key:-val for key,val in zip(b.keys(),b.values())}
        return self(self.add_dict(a,b))
    
    def __neg__(self):
        '''pbw取负'''
        c = self.element
        return self({key:-val for key,val in zip(c.keys(),c.values())})
    
    def __bool__(self):
        '''判断是否为0'''
        return bool(self.element)
    
    def __eq__(self,obj):
        '''判断相等'''
        assert isinstance(obj,PBWElement)
        return bool(self.element==obj.element)
    
    def __mul__(self,obj):
        '''pbw元素乘法，，支持数乘与元素乘(乘法指右乘obj)
            注意：数乘不检查数据类型！！'''
        if not isinstance(obj,PBWElement): #若obj为数
            c = {key:obj*val for key,val in zip(self.element.keys(),self.element.values())}
            return self(c)
        a,b = self.element,obj.element
        res = {}
        for k1 in a.keys():
            for k2 in b.keys():
                res = self.add_dict(res,self.tuple2dict(k1+k2,a[k1]*b[k2]))
        return self(res)
    
    def __pow__(self,k):
        '''元素的次幂'''
        assert isinstance(k,(int,Integer)),'仅支持整数次幂！'
        a = self
        for i in range(k-1):
            a = a*self
        return a
    
    def __contains__(self,key):
        return bool(key in self.element)
    
    def __getitem__(self,key):
        return self.element[key]
    
    def __str__(self):
        element = self.element
        element = {self.element_simplify(key):val for key,val in zip(element.keys(),element.values())}
        return str(element)

    def __repr__(self):
        return self.__str__()
    
    def linear_rep(self,elements):
        '''将elements用pbw基+矩阵线性表示，elements为元素或符号变量元组'''
        if isinstance(elements[0],tuple): #若输入元组，进行转格式
            elements = [self.tuple2element(element) for element in elements ]
        keys = set()
        for ele in elements:
            keys = keys.union(ele.keys)
        mat = matrix([[ele.element[key] if key in ele else 0 for ele in elements] for key in keys])
        basis = [self({key:1}) for key in keys]
        return mat,basis
    
    def tuple2dict(self,element,c=1):
        '''将符号变量元组转为字典，系数默认为1'''
        if not c:return {} #零元情形
        syms = self.syms
        n = len(element) #元素长度
        N = self.N #结构常数
        coefs = [syms.index(i) for i in element] #获取对应系数
        for i in range(n-1):
            x,y = coefs[i],coefs[i+1]
            if x>y: #发现反序  
                e = element[:i] + (element[i+1],element[i]) + element[i+2:] #逆序处理
                seq = [(c*coef,sym) for coef,sym in zip(N[x][y],syms) if coef] #李括号项
                seq = [(coef,element[:i]+(sym,)+element[i+2:]) for coef,sym in seq]
                break
        else: #一切正序
            return {element:c}
        res = self.tuple2dict(e,c)
        for coef,sym in seq:
            new = self.tuple2dict(sym,coef)
            res = self.add_dict(res,new)
        return res
    
    @property
    def identity(self):
        '''返回单位元'''
        return self({tuple():1})
    
    def tuple2element(self,t):
        '''将系列符号元组按顺序乘积转化为元素'''
        element = self({tuple():1})
        for i in t:
            element = element*self(i)
        return element
    
    @staticmethod
    def element_simplify(element):
        '''将元素key值的重复项合并'''
        if not len(element):return element
        res = [element[0]]
        k = 0 #位置校准
        for i in range(len(element)-1):
            if element[i] == element[i+1]:
                res[i-k] *= element[i]
                k += 1
            else:
                res.append(element[i+1])
        return tuple(res)
    
    ##定义一个Un基的函数
    @staticmethod
    def add_dict(a,b):
        '''字典相加合并'''
        c = {key:a[key] for key in a.keys() if key not in b}
        for key in b.keys():
            c[key] = (a[key]+b[key]) if key in a else b[key]
        return c
    
    def to_verma(self,weight,L,n):
        '''转化为Verma模元素，weight为权的正交基形式，n为李代数阶数，L为相应的李代数'''
        syms = self.syms
        dim = self.L_dim
        m = (dim-n)/2
        fs,hs,es = syms[:m],syms[m:m+n],syms[m+n:] #分三部分
        element = self.element #清除零键值
        element = {key:value for key,value in zip(element.keys(),element.values()) if not self.zero_key(key,es)}
        h_matrixs = L.h_matrixs
        res = {}
        for key in element.keys():
            num = 1
            for h in key:
                if h in hs:
                    i = hs.index(h)
                    num *= L.root_on_diag(weight,h_matrixs[i])
            if num:
                new_key = tuple(i for i in key if i in fs) #去掉h部分
                if new_key in res:
                    res[new_key] += num*element[key]
                else:
                    res.update({new_key:num*element[key]})
        return self(res)
    
    @staticmethod
    def zero_key(key,es):
        '''包含e则判断为零键值'''
        for s in key:
            if s in es:
                return True
        return False
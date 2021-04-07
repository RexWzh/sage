# var('q')
class Element():
    '''sp型数据
    数据形式为字典，代表基元的线性组合。
    键值形式为元组，代表基元。
    '''
    def __init__(self,element,n,ordering=False):
        '''输入元组或字典'''
        if isinstance(element,tuple):
            element = {element:1}
        assert isinstance(element,dict),'输入格式不对！'
        syms = var('y%d '*n%tuple(range(n,0,-1))+
                   'x%d '*n%tuple(range(1,n+1))) # 基元
        self.n = n
        self.syms = syms
        if ordering:
            element = self.standardize(element)
        self.element = element
        
        
    def key2element(self,element,c=1):
        '''将单个元素或键值化标准形式，系数默认为1'''
        if not c:return {} # 零元情形
        syms = self.syms
        m = len(element) # 元素长度
        n = self.n
        coefs = [syms.index(i) for i in element] # 获取序信息
        for i in range(m-1):
            x,y = coefs[i],coefs[i+1]
            if x > y: # 发现反序
                head,tail = element[:i],element[i+2:]
                element = head+(element[i+1],element[i])+tail
                if x+y != (2*n-1): # 非反号逆序
                    elements = [(element,c*q)]
                else: # 反号逆序
                    elements = [(element,c*q*q)]
                    for j in range(y):
                        elements.append([head+(syms[j],syms[2*n-1-j])+tail,c*(q^3-q)*q^(y-1-j)])
                break
        else: # 没有逆序
            return {element:c}
        res = {} # 递归调用
        for element,c in elements:
            new = self.key2element(element,c)
            self.add_dict(res,new,replace=True)
        return res
    
    def standardize(self,element):
        '''将元素化为标准基的形式'''
        res = {}
        for key in element.keys():
            self.add_dict(res,self.key2element(key,element[key]),replace=True)
        return res
    
    @staticmethod
    def vec2element(vec,syms=None):
        '''以向量方式，输入元素'''
        n = len(vec)
        assert n%2==0,'向量长度必须为偶数'
        n //= 2
        if syms is None: # 初始化基元符号
            syms = var('y%d '*n%tuple(range(n,0,-1))+
                       'x%d '*n%tuple(range(1,n+1)))
        element = []
        for x,n in zip(syms,vec):
            element.extend([x]*n)
        return tuple(element)

    
    @staticmethod
    def add_dict(a,b,replace=False):
        '''字典相加合并'''
        if not replace: a = a.copy()
        for key in b.keys():
            a[key] = (a[key]+b[key]) if key in a else b[key]
        return a
    
    
    def identity(self):
        '''返回单位元'''
        return ELement({tuple():1},self.n)
    
    @staticmethod
    def pow_simplify(element):
        '''为便于编写，计算中，元素幂次为打开状态
        本函数将幂次重新收起，用于打印显示'''
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
    
    def __add__(self,obj):
        '''pbw元素加法'''
        assert isinstance(obj,Element),'运算必须为Element类'
        a,b = self.element,obj.element
        return Element(self.add_dict(a,b),self.n)
    
    def __sub__(self,obj):
        '''pbw元素减法'''
        assert isinstance(obj,Element),'运算必须为Element类'
        a,b = self.element,obj.element
        b = {key:-val for key,val in zip(b.keys(),b.values())}
        return Element(self.add_dict(a,b),self.n)
    
    def __neg__(self):
        '''pbw取负'''
        c = self.element
        return self({key:-val for key,val in zip(c.keys(),c.values())})
    
    def __bool__(self):
        '''判断是否为0'''
        return bool(self.element)
    
    def __eq__(self,obj):
        '''判断相等'''
        assert isinstance(obj,Element),'比较必须为Element类'
        return bool(self.element==obj.element)
    
    def __mul__(self,obj):
        '''元素乘法，支持右元素或数值
            注意：乘法不检查数据类型！！'''
        if not isinstance(obj,Element): #若obj为数
            c = {key:obj*val for key,val in zip(self.element.keys(),self.element.values())}
            return ELement(c)
        a,b = self.element,obj.element
        res = {}
        for k1 in a.keys():
            for k2 in b.keys():
                self.add_dict(res,self.key2element(k1+k2,a[k1]*b[k2]),replace=True)
        return Element(res,self.n)
    
    def __pow__(self,k):
        '''元素的次幂（仅支持正整数）'''
        assert isinstance(k,(int,Integer)),'仅支持整数次幂！'
        a = self
        for i in range(k-1):
            a = a*self
        return a
    
    def __contains__(self,key):
        '''判断该项系数是否非0'''
        return bool(key in self.element)
    
    def __getitem__(self,key):
        '''获取该项系数'''
        return self.element[key]
    
    def __str__(self):
        element = self.element
        element = {self.pow_simplify(key):val for key,val in zip(element.keys(),element.values())}
        return str(element)

    def __repr__(self):
        return self.__str__()
    
    @classmethod
    def expand(cls,element):
        '''将系数展开'''
        return {cls.pow_simplify(key):expand(element[key]) for key in element.keys()}
    
    @classmethod
    def simplify(element):
        '''将系数化简'''
        return {cls.pow_simplify(key):simplify(element[key]) for key in element.keys()}
    
    @classmethod
    def factor(element):
        '''将系数因式分解'''
        return {cls.pow_simplify(key):factor(element[key]) for key in element.keys()}
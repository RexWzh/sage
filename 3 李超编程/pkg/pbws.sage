class SPBWElement(PBWElement):
    '''李超代数的pbw基（修改两个函数）
    '''
    
    def __init__(self,element,N,syms,odds):
        '''初始元素为字典或U1上的符号变量'''
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
        self.odds = odds #奇根对应的符号元
        self.keys = element.keys()
        
    def __call__(self,obj):
        '''元素转pbw类'''
        return SPBWElement(obj,self.N,self.syms,self.odds)
        
    def tuple2dict(self,element,c=1):
        '''将符号变量元组转为字典，系数默认为1'''
        if not c:return {} #零元情形
        syms = self.syms
        odds = self.odds
        n = len(element) #元素长度
        for i in range(n-1): #奇元素情形！
            if element[i] in odds and element[i]==element[i+1]:
                return {}
        N = self.N #结构常数
        coefs = [syms.index(i) for i in element] #获取对应系数
        for i in range(n-1):
            x,y = coefs[i],coefs[i+1]
            if x > y: #发现反序
                e = element[:i] + (element[i+1],element[i]) + element[i+2:] #逆序处理
                seq = [(c*coef,sym) for coef,sym in zip(N[x][y],syms) if coef] #李括号项
                seq = [(coef,element[:i]+(sym,)+element[i+2:]) for coef,sym in seq]
                break
        else: #一切正序
            return {element:c}
        if element[i] in odds and element[i+1] in odds:
            c = -c
        res = self.tuple2dict(e,c)
        for coef,sym in seq:
            new = self.tuple2dict(sym,coef)
            res = self.add_dict(res,new)
        return res
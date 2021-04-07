class CSA():
    '''cartan 子代数
    属性
        .cartan_type Cartan类型
        .cartan_matrix Cartan阵
        .gram_matrix 格拉姆矩阵
        .D 可对称化的对角阵
        .standard_form {h_i}上的标准内积
        .ID 双线性型的ID
    方法
        .simple_roots_check() 对偶单根基
        .simple_roots         单根基
        .fundamental_weights_check 对偶基本支配权
        .fundamental_weights       基本支配权
        .highest_root_check() 对偶最高根
        .highest_root()       最高根（等同与对偶最高根）
        .weyl_vector_check() 对偶weyl向量
        .weyl_vector()       weyl向量
        .alcove_vertices() alcove的顶点
        .alcove_center()   alcove的加权中心，通过check参数可设置类型
        .simple_reflections() 单反射
        .coef_c() 系数c_i
        .coef_a() 系数a_i
        .basis() 标准基（等同于.simple_roots_check）
        .ascending_chain() 元素到alcove的偏序链
        
    待添加功能：求给定基元下的线性表示系数
    '''
    def __init__(self,cartan_type):
        s,l = self.cartan_type = cartan_type = CartanType(cartan_type)
        assert cartan_type.is_finite(),'Cartan阵必须为有限型'
        
        # Cartan矩阵
        self.cartan_matrix = A = CartanMatrix(cartan_type)
        
        # 可对称化中的对角矩阵
        self._c,self._a = None,None
        self.D =  matrix.diagonal([i/j for i,j in zip(self.coef_a(),self.coef_c())])
        
        # {h_i}上的双线性型及gram矩阵
        self.gram_matrix = A*self.D
        self.standard_form = BilinearForm(dim=l, sym='h', gram=self.gram_matrix, ID=s+'%d'%l)
        
        # 初始化数据（隐藏属性）
        self._simple_roots = None # 单根
        self._simple_roots_check = None # 对偶单根
        
        self._fw_check = None # 对偶支配权
        self._fundamental_weights = None # 支配权
        
        self._rho = None # weyl向量
        self._rho_check = None # 对偶weyl向量
        
        self._theta = None # 最长根（同时也是对偶最长根）
        
        self._reflections = None # 单反射
        self._vertices = None # 顶点集
    
    @property
    def ID(self):
        '''双线性型计算时的ID'''
        return self.standard_form.ID
    
    def basis(self):
        '''空间标准基'''
        return self.simple_roots_check()
    
    def simple_roots_check(self):
        '''H上的对偶单根'''
        if self._simple_roots_check is None:
            self._simple_roots_check = self.standard_form.basis()
        return self._simple_roots_check
    
    def fundamental_weights_check(self):
        '''H上的对偶基本支配权'''
        if self._fw_check is None:
            A = self.cartan_matrix^-1
            self._fw_check = [sum([i*vec for vec,i in zip(self.basis(),col)]) for col in A]
        return self._fw_check
    
    def simple_roots(self):
        '''H上的单根'''
        if self._simple_roots is None:
            a,c = self.coef_a(),self.coef_c()
            hs = self.simple_roots_check()
            self._simple_roots = [c_i/a_i*h for h,a_i,c_i in zip(hs,a,c)]
        return self._simple_roots
    
    def fundamental_weights(self):
        '''H上的基本支配权'''
        if self._fundamental_weights is None:
            a,c = self.coef_a(),self.coef_c()
            fw_checks = self.fundamental_weights_check()
            self._fundamental_weights = [c_i/a_i*fw for fw,a_i,c_i in zip(fw_checks,a,c)]
        return self._fundamental_weights
    
    def alcove_vertices(self):
        '''alcove 顶点'''
        if self._vertices is None:
            fw_check = self.fundamental_weights_check()
            self._vertices = [fw/i for fw,i in zip(fw_check,self.coef_a())]
        return self._vertices
    
    def weyl_vector_check(self):
        '''H上的对偶weyl向量'''
        if self._rho_check is None:
            standard = self.standard_form
            self._rho_check = standard.vec2element(sum(self.cartan_matrix^-1)) # sum按列求和
        return self._rho_check
# A-D型直接构造
#         h = self.basis()
#         if s == "A":
#             rho = 1/2*sum([(l-i+1)*i*h[i-1] for i in range(1,l+1)])
#         elif s == "B":
#             rho = 1/2*sum([i*(2*l-i+1)*h[i-1] for i in range(1,l)])
#             rho += 1/4*l*(l+1)*h[-1]
#         elif s == "C":
#             rho = 1/2*sum([i*(2*l-i)*h[i-1] for i in range(1,l+1)])
#         elif s == "D":
#             rho = 1/2*sum([i*(2*l-i-1)*h[i-1] for i in range(1,l-1)])
#             rho += 1/4*l*(l-1)*(h[-1]+h[-2])
    def weyl_vector(self):
        '''H上的weyl向量'''
        if self._rho is None:
            vec = sum(self.cartan_matrix.T^-1)
            hs = self.simple_roots()
            self._rho = sum([i*h for i,h in zip(vec,hs)])
        return self._rho
    
    def alcove_center(self,dual=True):
        '''alcove 的加权中心，dual设置加权方式'''
        if dual:
            a = sum(self.coef_a())
            return self.weyl_vector_check()/(a+1)
        c = sum(self.coef_c())
        return self.weyl_vector()/(c+1)
    
    def highest_root(self):
        '''H上的最高根'''
        if self._theta is None:
            a = self.coef_a()
            hs = self.simple_roots()
            self._theta = sum([a_i*h for a_i,h in zip(a,hs)])
        return self._theta
        
    def highest_root_check(self):
        '''H上的对偶最高根'''
        return self.highest_root()
    
    def simple_reflections(self):
        '''单反射关于{h_i}的作用矩阵'''
        if self._reflections is None:
            self._reflections = []
            s,l = self.cartan_type
            A = self.cartan_matrix
            for i in range(l): # W0
                ref  = matrix.identity(l) 
                ref[i] = [-1 if i==j else -A[j,i] for j in range(l)]
                self._reflections.append(ref)
        return self._reflections
    
    def ascending_chain(self,x0,x_seq=False):
        '''从x0出发，作用到基本域的偏序链，x_seq指定是否返回x序列'''
        assert x0.parent.ID == self.ID,'数据类型不匹配！'
        refs = self.simple_reflections() # 单反射
        chain = [x0]
        actions = []
        while True:
            x0 = chain[-1] # 末端元素
            for i,r in enumerate(refs):
                x = r*x0
                if x>x0:
                    actions.append(i+1)
                    chain.append(x)
                    break
            else: # 没有发生break
                break
        return [chain,actions] if x_seq else actions
    
    def __str__(self):
        s,l = self.cartan_type
        return "Cartan subalgebra of type %s%d"%(s,l)
    def __repr__(self):
        return self.__str__()
    
    def coef_a(self):
        '''系数a_i'''
        if self._a is None:
            self._a = RootSystem(self.cartan_type).root_space().highest_root().coefficients()
        return self._a
        return RootSystem(cartan_type).root_space().highest_root().coefficients()
    
    def coef_c(self):
        '''系数c_i'''
        if self._c is None:
            s,l = self.cartan_type
            if s in "ADE":
                self._c = self.coef_a()
            elif s == "C":
                self._c = [1] * l
            elif s == "B":
                self._c = [2] * l
                self._c[0],self._c[-1] = 1,1
            elif s == "F":
                self._c = [2,3,2,1]
            elif s == "G":
                self._c = [1,2]
        return self._c
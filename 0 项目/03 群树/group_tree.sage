import time

def permutation_mat(sigma):
    '''输入置换序列，返回置换矩阵
    sigma为序列，标号可含0或不含
    length指定长度，指定后，标号从0开始，'''
    n = len(sigma)
    if 0 not in sigma: # 移位
        sigma = [i-1 for i in sigma]
    assert set(sigma)==set(range(n)),'必须为从0或1开始的序列'
    mat = matrix.zero(n)
    for i,j in enumerate(sigma):
        mat[j,i] = 1
    return mat


class GroupTree():
    '''群生成树，生成元为矩阵(哈希值用了myhash)'''
    def __init__(self, gens, length=-1, index=None, elements=True):
        '''输入信息：
        gens 为矩阵生成元，可逆且不能为空
        length 设置元素长度，默认-1生成整树，若为无限群将死循环
        index 设置生成元名称
        elements 默认不记录元素
        
        类属性：
        .n 矩阵列数
        .unit 单位矩阵
        .depth 当前遍历深度
        .bottom 当前遍历的最底层（列表）
        .elements 当前得到的群元素（列表）
        .expressions 群元关于生成元的表达元组，取字典序下最小（字典）
        '''
        assert len(gens),'生成元集不能为空！'
        assert all([m.is_invertible() for m in gens]),'生成元矩阵必须可逆！'
        if index is None: # 初始化生成元名称
            index = tuple([i for i in range(len(gens))])
        self.n = n = gens[0].ncols() # 矩阵列数
        self.unit = unit = matrix.identity(n) # 单位阵
        self.bottom = [unit] # 当前树底
        self.depth = 0 # 当前深度
        self.expressions = {self.myhash(unit):tuple()} # 群元关于生成元的表达（只给出字典序下最小的一种）
        self.elements = [unit] if elements else None
        self.gens,self.length,self.index = gens,length,index
        self._bool = False
        self.search(length) # 遍历群元素
    
    def search(self,length):
        '''遍历至 depth>=length 或达到最长'''
        if 0 <= length <= self.depth: # 已经达到所需深度
            return True
        upper = self.bottom # 当前最深一层
        gens,index = self.gens,self.index #生成元及名称
        exprs,elements = self.expressions,self.elements
        while length-self.depth: # 深度不够
            lower = [] # 开辟一层矩阵
            new_exprs = {} # 开辟一层表达
            for m1 in upper:
                h1 = self.myhash(m1)
                for i,m2 in zip(index,gens):
                    m = m1 * m2
                    hm = self.myhash(m) # 新矩阵哈希值
                    if (hm in exprs) or (hm in new_exprs): # 跳过已有元素
                        continue
                    expr = exprs[h1] + (i,) # 新元素的表达式
                    new_exprs[hm] = expr # 添加群表
                    lower.append(m)
            if not len(lower): # 没有新元素生成
                self._bool = True #已经结束
                break
            self.depth += 1 # 添加了一层
            if elements is not None:
                elements.extend(lower)
            upper = lower # 放到最高层
            exprs.update(new_exprs)
            self.bottom = upper
        return True
    
    def __iter__(self):
        '''设置遍历，直至得到终点'''
        return self
    
    def __next__(self):
        '''往下遍历一阶，遍历结束返回-1，否则返回深度'''
        if self._bool:
            raise StopIteration
        self.search(self.depth+1)
        return self.depth
    
    @property
    def card(self):
        '''返回当前群的元素个数'''
        return len(self.expressions)
    
    @staticmethod
    def myhash(mat):
        '''矩阵展平哈希化'''
        n = mat.nrows()
        mat = mat + 11*matrix.ones(n)
        res = tuple(i for line in mat for i in line)
        return hash(res)
    
    def __bool__(self):
        '''是否完整生成树（只有 research 确认过才会返回True）'''
        return self._bool
    
    def __str__(self):
        c = self.card
        if c == 1:
            return 'identity group'
        return 'Group of %d elements'%c
    def __repr__(self):
        return self.__str__()
        
def transform_mat(cartan_type,dx=0):
    '''B,C,D,E族的平移矩阵，D族有两个元素，通过dx确定，取l,l-1或0'''
    s,l = CartanType(cartan_type)
    assert s not in ['F',"G"],'仅定义了A-E族'
    Tx = matrix.identity(l+1) # 平移阵
    if s =="A":
        Tx[0] = [2] + [1]*l
        Tx[-1] = [-1]*l + [0]
    if s == "B":
        assert l>=3,'B族至少为3阶'
        Tx[0] = [ 2,1]+[ 2]*(l-2)+[ 1]
        Tx[1] = [-1,0]+[-2]*(l-2)+[-1]
    if s == "C":
        assert l>=2,'C族至少为2阶'
        return transform_mat(cartan_type=["A",l]) # 二者平移阵相同
    if s == "D":
        assert l>=4,'D族至少为4阶'
        assert dx in [l,l-1,0,-1,-2]
        if dx==0: # 平移元素x
            Tx[l-1] = [-1,-1]+[-2]*(l-3)+[0,-1]
            Tx[l]   = [ 1, 1]+[ 2]*(l-3)+[1, 2]
        else:
            Tx[0] =  [ 2, 1]+[ 2]*(l-3)+[ 1, 1]
            Tx[dx] = [-1,-1]+[-2]*(l-3)+([-1, 0] if (dx in [-1,l]) else [0,-1])
    if s == "E":
        assert 6<=l<=7,'E族必须为6或7阶'
        if l==6:
            Tx[0] = [2,1,2,2,3,2,1]
            Tx[6] = [-1,-1,-2,-2,-3,-2,0]
        else:
            Tx[0] = [2,2,2,3,4,3,2,1]
            Tx[7] = [-1,-2,-2,-3,-4,-3,-2,0]
    return Tx.T

# 调试代码
# Lie = [["A",range(1,9)],["D",range(4,9)],["E",range(6,8)],
#        ["B",range(3,9)],["C",range(2,9)]]
# for s,l in Lie:
#     for n in l:
#         print(s+'%d'%n)
#         print(transform_mat(cartan_type=[s,n]))
    
def affine_graph_automorphism(cartan_type):
    '''仿射1型的图自同构群'''
    s,l = CartanType(cartan_type)
    if (s,l) in (("E",8),("G",2),("F",4)): # 平凡群
        return GroupTree([matrix.identity(l+1)])
    
    if s == "A":
        rotate = permutation_mat(list(range(1,l+1))+[0]) # 旋转元
        if l == 1:
            return GroupTree([rotate])
        ref = permutation_mat(list(range(l+1)[::-1]))# 反射元
        return GroupTree([rotate,ref])
    
    if s == "B":
        assert l>=3,'B族至少为3阶'
        return GroupTree([permutation_mat([1,0]+list(range(2,l+1)))])
    
    if s == "C":
        assert l>=2,'C族至少为2阶'
        ref = permutation_mat(list(range(l+1)[::-1]))
        return GroupTree([ref])
    
    if s == "D":
        assert l>=4,'D族至少为4阶'
        s1 = permutation_mat([1,0]+list(range(2,l+1)))
        if l==4:
            s2 = permutation_mat((1,3,2,4,0))
            return GroupTree([s1,s2])
        s2 = permutation_mat(list(range(l-1))+[l,l-1])
        s3 = permutation_mat(list(range(l+1))[::-1])
        return GroupTree([s1,s2,s3])
    
    if s == "E":
        assert l in [6,7],'E族必须为6-8阶'
        if l==6:
            s1 = permutation_mat((6,0,5,2,4,3,1))
            s2 = permutation_mat((1,0,3,2,4,5,6))
            return GroupTree([s1,s2])
        s1 = permutation_mat((7,6,2,5,4,3,1,0))
        return GroupTree([s1])
        
def save_vari(var,name):
    '''Python数据存储'''
    import pickle
    fout = open(name,'wb')
    pickle.dump(var,fout)
    fout.close()
    
def read_vari(name):
    '''Python数据读取'''
    import pickle
    fin = open(name,'rb')
    a = pickle.load(fin)
    fin.close()
    return a
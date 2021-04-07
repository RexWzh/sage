'''文件包括下列函数
    line2root(line,positive_roots)  线段型解转符号根（有序，不合并系数）
    dec2root(dec,positive_roots)   正根分解转符号根（无序，合并系数）
    highest_root(cartan_type)     最高根关于单根的系数
    positive_decomposition_recursive(cartan_type,lam=None) 基于线段切割的正根分解算法（递归算法，针对A-D的线段情形）
    positive_decomposition(M,lam)   基于正根性质的分解算法（BFS搜索算法，针对一般情形）
    group_table_product(M,elements)  乘法表
    root_system_matrix(cartan_type,non_root=0) 根系信息矩阵
'''

load('kostant.sage')
# 将线段型解转化为符号根形式
line2root = lambda line,positive_roots:tuple(positive_roots[i] for i in line)
# 正根分解转符号根形式
dec2root = lambda dec,positive_roots:tuple(i*root for i,root in zip(dec,positive_roots) if i)
highest_root = lambda cartan_type: RootSystem(cartan_type).root_space().highest_root().coefficients()


def positive_decomposition_recursive(cartan_type, lam=None, line=None, M=None):
    '''基于线段切割的正根分解算法，通过递归实现
基本使用
    输入必要参数 cartan_type
    lam 默认取最高根，可赋其他值
    
参数信息
    cartan_type Cartan类型，lam指定时，阶数l尽量取小
    lam 需要分解的向量(s,k)，默认为最高根
    line 解线段
    M 根系信息矩阵，由 cartan_type 参数自动生成
    
注意事项
    line 与 M 为递归过程变量，外部调用时不应修改
    cartan_type 的参数 l
        l 用于初始化矩阵 M ，最高根 lam ，以及解线段 line
        l 的取值会影响整体的标号信息，不影响解的总数和性质
        当这三个变量都指定时，参数 l 将被忽略
    标号沿用内置的单根系标号（递推标号和内置标号相反，所以讨论要分类）
Ps：算法主要针对最高根，很多内容可以不用这么麻烦地分类讨论
    '''
    # 检查输入
    s,l = cartan_type = CartanType(cartan_type)
    assert s in "ABCD",'仅支持A-D型计算'
    
    #### 初始化数据 ####
    if lam is None:# 初始化lam（默认作为最长根）
        lam = {"A":(l,0), "B":(l,l-1), "C":(l-1,l-1), "D":(l-2,l-3)}[s]
    left,right = lam
    assert left>=right>=0,'左侧下标需大于右侧，且下标非负'
    
    if M is None: # 初始化根系矩阵
        M = root_system_matrix(cartan_type,non_root=None)
    n = M.ncols() # 正根数目
    
    if line is None: # 初始化解线段
        line_left = tuple(i for i in range(l-left,l))
        line_right = tuple(i for i in range(l-right,l))[::-1]
        if s=="C":
            line_left = (l-left-1,) + line_left
            line_right = line_right[1:] + (l-right-1,)
        elif s=="D": # D 型
            line_left = tuple(i for i in range(l-left-2,l))
            line_right = tuple(i for i in range(l-right-2,l-2))[::-1]
        line = line_left + line_right
    
    #### 递归算法 ####
    sols = []
    # A型
    if s=="A":
        assert 0 in lam,'A型必须有一侧为0'
        if 1 in lam: return [line]
        # 拆末尾，补末尾
        new_line = line[:-1] 
        for new_line in  positive_decomposition_recursive(\
                          cartan_type, lam=(left-1,0), M=M, line=new_line):
            sols.append(new_line+(line[-1],)) 
        # 绑末尾
        i,j = line[-2:]
        new_line = line[:-2] + (M[i][j],) 
        for new_line in  positive_decomposition_recursive(\
                          cartan_type, lam=(left-1,0), M=M, line=new_line):
            sols.append(new_line)
        return sols
    # BCD 型
    adj = 0 if s=="B" else 1 # BCD型公式相近，adj为索引调整参数
    if right==0: # 初值
        if s in "BC":
            return positive_decomposition_recursive(["A",l], lam=(left+adj,0), M=M, line=line)
        l1 = line[:-3] + (group_table_product(M,line[-3:]),)
        i,j,k = line[-1],line[-2],line[-3]
        l2,l3,l4 = line[:-3]+(M[i,k],), line[:-3]+(M[j,k],), line[:-2]
        for new_line in positive_decomposition_recursive(\
                  ["A",l], lam=(left,0), M=M, line=l1):
            sols.append(new_line)
        for new_line in positive_decomposition_recursive(\
                  ["A",l], lam=(left,0), M=M, line=l2):
            sols.append(new_line+(j,))
        for new_line in positive_decomposition_recursive(\
                  ["A",l], lam=(left,0), M=M, line=l3):
            sols.append(new_line+(i,))
        for new_line in positive_decomposition_recursive(\
                  ["A",l], lam=(left,0), M=M, line=l4):
            sols.append(new_line+(i,j))
        return sols
    if left==right:
        if left==1: # 初值
            if s=="B": return [line]
            if s=="C":
                i,j = line[:2]
                sols = [line,(M[i,j],i),(group_table_product(M,line),)]
                return sols
            a,b,c = line[:-1] # D 型
            sols =  [line,(a,b,M[a][c]),(a,c,M[a,b]),(M[a,b],M[a,c]),
                          (a,group_table_product(M,[a,b,c])),]
            return sols
        # 首尾两处绑定
        i1,i2 = line[:2]
        j1,j2 = line[-2:]
        new_line = (M[i1][i2],)+line[2:-2]+(M[j1][j2],)
        for new_line in positive_decomposition_recursive(\
                      cartan_type, lam=(left-1,right-1), M=M, line=new_line):
            sols.append(new_line)
        # 拆开右侧
        new_line = line[:-1]
        for new_line in positive_decomposition_recursive(\
                      cartan_type, lam=(left,right-1), M=M, line=new_line):
            sols.append(new_line+(line[-1],))
        return sols
    if left>right+1:
        # 绑定左侧两根
        i1,i2 = line[:2]
        new_line = (M[i1][i2],)+line[2:]
        for new_line in positive_decomposition_recursive(\
                      cartan_type, lam=(left-1,right), M=M, line=new_line):
            sols.append(new_line)
        # 拆开左侧
        new_line = line[1:]
        for new_line in positive_decomposition_recursive(\
                      cartan_type, lam=(left-1,right), M=M, line=new_line):
            sols.append((line[0],)+new_line)
        return sols
    if left==right+1: # left=right+1
        # 逐级绑定左侧根
        for i in range(1,left+adj): # 左侧逐级打开
            ele = group_table_product(M,line[:i]) # 绑定左侧i个根
            new_line = line[i:][::-1] # 反序，使得左侧保持最长
            for new_line in positive_decomposition_recursive(\
                          cartan_type, lam=(right,left-i), M=M, line=new_line):
                sols.append((ele,)+new_line[::-1])
        if s in "BC":
            ele = group_table_product(M,line[:left+adj])
            new_line = (ele,) + line[left+adj:]
            sols.extend(positive_decomposition_recursive(["A",l], lam=(right+1,0), M=M, line=new_line))
            return sols
        e0 = group_table_product(M,line[:left]) # 左侧所有根之和
        i,j = line[left:left+2] # - + 根
        e1,e2 = M[e0,i],M[e0,j]
        l1,l2,l3 = line[left+1:], (line[left],)+line[left+2:], (M[e1,j],)+line[left+2:]
        for new_line in positive_decomposition_recursive(["A",l],
                                    lam=(right+1,0), M=M, line=l3[::-1]):
            sols.append(new_line[::-1]) # 绑定左侧元素
        for new_line in positive_decomposition_recursive(["A",l],
                                    lam=(right+1,0), M=M, line=l1[::-1]):
            sols.append((e1,)+new_line[::-1])
        for new_line in positive_decomposition_recursive(["A",l],
                                    lam=(right+1,0), M=M, line=l2[::-1]):
            sols.append((e2,)+new_line[::-1])
        return sols


def group_table_product(M,elements):
    '''通过群乘法表给出的元素乘积，elements 由索引构成
    mij = mi*mj
    '''
    assert len(elements),'乘积个数不能为0'
    if len(elements)==1:return elements[0]
    M = matrix(M)
    ele = elements[0]
    for i in elements[1:]:
        ele = M[ele][i]
    return ele

def root_system_matrix(cartan_type,non_root=0):
    r'''记录根系信息的矩阵，前l个索引对应为单根
    M[i][j] = k        if \beta_i+\beta_j=\beta_k \in \Fai
            = non_root if \beta_i+\beta_j not \in \Fai
    non_root 默认取0，注意到单根不能分解为正根和，取0与索引0不会有歧义
    non_root 设置为 None 时，取值将是 M.cols，即溢出索引，调试时可使用。
    '''
    s,l = cartan_type = CartanType(cartan_type)
    positive_roots = list(RootSystem(cartan_type).root_space().positive_roots())
    n = len(positive_roots)
    if non_root == None:
        non_root = n
    M = matrix.zero(n) # 信息矩阵
    for i,ri in enumerate(positive_roots):
        for j,rj in enumerate(positive_roots):
            if j < i: 
                M[i,j] = M[j,i]
            elif ri+rj not in positive_roots:
                M[i,j] = non_root
            else:
                M[i,j] = positive_roots.index(ri+rj)
    return M
"""测试代码
cartan_type = "B3"
print(root_system_matrix(cartan_type))
print(list(RootSystem(cartan_type).root_space().positive_roots()))""";


def positive_decomposition(M,lam):
    '''生成向量lam对应的正根分解，lam的维数为单根的个数
    M 根系信息的矩阵，非根取数值0
    lam 关于单根集的系数，数据形式为元组
    '''
    lam = tuple(lam) # 元组形式
    M = matrix(M) # 数据为矩阵形式
    n,m = M.ncols(),M.nrows()
    assert all([i>=0 for i in lam]),'向量系数不能为负！'
    assert n==m>=len(lam),'矩阵维数不能小于向量维数'
    
    sol = lam + (0,)*(n-len(lam)) # 初始解
    sols = [[sol]] # 解集
    sols_hash = {hash(sol)} # 哈希表
    while True:
        new = [] # 新一层解
        for root in sols[-1]: # 上层解
            for i in range(n): # 抓一根
                if not root[i]:continue
                for j in range(i,n): # 再抓一根
                    if any([not root[j], not M[i][j], i==j and root[i]==1]):
                        continue # 跳过三种情形：系数0，相加非根，同一根相加但系数1
                    k = M[i][j]
                    sol_new = list(root)
                    sol_new[i] -= 1
                    sol_new[j] -= 1
                    sol_new[k] += 1
                    sol_new = tuple(sol_new)
                    sol_hash = hash(sol_new)
                    if sol_hash not in sols_hash:
                        sols_hash.add(sol_hash)
                        new.append(sol_new)
        if not len(new): break # 没有新解
        sols.append(new)
    return sols
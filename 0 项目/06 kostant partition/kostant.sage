def kostant_partition_fun(cartan_type, lam=None, method='formula', reduced=True):
    '''A-D型kostant partition函数的计算，支持递归法和公式法，默认使用公式法
    reduced 为 False 时，输入非简约根系'''
    # 检查输入，指定lam时，l将被忽略
    s,l = cartan_type = CartanType(cartan_type)
    assert s in "ABCD",'仅支持A-D型计算'
    assert reduced or s=='B','非简约型只有B型'
    
    if lam is None: # 初始化lam（默认作为最长根）
        lam = {"A":(l,0), "B":(l,l-1), "C":(l-1,l-1), "D":(l-2,l-3)}[s] if reduced else (l,l)
    r,k = max(lam),min(lam)
    assert k>=0,'数列下标不能为负'
    
    # 数列初值 k=0
    if k == 0:
        if s == "A":
            assert r!=0,'A型 r 必须为正整数'
            return 2^(r-1)
        if s == "B":
            if r>0: return 2^(r-1)
            return 0 if reduced else 1
        if s == "C": return 2^r
        if s == "D": return 2^(r+1) if r else 1
    assert s!="A",'A型另一参数必须为0'
    
    # 递归算法
    if method == 'recursive':
        if r>k+1:
            return 2^(r-k-1)*kostant_partition_fun(cartan_type,lam=(k+1,k),method='recursive',reduced=reduced)
        if r==k: # r=k情形
            res = kostant_partition_fun(cartan_type,lam=(k-1,k-1),method='recursive',reduced=reduced)
            return res + kostant_partition_fun(cartan_type,lam=(k,k-1),method='recursive',reduced=reduced)
        # r=k+1情形
        res = 0
        for i in range(r):
            res += kostant_partition_fun(cartan_type,lam=(k,i),method='recursive',reduced=reduced)
        res += {"B":2^(k-1), "C":2^k, "D":3*2^k}[s]
        return res
    
    # 构造性算法 (k>0)
    s5,x = sqrt(5),var('x')
    expr_c = (((5-s5)/2)^(x+1) + ((5+s5)/2)^(x+1))/5
    expr_d = (1-s5)/2*((5-s5)/2)^x + (1+s5)/2*((5+s5)/2)^x
    if reduced:
        expr = {"C":expr_c,"D":expr_d,"B":expr_d/5}[s]
        if r==k: return int(expr.subs(x=k))
        res = int(expr.subs(x=k+1)-expr.subs(x=k)) # C_{k+1,k}的值
        return 2^(r-k-1)*res
    expr = expr_c*2 # BC型根系
    if r==k:return int(expr.subs(x=k-1))
    res = int(expr.subs(x=k)-expr.subs(x=k-1))
    return 2^(r-k-1)*res
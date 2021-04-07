
def Coordination(base,vec,field=QQ):
    """求vec关于base的坐标"""
    V = VectorSpace(field,vec.degree()) #生成向量空间
    U = V.subspace_with_basis(base) #以base为基
    return U.coordinate_vector(vec)

def alpha_chain(alpha,beta,roots):
    """计算beta上alpha-chain的p,q值"""
    p,q = 0,0
    root = beta + alpha
    while root in roots: #计算q值
        root += alpha
        q += 1
    root = beta - alpha
    while root in roots:
        root -= alpha
        p += 1
    return p,q
### 考虑计算效率，使用gap进行交互 ###
LieExp := function(mat) 
    ## 求导子的指数函数,第二参数默认取0 ##
    local i,exp_mat,new_mat;
    i := 1;
    exp_mat := IdentityMat(Size(mat),Rationals);
    new_mat := mat;
    while not IsZero(new_mat) do #循环至添加矩阵为0
        exp_mat := exp_mat + 1/Factorial(i)*new_mat;
        i := i + 1;
        new_mat := mat^i;
    od;
    return exp_mat;
end;;

LieThetas := function(s,n)
    ## 求thetai=exp(adei)*exp(-adfi)*exp(adei) 的矩阵表达 ##
    local L,b,m,ade,adf,theta;
    L := SimpleLieAlgebra(s,n,Rationals); #新建李代数
    b := Basis(L); #获取李代数一组基
    ade := List([1..n],i->AdjointMatrix(b,b[i])); #初始化ei
    m := (Size(b)-n)/2; #正根总数
    adf := List([1..n],i->AdjointMatrix(b,b[m+i])); #初始化fi
    theta := List([1..n],i->LieExp(ade[i])*LieExp(-adf[i])*LieExp(ade[i])); #初始化thetai
    return theta;
end;
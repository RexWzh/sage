import os
import subprocess as sp
Join = lambda *l:'\\'.join([str(i) for i in l])
def latex2png(txt, filename='rex', path='.', density=None,
              packages=(), outputtex=False, outputpdf=False):
    '''latex代码转png
    1.filename 设置输出文件名，后缀自动过滤
    2.path 设置保存路径
    3.density 设置图片分辨率
    4.packages 新增宏包
    5.tex 是否同时生成tex
    6.pdf 是否同时生成pdf
    需要函数库：subprocess as sp
    需要函数：Join, writetpp
    '''
    assert txt,'输入不能为空！'
    # 初始化参数
    filename = filename.split('.')[0]
    png = filename + '.png'
    pdf = filename + '.pdf'
    tex = filename + '.tex'
    packages = ['\\usepackage{%s}'%i for i in packages]
    packages = '\n'.join(packages)
    latex_main = r"""\documentclass[UTF8]{ctexart}
\pagestyle{empty}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{tikz}
%s
\begin{document}
%s
\end{document}"""%(packages,txt)
    
    # 生成文件
    try:
        writetpp(latex_main,filename,path,density)
    except:
        print(latex_main)
        raise Exception('tex内容编写有误')
        
    # 移动文件
    if not outputpdf:
        os.remove(Join(path,filename)+'.pdf')
    if not outputtex:
        os.remove(Join(path,filename)+'.tex')
    return True


def writetpp(txt, filename='rex', path='.', density=None, remove=True):
    '''将 txt 文本写入 tex 文件，并导出 pdf 和 png 文件
    1.filename 设置输出文件名，后缀自动过滤
    2.path 设置文件导出路径
    3.density 设置图片分辨率
    4.remove 删除 xelatex 编译时，产生的 .aux, .log 文件
    需要函数库：os（remove指定为False时不需要）, subprocess as sp
    需要函数：Join
    '''
    # 参数初始化
    filename = filename.split('.')[0] # 去后缀
    tex = filename + '.tex' 
    pdf = filename + '.pdf'
    png = filename + '.png'
    density = '' if (density is None) else '-density %s'%density # 设置分辨率
    
    # 写入tex
    with open(Join(path,tex), 'w', encoding='utf-8') as f:
        f.write(txt)
        
    # 写入pdf
    sp.check_output('xelatex %s'%tex,cwd=path, shell=True)
    
    # 写入png
    sp.check_output('magick convert  -trim +repage {density} {pdf} {png}'\
                    .format(density=density, pdf=pdf, png=png), cwd=path, shell=True)
    
    # 删除多余文件
    if remove:
        os.remove(Join(path,filename)+'.aux')
        os.remove(Join(path,filename)+'.log')
    return True
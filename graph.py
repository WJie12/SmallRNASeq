import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            argv = argv[1:]
    return opts


def main():
    args = sys.argv[1:]

    try:
        opts = getopts(args)
    except IndexError:
        print("Usage:")
        print(" -i        Input file")
        print(" -o        Output file")
        return 0

    output_path = opts.get("-o")
    if output_path is None:
        print("No output file specified.")
        return -1

    input_path = opts.get("-i")
    if input_path is None:
        print("No input file specified.")
        return -2

    rs_csv = pd.read_csv(input_path + 'y_rnadb_y_genome_result.csv')

    count_list = []
    rna_dbs = ['mature-hg19-tRNAs', 'hg19-tRNAs', 'human_rRNA_5.8S', 'human_rRNA_5S', 'human_rRNA_12S', 'human_rRNA_16S', 'human_rRNA_18S',
               'human_rRNA_28S', 'human_rRNA_45S', 'human_rRNA_other', 'miRBase_21-hsa', 'piR_human','Rfam-12.3-human']
    for idx in range(0, len(rna_dbs)):
        no_nan = rs_csv[~rs_csv[rna_dbs[idx]].isin(['*', 'NA'])]
        no_nan = no_nan.dropna()
        no_nan_name = no_nan['name']
        count_list.append(no_nan_name.astype(int).sum())
    df = pd.DataFrame(columns=('db', 'count'))
    df['db'] = ['tRNA', 'rRNA', 'miRNA', 'piRNA', 'Rfam']
    df['count'] = [sum(count_list[0:1]), sum(count_list[2:9]), count_list[10], count_list[11], count_list[12]]

    # 调节图形大小，宽，高
    plt.figure(figsize=(12, 6))
    # 定义饼状图的标签，标签是列表
    labels = df['db']
    # 每个标签占多大，会自动去算百分比
    sizes = df['count']
    # colors = ['red','yellowgreen','lightskyblue']
    # 将某部分爆炸出来， 使用括号，数值的大小是分割出来的与其他两块的间隙
    explode = (0, 0.05, 0, 0.05, 0)
    colors = cm.rainbow(np.arange(len(sizes)) / len(sizes))

    patches, l_text, p_text = plt.pie(sizes, explode=explode, colors=colors,
                                      labeldistance=1.2, autopct='%3.2f%%', shadow=False,
                                      startangle=90, pctdistance=1.1)

    # labeldistance，文本的位置离远点有多远，1.1指1.1倍半径的位置
    # autopct，圆里面的文本格式，%3.1f%%表示小数有三位，整数有一位的浮点数
    # shadow，饼是否有阴影
    # startangle，起始角度，0，表示从0开始逆时针转，为第一块。一般选择从90度开始比较好看
    # pctdistance，百分比的text离圆心的距离
    # patches, l_texts, p_texts，为了得到饼图的返回值，p_texts饼图内部文本的，l_texts饼图外label的文本

    # 改变文本的大小
    # 方法是把每一个text遍历。调用set_size方法设置它的属性
    for t in l_text:
        t.set_size(7)
    for t in p_text:
        t.set_size(9)
    # 设置x，y轴刻度一致，这样饼图才能是圆的
    plt.axis('equal')
    plt.legend(labels=labels)
    plt.title("small RNA database mapping result")
    plt.savefig(output_path + 'y_rnadb_y_genome_result.png')
    plt.show()
    plt.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
从deepBGC输出的GBK文件中提取BGC区域的CDS蛋白序列

功能：
1. 解析GBK文件中的BGC区域
2. 提取BGC区域内的CDS蛋白序列
3. 输出为FASTA格式

用法：
    python extract_deepbgc.py <input_gbk> <output_faa>
"""

import sys
from Bio import SeqIO


def extract_bgc_proteins(gbk_file: str, output_faa: str) -> dict:
    """
    从GBK文件中提取BGC区域的CDS蛋白序列

    Args:
        gbk_file: 输入GBK文件
        output_faa: 输出FASTA文件

    Returns:
        dict: 统计信息
    """
    extracted_count = 0
    bgc_region_count = 0

    try:
        with open(output_faa, "w") as out_f:
            for record in SeqIO.parse(gbk_file, "genbank"):
                bgc_regions = []

                # 第一轮遍历：精准定位被工具认定为候选 BGC 的区间
                for feature in record.features:
                    if feature.type.lower() in ['cluster', 'region', 'bgc', 'protocluster', 'candidate']:
                        bgc_regions.append((int(feature.location.start), int(feature.location.end)))
                        bgc_region_count += 1

                # 如果这条序列里压根没有 BGC，直接跳过当前 record 以提升速度
                if not bgc_regions:
                    continue

                # 第二轮遍历：寻找落在 BGC 区间内的 CDS
                for feature in record.features:
                    if feature.type == 'CDS':
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        in_bgc = False

                        # 核心修正：只依赖严谨的物理坐标交集判断
                        for bgc_start, bgc_end in bgc_regions:
                            # 只要有交集（哪怕是部分重叠）就认为该基因属于此 BGC
                            if max(start, bgc_start) <= min(end, bgc_end):
                                in_bgc = True
                                break

                        if in_bgc:
                            # 获取基因名，把任何空格强行转为下划线
                            raw_id = feature.qualifiers.get('locus_tag', [f"{record.id}_{start}"])[0]
                            gene_id = str(raw_id).replace(" ", "_")

                            translation = ""
                            if 'translation' in feature.qualifiers:
                                translation = feature.qualifiers['translation'][0]
                            else:
                                try:
                                    seq_dna = feature.extract(record.seq)
                                    translation = str(seq_dna.translate(to_stop=True))
                                except Exception:
                                    pass

                            if translation:
                                out_f.write(f">{gene_id} {record.id}:{start}-{end}\n{translation}\n")
                                extracted_count += 1

        return {
            'bgc_regions': bgc_region_count,
            'extracted': extracted_count,
            'success': True
        }

    except Exception as e:
        return {
            'error': str(e),
            'success': False
        }


def main():
    if len(sys.argv) != 3:
        print("用法: python extract_deepbgc.py <input_gbk> <output_faa>")
        print("  input_gbk:  deepBGC输出的GBK文件")
        print("  output_faa: 输出蛋白序列FASTA文件")
        sys.exit(1)

    input_gbk = sys.argv[1]
    output_faa = sys.argv[2]

    if not input_gbk.endswith('.gbk'):
        print("警告: 输入文件可能不是GBK格式")

    result = extract_bgc_proteins(input_gbk, output_faa)

    if result['success']:
        print(f"扫描结束！共发现 {result['bgc_regions']} 个 BGC 区域，成功提取 BGC 蛋白序列: {result['extracted']} 条")
        sys.exit(0)
    else:
        print(f"运行出错: {result['error']}")
        sys.exit(1)


if __name__ == "__main__":
    main()
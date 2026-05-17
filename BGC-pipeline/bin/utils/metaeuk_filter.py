#!/usr/bin/env python3
"""
处理metaEUK输出结果

功能：
1. 从metaEUK输出的GFF3文件中提取去重的CDS
2. 去重规则：染色体名 + 起始位置 + 终止位置 + 链
3. 生成可用于antiSMASH的GFF3文件

用法：
    python metaeuk_filter.py <input_gff> <output_gff>
"""

import sys
from pathlib import Path


def deduplicate_cds_gff(input_gff: str, output_gff: str) -> dict:
    """
    去重CDS GFF3文件

    Args:
        input_gff: 输入GFF3文件
        output_gff: 输出GFF3文件

    Returns:
        dict: 统计信息
    """
    input_path = Path(input_gff)
    output_path = Path(output_gff)

    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_gff}")

    seen = set()
    stats = {
        'total_lines': 0,
        'cds_lines': 0,
        'header_lines': 0,
        'unique_cds': 0,
        'duplicate_cds': 0,
        'other_lines': 0
    }

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            stats['total_lines'] += 1

            # 保留注释行
            if line.startswith('#'):
                outfile.write(line)
                stats['header_lines'] += 1
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]

            if feature_type == 'CDS':
                stats['cds_lines'] += 1
                # 去重key: 染色体 + 起始 + 终止 + 链
                key = f"{chrom}_{start}_{end}_{strand}"

                if key not in seen:
                    seen.add(key)
                    outfile.write(line)
                    stats['unique_cds'] += 1
                else:
                    stats['duplicate_cds'] += 1
            else:
                # 保留非CDS特征
                outfile.write(line)
                stats['other_lines'] += 1

    return stats


def main():
    if len(sys.argv) != 3:
        print("用法: python metaeuk_filter.py <input_gff> <output_gff>")
        print("  input_gff:  metaEUK输出的GFF3文件")
        print("  output_gff: 去重后的GFF3文件")
        sys.exit(1)

    input_gff = sys.argv[1]
    output_gff = sys.argv[2]

    print("=" * 50)
    print("处理metaEUK GFF3文件...")
    print("=" * 50)

    try:
        stats = deduplicate_cds_gff(input_gff, output_gff)

        print("\n处理完成！")
        print(f"总行数: {stats['total_lines']}")
        print(f"注释行: {stats['header_lines']}")
        print(f"CDS行数: {stats['cds_lines']}")
        print(f"  唯一CDS: {stats['unique_cds']}")
        print(f"  重复CDS: {stats['duplicate_cds']}")
        print(f"其他特征: {stats['other_lines']}")
        print(f"\n输出文件: {output_gff}")

    except FileNotFoundError as e:
        print(f"错误: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"处理失败: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
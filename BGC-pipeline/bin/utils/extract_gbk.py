#!/usr/bin/env python3
"""
提取antiSMASH输出的region gbk文件

功能：
1. 从antiSMASH输出目录提取所有 *.region*.gbk 文件
2. 添加样本名前缀防止同名覆盖
3. 统计提取的文件数量

用法：
    python extract_gbk.py <input_dir> <output_dir>
"""

import os
import sys
import shutil
from pathlib import Path


def extract_region_gbks(input_dir: str, output_dir: str) -> dict:
    """
    提取region gbk文件

    Args:
        input_dir: antiSMASH输出目录
        output_dir: 提取文件的目标目录

    Returns:
        dict: 统计信息
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # 创建输出目录
    output_path.mkdir(parents=True, exist_ok=True)

    stats = {
        'total_samples': 0,
        'total_regions': 0,
        'samples_with_regions': [],
        'samples_without_regions': []
    }

    # 遍历所有样本目录
    for sample_dir in input_path.iterdir():
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name
        region_count = 0

        # 查找region gbk文件
        for gbk_file in sample_dir.glob("*.region*.gbk"):
            if gbk_file.is_file():
                original_name = gbk_file.name
                new_name = f"{sample_name}_{original_name}"

                # 复制文件
                dest_file = output_path / new_name
                shutil.copy2(gbk_file, dest_file)

                region_count += 1
                stats['total_regions'] += 1

        stats['total_samples'] += 1

        if region_count > 0:
            stats['samples_with_regions'].append({
                'sample': sample_name,
                'regions': region_count
            })
        else:
            stats['samples_without_regions'].append(sample_name)

    return stats


def main():
    if len(sys.argv) != 3:
        print("用法: python extract_gbk.py <input_dir> <output_dir>")
        print("  input_dir:  antiSMASH输出目录")
        print("  output_dir: 提取文件的目标目录")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)

    print("=" * 50)
    print("开始提取 region gbk 文件...")
    print("=" * 50)

    stats = extract_region_gbks(input_dir, output_dir)

    print("\n" + "=" * 50)
    print("提取完成！")
    print("=" * 50)
    print(f"总样本数: {stats['total_samples']}")
    print(f"提取region数: {stats['total_regions']}")
    print(f"有BGC的样本: {len(stats['samples_with_regions'])}")
    print(f"无BGC的样本: {len(stats['samples_without_regions'])}")
    print(f"\n输出目录: {output_dir}")

    if stats['samples_without_regions']:
        print("\n以下样本未检测到BGC:")
        for sample in stats['samples_without_regions']:
            print(f"  - {sample}")


if __name__ == "__main__":
    main()
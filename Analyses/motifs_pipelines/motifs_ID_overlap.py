import sys

def main(file1_path, file2_path, output_path):
    # 读取两个motif文件
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        motifs1 = set(file1.read().splitlines())
        motifs2 = set(file2.read().splitlines())

    # 计算交集
    intersection = motifs1 & motifs2

    # 计算不在交集内的motifs
    non_intersection = (motifs1 | motifs2) - intersection

    # 将不在交集内的motifs写入新文件
    with open(output_path, 'w') as output_file:
        for motif in non_intersection:
            output_file.write(f"{motif}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1_path> <file2_path> <output_path>")
    else:
        file1_path = sys.argv[1]
        file2_path = sys.argv[2]
        output_path = sys.argv[3]
        main(file1_path, file2_path, output_path)


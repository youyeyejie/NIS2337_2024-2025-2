import os
import re

# 获取当前目录
root_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
current_directories = [root_directory, os.path.join(root_directory, "pdf")]

def convert_chinese_to_arabic(chinese_number):
    chinese_to_arabic_map = {
        '一': 1, '二': 2, '三': 3, '四': 4, '五': 5,
        '六': 6, '七': 7, '八': 8, '九': 9, '十': 10,
        '百': 100
    }
    result = 0
    temp = 0
    for char in chinese_number:
        if char in chinese_to_arabic_map:
            num = chinese_to_arabic_map[char]
            if num == 10 or num == 100:  # 处理“十”或“百”
                if temp == 0:
                    temp = 1
                result += temp * num
                temp = 0
            else:
                temp = num
    result += temp
    return result

# 遍历当前目录下的文件
for current_directory in current_directories:
    for filename in os.listdir(current_directory):
        # 匹配以“第x章”格式开头的文件名
        match = re.match(r"第([一二三四五六七八九十百]+)章(补充)?", filename)
        if match:
            # 提取章节号并将中文数字转换为阿拉伯数字
            chinese_number = match.group(1)
            chapter_number = str(convert_chinese_to_arabic(chinese_number))
            if match.group(2) == "补充":  # 如果匹配到“补充”
                new_filename = re.sub(r"第[一二三四五六七八九十百]+章补充", f"Ch{chapter_number}sup_", filename)
            else:
                new_filename = re.sub(r"第[一二三四五六七八九十百]+章", f"Ch{chapter_number}_", filename)
            # 重命名文件
            os.rename(os.path.join(current_directory, filename), os.path.join(current_directory, new_filename))
            print(f"重命名: {filename} -> {new_filename}")
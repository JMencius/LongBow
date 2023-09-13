import sys
import re
from multiprocessing import Pool


S = re.compile(r"(\d+)S")


def get_NMID(filename : str, t : int) -> tuple:
    NM_sum = 0
    aligned_length_sum = 0
    linecount = 0
    with open(filename, 'r') as f:
        for line in f:
            if (linecount % cpucount) == t:
                # print(linecount)
                if line[0] != '@':
                    m = line.split()
                    aligned_length_sum += (len(m[9]) - sum([int(i) for i in S.findall(m[5])]))
                    nm_score = int(m[11].split(':')[-1])
                    NM_sum += nm_score
            linecount += 1

    return (NM_sum, aligned_length_sum)

if __name__ == "__main__":
    samfile = str(sys.argv[1])
    cpucount = 11
    
    with Pool(cpucount) as p:
        output = p.starmap(get_NMID, [(samfile, i) for i in range(cpucount)])

    total_error = 0
    total_length = 0
    for i in output:
        total_error += i[0]
        total_length += i[1]
    print(total_error / total_length)




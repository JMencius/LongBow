import sys
from multiprocessing import Pool

def get_quality_score(qstring : str) -> float:
    assert type(qstring) == str
    score = 0
    count = 0
    for i in qstring:
        score += ord(i) - 33
        count += 1
    
    return score / count


def get_aqs(filename : str, t : int) -> float:
    quality_string = ""
    max_quality = 0
    linecount = 0
    with open(filename, 'r') as f:
        for line in f:
            if (linecount % cpucount) == t:
                # print(linecount)
                if line[0] != '@':
                    m = line.split()
                    quality_string = m[10]
                    
                    max_quality = max(max_quality, get_quality_score(quality_string))
            linecount += 1
    # print(f"Parallel {t} completed.")

    return max_quality

if __name__ == "__main__":
    samfile = str(sys.argv[1])
    cpucount = 11

    with Pool(cpucount) as p:
        output = p.starmap(get_aqs, [(samfile, i) for i in range(cpucount)])
    
    print(max([i for i in output]))







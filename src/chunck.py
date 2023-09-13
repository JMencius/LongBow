if __name__ == "__main__":
    import sys
    import random
    
    # print(sys.argv)

    input_fastq = str(sys.argv[1])
    num = int(sys.argv[2])
    
    input_dict = dict()
    i = 0
    wrap = list()
    with open(input_fastq, 'r') as g:
        ori = g.readlines()
        while i < (len(ori) - 3):
            wrap.append((ori[i], ori[i + 1], ori[i + 2], ori[i + 3]))
            i += 4
        del ori
        del i

    random.shuffle(wrap)
    wrap = wrap[ : num]
    with open("random.fastq", 'w') as f:
        for i in wrap:
            for j in i:
                f.write(j)

